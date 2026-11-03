use std::fs::File;
use std::io::{self, Cursor, Read, Write};
use std::path::{Path, PathBuf};

use clap::{Args, Parser, Subcommand, ValueEnum};
use rsomics_common::{
    Context, OutputArgs, Result, RsomicsError, ThreadArgs, ToolMeta, reject_output_alias,
    write_atomic_pair, write_output,
};
use rsomics_phylo::{
    AlignedRecord, Alignment, CompareMetric, CompareOptions, CopheneticMetric, DistanceMatrix,
    HommolaOptions, InferenceMethod, InteractionMatrix, NegativeBranchPolicy, NodeAge,
    SequenceMetric, TreeMeasure, compare_trees, hommola_cospeciation, infer_tree, node_ages,
    sequence_distances, tip_distances, tree_measure,
};
use rsomics_phylo_tree::{ParseOptions, Tree};
use serde::Serialize;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Debug, Parser)]
#[command(
    name = "rsomics-phylo",
    version,
    about = "Run coherent distance-based phylogenetic workflows",
    arg_required_else_help = true
)]
pub struct Cli {
    #[command(subcommand)]
    command: Command,

    #[command(flatten)]
    threads: ThreadArgs,

    #[command(flatten)]
    pub output: OutputArgs,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Remove alignment columns exceeding a gap-fraction threshold.
    Trim(TrimArgs),
    /// Calculate a labeled sequence-distance matrix.
    Distance(DistanceArgs),
    /// Infer a distance-based phylogenetic tree.
    Infer(InferArgs),
    /// Compare two trees using topology, branch, or path metrics.
    Compare(CompareArgs),
    /// Calculate all named-tip path distances.
    TipDistance(TipDistanceArgs),
    /// Calculate selected tree timing and balance measures.
    Summarize(SummarizeArgs),
    /// Run the Hommola host-associate cospeciation test.
    Cospeciation(CospeciationArgs),
}

impl Command {
    fn writes_data_to_stdout(&self) -> bool {
        let output = match self {
            Self::Trim(args) => args.io.output.as_deref(),
            Self::Distance(args) => args.io.output.as_deref(),
            Self::Infer(args) => args.io.output.as_deref(),
            Self::Compare(args) => args.output.as_deref(),
            Self::TipDistance(args) => args.io.output.as_deref(),
            Self::Summarize(args) => args.io.output.as_deref(),
            Self::Cospeciation(args) => args.output.as_deref(),
        };
        output.is_none() || output == Some(Path::new("-"))
    }

    fn supports_threads(&self) -> bool {
        matches!(self, Self::Distance(_))
    }
}

#[derive(Debug, Args)]
struct InputOutput {
    /// Input path; use - for standard input.
    #[arg(short, long)]
    input: PathBuf,

    /// Transactional output path; omit or use - for standard output.
    #[arg(short, long)]
    output: Option<PathBuf>,
}

#[derive(Debug, Args)]
struct TrimArgs {
    #[command(flatten)]
    io: InputOutput,

    /// Maximum retained-column gap fraction, including '-' and '.'.
    #[arg(long, default_value_t = 0.5)]
    maximum_gap_fraction: f64,

    /// Optional transactional zero-based retained-column map.
    #[arg(long, value_name = "TSV")]
    column_map: Option<PathBuf>,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum SequenceMetricArg {
    Hamming,
    Jc69,
    K2p,
}

impl From<SequenceMetricArg> for SequenceMetric {
    fn from(value: SequenceMetricArg) -> Self {
        match value {
            SequenceMetricArg::Hamming => Self::Hamming,
            SequenceMetricArg::Jc69 => Self::Jc69,
            SequenceMetricArg::K2p => Self::K2p,
        }
    }
}

#[derive(Debug, Args)]
struct DistanceArgs {
    #[command(flatten)]
    io: InputOutput,

    /// Sequence-distance model.
    #[arg(long, value_enum, default_value_t = SequenceMetricArg::Hamming)]
    metric: SequenceMetricArg,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum InferenceMethodArg {
    Nj,
    Upgma,
    Wpgma,
}

impl From<InferenceMethodArg> for InferenceMethod {
    fn from(value: InferenceMethodArg) -> Self {
        match value {
            InferenceMethodArg::Nj => Self::NeighborJoining,
            InferenceMethodArg::Upgma => Self::Upgma,
            InferenceMethodArg::Wpgma => Self::Wpgma,
        }
    }
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum NegativeBranchArg {
    Zero,
    Allow,
}

impl From<NegativeBranchArg> for NegativeBranchPolicy {
    fn from(value: NegativeBranchArg) -> Self {
        match value {
            NegativeBranchArg::Zero => Self::Zero,
            NegativeBranchArg::Allow => Self::Allow,
        }
    }
}

#[derive(Debug, Args)]
struct InferArgs {
    #[command(flatten)]
    io: InputOutput,

    /// Tree inference algorithm.
    #[arg(long, value_enum, default_value_t = InferenceMethodArg::Nj)]
    method: InferenceMethodArg,

    /// Neighbor-joining policy for negative inferred branches.
    #[arg(long, value_enum, default_value_t = NegativeBranchArg::Zero)]
    negative_branches: NegativeBranchArg,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum CompareMetricArg {
    Rf,
    Wrf,
    Kf,
    FpFn,
    Cophenetic,
}

impl From<CompareMetricArg> for CompareMetric {
    fn from(value: CompareMetricArg) -> Self {
        match value {
            CompareMetricArg::Rf => Self::RobinsonFoulds,
            CompareMetricArg::Wrf => Self::WeightedRobinsonFoulds,
            CompareMetricArg::Kf => Self::BranchScore,
            CompareMetricArg::FpFn => Self::FalsePositiveNegative,
            CompareMetricArg::Cophenetic => Self::Cophenetic,
        }
    }
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum CopheneticMetricArg {
    UnitCorrelation,
    Euclidean,
    CityBlock,
}

impl From<CopheneticMetricArg> for CopheneticMetric {
    fn from(value: CopheneticMetricArg) -> Self {
        match value {
            CopheneticMetricArg::UnitCorrelation => Self::UnitCorrelation,
            CopheneticMetricArg::Euclidean => Self::Euclidean,
            CopheneticMetricArg::CityBlock => Self::CityBlock,
        }
    }
}

#[derive(Debug, Args)]
struct CompareArgs {
    /// Reference Newick tree; use - for standard input.
    #[arg(long)]
    reference: PathBuf,

    /// Comparison Newick tree.
    #[arg(long)]
    comparison: PathBuf,

    /// Tree comparison metric.
    #[arg(long, value_enum, default_value_t = CompareMetricArg::Rf)]
    metric: CompareMetricArg,

    /// Interpret splits as rooted clades.
    #[arg(long)]
    rooted: bool,

    /// Exclude terminal branches from weighted metrics.
    #[arg(long)]
    exclude_tips: bool,

    /// Normalize RF count by the total number of splits.
    #[arg(long)]
    proportion: bool,

    /// Vector distance used by cophenetic comparison.
    #[arg(long, value_enum, default_value_t = CopheneticMetricArg::UnitCorrelation)]
    cophenetic_metric: CopheneticMetricArg,

    /// Convert unquoted Newick underscores to spaces.
    #[arg(long)]
    convert_underscores: bool,

    /// Transactional comparison result; omit or use - for standard output.
    #[arg(short, long)]
    output: Option<PathBuf>,
}

#[derive(Debug, Args)]
struct TipDistanceArgs {
    #[command(flatten)]
    io: InputOutput,

    /// Count edges instead of summing branch lengths.
    #[arg(long)]
    edge_count: bool,

    /// Convert unquoted Newick underscores to spaces.
    #[arg(long)]
    convert_underscores: bool,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum MeasureArg {
    Sackin,
    NormalizedSackin,
    NBar,
    Colless,
    NormalizedColless,
    B1,
    Treeness,
    Gamma,
    NodeAges,
    CoalescenceAges,
}

#[derive(Debug, Args)]
struct SummarizeArgs {
    #[command(flatten)]
    io: InputOutput,

    /// Measures to calculate, as a comma-separated list.
    #[arg(long, value_enum, value_delimiter = ',', required = true)]
    measure: Vec<MeasureArg>,

    /// Maximum allowed root-to-tip depth span for timing measures.
    #[arg(long, default_value_t = 1e-5)]
    ultrametric_tolerance: f64,

    /// Convert unquoted Newick underscores to spaces.
    #[arg(long)]
    convert_underscores: bool,
}

#[derive(Debug, Args)]
struct CospeciationArgs {
    /// Labeled host distance matrix.
    #[arg(long)]
    hosts: PathBuf,

    /// Labeled associate distance matrix.
    #[arg(long)]
    associates: PathBuf,

    /// Associate-by-host binary interaction matrix.
    #[arg(long)]
    interactions: PathBuf,

    /// Number of one-sided permutations; zero calculates only the statistic.
    #[arg(long, default_value_t = 999)]
    permutations: usize,

    /// SplitMix64 seed for the native deterministic stream.
    #[arg(long, default_value_t = 0)]
    seed: u64,

    /// Optional transactional null-distribution table.
    #[arg(long, value_name = "TSV")]
    null_output: Option<PathBuf>,

    /// Transactional result table; omit or use - for standard output.
    #[arg(short, long)]
    output: Option<PathBuf>,
}

#[derive(Debug, Serialize)]
pub struct ExecutionReport {
    operation: &'static str,
    items: usize,
    output: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    detail: Option<String>,
}

impl Cli {
    pub fn execute(self) -> Result<ExecutionReport> {
        if self.output.json && self.command.writes_data_to_stdout() {
            return Err(RsomicsError::ConfigError(
                "--json requires a named data output".into(),
            ));
        }
        let requested_threads = self.threads.requested();
        if requested_threads.is_some() && !self.command.supports_threads() {
            return Err(RsomicsError::ConfigError(
                "--threads applies only to distance".into(),
            ));
        }
        if let Some(threads) = requested_threads {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads.get())
                .build_global()
                .map_err(|error| RsomicsError::ConfigError(error.to_string()))?;
        }
        match self.command {
            Command::Trim(args) => execute_trim(args),
            Command::Distance(args) => execute_distance(args),
            Command::Infer(args) => execute_infer(args),
            Command::Compare(args) => execute_compare(args),
            Command::TipDistance(args) => execute_tip_distance(args),
            Command::Summarize(args) => execute_summarize(args),
            Command::Cospeciation(args) => execute_cospeciation(args),
        }
    }
}

fn execute_trim(args: TrimArgs) -> Result<ExecutionReport> {
    reject_secondary_stdout(args.column_map.as_deref(), "--column-map")?;
    reject_aliases(args.io.output.as_deref(), [args.io.input.as_path()])?;
    reject_aliases(args.column_map.as_deref(), [args.io.input.as_path()])?;
    let alignment = read_alignment(&args.io.input)?;
    let original_width = alignment.width();
    let trimmed = alignment
        .trim_by_gap(args.maximum_gap_fraction)
        .map_err(phylo_error)?;
    let retained = trimmed.kept_columns().len();
    match (args.io.output.as_deref(), args.column_map.as_deref()) {
        (Some(output), Some(map)) if output != Path::new("-") && map != Path::new("-") => {
            write_atomic_pair(output, map, |output, map| {
                write_alignment(trimmed.alignment(), output)?;
                write_column_map(trimmed.kept_columns(), map)
            })?;
        }
        (output, map) => {
            write_output(output, |writer| {
                write_alignment(trimmed.alignment(), writer)
            })?;
            if let Some(map) = map {
                write_output(Some(map), |writer| {
                    write_column_map(trimmed.kept_columns(), writer)
                })?;
            }
        }
    }
    Ok(ExecutionReport {
        operation: "trim",
        items: retained,
        output: output_label(args.io.output.as_deref()),
        detail: Some(format!("retained {retained} of {original_width} columns")),
    })
}

fn execute_distance(args: DistanceArgs) -> Result<ExecutionReport> {
    reject_aliases(args.io.output.as_deref(), [args.io.input.as_path()])?;
    let alignment = read_alignment(&args.io.input)?;
    let matrix = sequence_distances(&alignment, args.metric.into()).map_err(phylo_error)?;
    write_matrix(&matrix, args.io.output.as_deref())?;
    Ok(ExecutionReport {
        operation: "distance",
        items: matrix.size(),
        output: output_label(args.io.output.as_deref()),
        detail: Some(format!("metric={:?}", args.metric)),
    })
}

fn execute_infer(args: InferArgs) -> Result<ExecutionReport> {
    if !matches!(args.method, InferenceMethodArg::Nj)
        && matches!(args.negative_branches, NegativeBranchArg::Allow)
    {
        return Err(RsomicsError::ConfigError(
            "--negative-branches applies only to neighbor joining".into(),
        ));
    }
    reject_aliases(args.io.output.as_deref(), [args.io.input.as_path()])?;
    let matrix = read_matrix(&args.io.input)?;
    let tree = infer_tree(&matrix, args.method.into(), args.negative_branches.into())
        .map_err(phylo_error)?;
    let newick = tree.to_newick();
    write_text(&newick, args.io.output.as_deref())?;
    Ok(ExecutionReport {
        operation: "infer",
        items: tree.tip_count(),
        output: output_label(args.io.output.as_deref()),
        detail: Some(format!("method={:?}", args.method)),
    })
}

fn execute_compare(args: CompareArgs) -> Result<ExecutionReport> {
    if args.proportion && !matches!(args.metric, CompareMetricArg::Rf) {
        return Err(RsomicsError::ConfigError(
            "--proportion applies only to --metric rf".into(),
        ));
    }
    reject_two_stdin(&args.reference, &args.comparison)?;
    reject_aliases(
        args.output.as_deref(),
        [args.reference.as_path(), args.comparison.as_path()],
    )?;
    let reference = read_tree(&args.reference, args.convert_underscores)?;
    let comparison = read_tree(&args.comparison, args.convert_underscores)?;
    let metric = args.metric.into();
    let result = compare_trees(
        &reference,
        &comparison,
        metric,
        CompareOptions {
            rooted: args.rooted,
            include_tips: !args.exclude_tips,
            proportion: args.proportion,
            cophenetic_metric: args.cophenetic_metric.into(),
        },
    )
    .map_err(phylo_error)?;
    write_output(args.output.as_deref(), |writer| {
        writeln!(
            writer,
            "metric\tvalue\tfalse_positives\tfalse_negatives\ttaxa\trooted\tinclude_tips\tproportion"
        )?;
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            compare_metric_name(metric),
            result.value,
            result
                .false_positives
                .map_or_else(String::new, |value| value.to_string()),
            result
                .false_negatives
                .map_or_else(String::new, |value| value.to_string()),
            result.taxa,
            result.rooted,
            result.include_tips,
            result.proportion,
        )?;
        Ok(())
    })?;
    Ok(ExecutionReport {
        operation: "compare",
        items: result.taxa,
        output: output_label(args.output.as_deref()),
        detail: Some(format!(
            "metric={}, value={}",
            compare_metric_name(metric),
            result.value
        )),
    })
}

fn execute_tip_distance(args: TipDistanceArgs) -> Result<ExecutionReport> {
    reject_aliases(args.io.output.as_deref(), [args.io.input.as_path()])?;
    let tree = read_tree(&args.io.input, args.convert_underscores)?;
    let matrix = tip_distances(&tree, !args.edge_count).map_err(phylo_error)?;
    write_matrix(&matrix, args.io.output.as_deref())?;
    Ok(ExecutionReport {
        operation: "tip-distance",
        items: matrix.size(),
        output: output_label(args.io.output.as_deref()),
        detail: Some(
            if args.edge_count {
                "edge-count"
            } else {
                "branch-length"
            }
            .into(),
        ),
    })
}

fn execute_summarize(args: SummarizeArgs) -> Result<ExecutionReport> {
    reject_aliases(args.io.output.as_deref(), [args.io.input.as_path()])?;
    let tree = read_tree(&args.io.input, args.convert_underscores)?;
    let mut rows = Vec::<(String, String, f64)>::new();
    for measure in args.measure {
        match measure {
            MeasureArg::NodeAges => {
                rows.extend(
                    node_ages(&tree, args.ultrametric_tolerance)
                        .map_err(phylo_error)?
                        .into_iter()
                        .map(|age| node_age_row("node-age", age)),
                );
            }
            MeasureArg::CoalescenceAges => {
                let mut ages = node_ages(&tree, args.ultrametric_tolerance)
                    .map_err(phylo_error)?
                    .into_iter()
                    .filter(NodeAge::is_internal)
                    .map(|age| age.age())
                    .collect::<Vec<_>>();
                ages.sort_by(f64::total_cmp);
                rows.extend(
                    ages.into_iter()
                        .enumerate()
                        .map(|(index, age)| ("coalescence-age".into(), index.to_string(), age)),
                );
            }
            scalar => {
                let (name, measure) = scalar_measure(scalar);
                rows.push((
                    name.into(),
                    String::new(),
                    tree_measure(&tree, measure, args.ultrametric_tolerance)
                        .map_err(phylo_error)?,
                ));
            }
        }
    }
    write_output(args.io.output.as_deref(), |writer| {
        writeln!(writer, "measure\tentity\tvalue")?;
        for (measure, entity, value) in &rows {
            writeln!(writer, "{measure}\t{entity}\t{value}")?;
        }
        Ok(())
    })?;
    Ok(ExecutionReport {
        operation: "summarize",
        items: rows.len(),
        output: output_label(args.io.output.as_deref()),
        detail: Some(format!("tips={}", tree.tip_count())),
    })
}

fn execute_cospeciation(args: CospeciationArgs) -> Result<ExecutionReport> {
    reject_secondary_stdout(args.null_output.as_deref(), "--null-output")?;
    reject_three_stdin(&args.hosts, &args.associates, &args.interactions)?;
    let inputs = [
        args.hosts.as_path(),
        args.associates.as_path(),
        args.interactions.as_path(),
    ];
    reject_aliases(args.output.as_deref(), inputs)?;
    reject_aliases(args.null_output.as_deref(), inputs)?;
    let hosts = read_matrix(&args.hosts)?;
    let associates = read_matrix(&args.associates)?;
    let interactions = read_interactions(&args.interactions)?;
    let result = hommola_cospeciation(
        &hosts,
        &associates,
        &interactions,
        HommolaOptions {
            permutations: args.permutations,
            seed: args.seed,
            store_null: args.null_output.is_some(),
        },
    )
    .map_err(phylo_error)?;
    match (args.output.as_deref(), args.null_output.as_deref()) {
        (Some(output), Some(null)) if output != Path::new("-") && null != Path::new("-") => {
            write_atomic_pair(output, null, |output, null| {
                write_hommola(&result, output)?;
                write_null(&result.null_distribution, null)
            })?;
        }
        (output, null) => {
            write_output(output, |writer| write_hommola(&result, writer))?;
            if let Some(null) = null {
                write_output(Some(null), |writer| {
                    write_null(&result.null_distribution, writer)
                })?;
            }
        }
    }
    Ok(ExecutionReport {
        operation: "cospeciation",
        items: result.permutations,
        output: output_label(args.output.as_deref()),
        detail: Some(format!("correlation={}", result.correlation)),
    })
}

fn read_alignment(path: &Path) -> Result<Alignment> {
    let bytes = read_bytes(path)?;
    let mut reader = rsomics_seqio::open_reader(Cursor::new(bytes))?;
    if reader.format() != rsomics_seqio::Format::Fasta {
        return Err(RsomicsError::InvalidInput(
            "aligned sequence input must be FASTA".into(),
        ));
    }
    let mut records = Vec::new();
    while let Some(record) = reader.read_record()? {
        let header = std::str::from_utf8(record.id)?.to_owned();
        records.push(AlignedRecord::new(header, record.seq.to_vec()).map_err(phylo_error)?);
    }
    Alignment::new(records).map_err(phylo_error)
}

fn write_alignment(alignment: &Alignment, writer: &mut dyn Write) -> Result<()> {
    let mut writer = rsomics_seqio::Writer::new(writer, rsomics_seqio::Format::Fasta);
    for record in alignment.records() {
        writer.write_record(rsomics_seqio::Record {
            id: record.header().as_bytes(),
            seq: record.sequence(),
            qual: None,
        })?;
    }
    writer.finish()
}

fn write_column_map(columns: &[usize], writer: &mut dyn Write) -> Result<()> {
    writeln!(writer, "retained_column\toriginal_column")?;
    for (retained, original) in columns.iter().enumerate() {
        writeln!(writer, "{retained}\t{original}")?;
    }
    Ok(())
}

fn read_matrix(path: &Path) -> Result<DistanceMatrix> {
    DistanceMatrix::from_tsv(read_bytes(path)?.as_slice()).map_err(phylo_error)
}

fn write_matrix(matrix: &DistanceMatrix, output: Option<&Path>) -> Result<()> {
    write_output(output, |writer| {
        matrix.write_tsv(writer).map_err(phylo_error)
    })
}

fn read_tree(path: &Path, convert_underscores: bool) -> Result<Tree> {
    let bytes = read_bytes(path)?;
    let source = std::str::from_utf8(&bytes)?;
    Tree::from_newick_with_options(
        source,
        ParseOptions {
            convert_underscores,
        },
    )
    .map_err(|error| RsomicsError::InvalidInput(error.to_string()))
}

fn read_interactions(path: &Path) -> Result<InteractionMatrix> {
    let bytes = read_bytes(path)?;
    let mut input = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .flexible(false)
        .from_reader(bytes.as_slice());
    let header = input
        .headers()
        .map_err(csv_error)?
        .iter()
        .map(str::to_owned)
        .collect::<Vec<_>>();
    if header.len() < 4 || header[0] != "associate" {
        return Err(RsomicsError::InvalidInput(
            "interaction header must be associate followed by at least three host identities"
                .into(),
        ));
    }
    let host_ids = header[1..].to_vec();
    let mut associate_ids = Vec::new();
    let mut values = Vec::new();
    for record in input.records() {
        let record = record.map_err(csv_error)?;
        associate_ids.push(record[0].to_owned());
        for token in record.iter().skip(1) {
            values.push(match token {
                "0" => false,
                "1" => true,
                _ => {
                    return Err(RsomicsError::InvalidInput(format!(
                        "interaction values must be 0 or 1, got {token:?}"
                    )));
                }
            });
        }
    }
    InteractionMatrix::new(associate_ids, host_ids, values).map_err(phylo_error)
}

fn read_bytes(path: &Path) -> Result<Vec<u8>> {
    let mut bytes = Vec::new();
    if path == Path::new("-") {
        io::stdin()
            .lock()
            .read_to_end(&mut bytes)
            .rs_context("reading standard input")?;
    } else {
        File::open(path)
            .rs_with_context(|| format!("opening {}", path.display()))?
            .read_to_end(&mut bytes)
            .rs_with_context(|| format!("reading {}", path.display()))?;
    }
    Ok(bytes)
}

fn write_text(text: &str, output: Option<&Path>) -> Result<()> {
    write_output(output, |writer| {
        writer.write_all(text.as_bytes())?;
        writer.write_all(b"\n")?;
        Ok(())
    })
}

fn write_hommola(result: &rsomics_phylo::HommolaResult, writer: &mut dyn Write) -> Result<()> {
    writeln!(writer, "correlation\tp_value\tpermutations\tseed\trng")?;
    writeln!(
        writer,
        "{}\t{}\t{}\t{}\tsplitmix64",
        result.correlation,
        result
            .p_value
            .map_or_else(String::new, |value| value.to_string()),
        result.permutations,
        result.seed,
    )?;
    Ok(())
}

fn write_null(values: &[f64], writer: &mut dyn Write) -> Result<()> {
    writeln!(writer, "permutation\tcorrelation")?;
    for (index, value) in values.iter().enumerate() {
        writeln!(writer, "{}\t{}", index + 1, value)?;
    }
    Ok(())
}

fn scalar_measure(measure: MeasureArg) -> (&'static str, TreeMeasure) {
    match measure {
        MeasureArg::Sackin => ("sackin", TreeMeasure::Sackin),
        MeasureArg::NormalizedSackin => ("normalized-sackin", TreeMeasure::NormalizedSackin),
        MeasureArg::NBar => ("n-bar", TreeMeasure::NBar),
        MeasureArg::Colless => ("colless", TreeMeasure::Colless),
        MeasureArg::NormalizedColless => ("normalized-colless", TreeMeasure::NormalizedColless),
        MeasureArg::B1 => ("b1", TreeMeasure::B1),
        MeasureArg::Treeness => ("treeness", TreeMeasure::Treeness),
        MeasureArg::Gamma => ("gamma", TreeMeasure::Gamma),
        MeasureArg::NodeAges | MeasureArg::CoalescenceAges => {
            unreachable!("timing measures are handled separately")
        }
    }
}

fn node_age_row(measure: &str, age: NodeAge) -> (String, String, f64) {
    let entity = age
        .name()
        .map_or_else(|| format!("node:{}", age.node()), str::to_owned);
    (measure.into(), entity, age.age())
}

fn compare_metric_name(metric: CompareMetric) -> &'static str {
    match metric {
        CompareMetric::RobinsonFoulds => "rf",
        CompareMetric::WeightedRobinsonFoulds => "wrf",
        CompareMetric::BranchScore => "kf",
        CompareMetric::FalsePositiveNegative => "fp-fn",
        CompareMetric::Cophenetic => "cophenetic",
    }
}

fn reject_aliases<'a>(
    output: Option<&Path>,
    inputs: impl IntoIterator<Item = &'a Path>,
) -> Result<()> {
    if let Some(output) = output {
        reject_output_alias(output, inputs)?;
    }
    Ok(())
}

fn reject_two_stdin(first: &Path, second: &Path) -> Result<()> {
    if first == Path::new("-") && second == Path::new("-") {
        return Err(RsomicsError::ConfigError(
            "two inputs cannot both read standard input".into(),
        ));
    }
    Ok(())
}

fn reject_three_stdin(first: &Path, second: &Path, third: &Path) -> Result<()> {
    let count = [first, second, third]
        .into_iter()
        .filter(|path| *path == Path::new("-"))
        .count();
    if count > 1 {
        return Err(RsomicsError::ConfigError(
            "at most one input may read standard input".into(),
        ));
    }
    Ok(())
}

fn reject_secondary_stdout(path: Option<&Path>, option: &str) -> Result<()> {
    if path == Some(Path::new("-")) {
        return Err(RsomicsError::ConfigError(format!(
            "{option} must be a named file so data streams cannot be mixed"
        )));
    }
    Ok(())
}

fn output_label(path: Option<&Path>) -> String {
    path.filter(|path| *path != Path::new("-"))
        .map_or_else(|| "stdout".into(), |path| path.display().to_string())
}

fn phylo_error(error: rsomics_phylo::PhyloError) -> RsomicsError {
    RsomicsError::InvalidInput(error.to_string())
}

fn csv_error(error: csv::Error) -> RsomicsError {
    RsomicsError::InvalidInput(error.to_string())
}
