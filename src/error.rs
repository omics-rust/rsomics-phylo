#[derive(Debug, thiserror::Error, PartialEq)]
#[non_exhaustive]
pub enum PhyloError {
    #[error("alignment must contain at least two records")]
    TooFewSequences,
    #[error("alignment sequence width must be greater than zero")]
    EmptySequence,
    #[error("record header has no sequence identity")]
    EmptyIdentity,
    #[error("duplicate identity {0:?}")]
    DuplicateIdentity(String),
    #[error("record {id:?} has width {observed}, expected {expected}")]
    RaggedAlignment {
        id: String,
        observed: usize,
        expected: usize,
    },
    #[error("record {id:?} contains invalid alignment byte 0x{byte:02x}")]
    InvalidAlignmentByte { id: String, byte: u8 },
    #[error("gap fraction threshold must be finite and within [0, 1], got {0}")]
    InvalidGapThreshold(f64),
    #[error("trimming removed every alignment column")]
    EmptyTrimResult,
    #[error("distance matrix must contain at least two identities")]
    TooFewMatrixRows,
    #[error("distance matrix shape is {observed}, expected {expected} cells")]
    MatrixShape { observed: usize, expected: usize },
    #[error("distance matrix diagonal at row {0} must be zero")]
    NonzeroDiagonal(usize),
    #[error("distance matrix contains infinity at row {row}, column {column}")]
    InfiniteDistance { row: usize, column: usize },
    #[error("distance matrix is asymmetric at rows {left} and {right}")]
    AsymmetricDistance { left: usize, right: usize },
    #[error("tree inference requires a finite distance at rows {row} and {column}")]
    NonFiniteInferenceDistance { row: usize, column: usize },
    #[error(
        "tree inference requires nonnegative distances, got {value} at rows {row} and {column}"
    )]
    NegativeInferenceDistance {
        row: usize,
        column: usize,
        value: f64,
    },
    #[error("invalid delimited table: {0}")]
    InvalidTable(String),
    #[error("tree error: {0}")]
    Tree(String),
    #[error(
        "tree tip sets differ: only reference={reference_only:?}, only comparison={comparison_only:?}"
    )]
    TaxonMismatch {
        reference_only: Vec<String>,
        comparison_only: Vec<String>,
    },
    #[error("node {0} is missing a branch length")]
    MissingBranchLength(usize),
    #[error("tree metric requires nonnegative branch lengths, got {value} at node {node}")]
    NegativeBranchLength { node: usize, value: f64 },
    #[error("metric is undefined: {0}")]
    UndefinedMetric(String),
    #[error("strict bifurcation is required at node {node}, which has {children} children")]
    NotBifurcating { node: usize, children: usize },
    #[error("tree is not ultrametric: tip depths span {minimum} to {maximum}")]
    NotUltrametric { minimum: f64, maximum: f64 },
    #[error("Hommola test requires at least three hosts, three associates, and three interactions")]
    TooFewInteractions,
    #[error(
        "interaction matrix shape is {rows}x{columns}, expected {expected_rows}x{expected_columns}"
    )]
    InteractionShape {
        rows: usize,
        columns: usize,
        expected_rows: usize,
        expected_columns: usize,
    },
    #[error("interaction identities do not match their distance matrices")]
    InteractionIdentityMismatch,
    #[error("permutation count exceeds addressable capacity")]
    TooManyPermutations,
}

pub type Result<T> = std::result::Result<T, PhyloError>;

impl From<rsomics_phylo_tree::BuildError> for PhyloError {
    fn from(error: rsomics_phylo_tree::BuildError) -> Self {
        Self::Tree(error.to_string())
    }
}

impl From<rsomics_phylo_tree::TipIndexError> for PhyloError {
    fn from(error: rsomics_phylo_tree::TipIndexError) -> Self {
        Self::Tree(error.to_string())
    }
}
