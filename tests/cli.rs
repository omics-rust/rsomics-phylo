use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

fn run(arguments: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_rsomics-phylo"))
        .args(arguments)
        .output()
        .unwrap()
}

fn assert_success(output: &Output) {
    assert!(
        output.status.success(),
        "status={}\nstdout={}\nstderr={}",
        output.status,
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn wpgma_weights_clusters_equally_and_preserves_failed_output() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("distance.tsv");
    let output = directory.path().join("tree.nwk");
    fs::write(
        &input,
        "id\ta\tb\tc\td\na\t0\t2\t4\t10\nb\t2\t0\t6\t14\nc\t4\t6\t0\t20\nd\t10\t14\t20\t0\n",
    )
    .unwrap();
    let observed = run(&[
        "infer",
        "-i",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--method",
        "wpgma",
    ]);
    assert_success(&observed);
    let weighted = fs::read_to_string(&output).unwrap();
    assert_eq!(weighted.trim(), "(d:8,(c:2.5,(a:1,b:1):1.5):5.5);");
    let observed = run(&[
        "infer",
        "-i",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--method",
        "upgma",
    ]);
    assert_success(&observed);
    assert_ne!(fs::read_to_string(&output).unwrap(), weighted);
    fs::write(&output, b"retained output\n").unwrap();
    let observed = run(&[
        "infer",
        "-i",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--method",
        "wpgma",
        "--negative-branches",
        "allow",
    ]);
    assert!(!observed.status.success());
    assert_eq!(fs::read(&output).unwrap(), b"retained output\n");
    fs::write(&input, b"id\ta\tb\na\t0\t2\nb\t2\t0\n").unwrap();
    let observed = run(&[
        "infer",
        "-i",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--method",
        "wpgma",
    ]);
    assert!(!observed.status.success());
    assert_eq!(fs::read(&output).unwrap(), b"retained output\n");
    fs::write(&input, b"id\ta\tb\na\t0\t-1\nb\t-1\t0\n").unwrap();
    let observed = run(&[
        "infer",
        "-i",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "--method",
        "wpgma",
    ]);
    assert!(!observed.status.success());
    assert_eq!(fs::read(&output).unwrap(), b"retained output\n");
}

#[test]
fn every_initial_subcommand_executes_real_work() {
    let directory = tempfile::tempdir().unwrap();
    let trim = directory.path().join("trim.fa");
    let map = directory.path().join("map.tsv");
    let output = run(&[
        "trim",
        "--input",
        fixture("aligned.fasta").to_str().unwrap(),
        "--output",
        trim.to_str().unwrap(),
        "--maximum-gap-fraction",
        "0",
        "--column-map",
        map.to_str().unwrap(),
    ]);
    assert_success(&output);
    assert!(
        fs::read_to_string(&trim)
            .unwrap()
            .contains(">a alpha\nACTACGT\n")
    );
    assert_eq!(fs::read_to_string(&map).unwrap().lines().count(), 8);

    let distances = directory.path().join("distance.tsv");
    let output = run(&[
        "distance",
        "--input",
        fixture("aligned.fasta").to_str().unwrap(),
        "--output",
        distances.to_str().unwrap(),
        "--metric",
        "hamming",
    ]);
    assert_success(&output);
    assert!(
        fs::read_to_string(&distances)
            .unwrap()
            .starts_with("id\ta\tb\tc\td\n")
    );

    let tree = directory.path().join("tree.nwk");
    let output = run(&[
        "infer",
        "--input",
        fixture("distance.tsv").to_str().unwrap(),
        "--output",
        tree.to_str().unwrap(),
        "--method",
        "nj",
    ]);
    assert_success(&output);
    assert!(fs::read_to_string(&tree).unwrap().trim_end().ends_with(';'));

    let comparison = directory.path().join("comparison.tsv");
    let output = run(&[
        "compare",
        "--reference",
        fixture("tree1.nwk").to_str().unwrap(),
        "--comparison",
        fixture("tree2.nwk").to_str().unwrap(),
        "--metric",
        "wrf",
        "--output",
        comparison.to_str().unwrap(),
    ]);
    assert_success(&output);
    assert!(
        fs::read_to_string(&comparison)
            .unwrap()
            .contains("wrf\t16\t")
    );

    let tips = directory.path().join("tips.tsv");
    let output = run(&[
        "tip-distance",
        "--input",
        fixture("tree1.nwk").to_str().unwrap(),
        "--output",
        tips.to_str().unwrap(),
    ]);
    assert_success(&output);
    assert!(
        fs::read_to_string(&tips)
            .unwrap()
            .starts_with("id\ta\tb\tc\td\te\tf\n")
    );

    let summary = directory.path().join("summary.tsv");
    let mut child = Command::new(env!("CARGO_BIN_EXE_rsomics-phylo"))
        .args([
            "summarize",
            "--input",
            "-",
            "--output",
            summary.to_str().unwrap(),
            "--measure",
            "sackin,b1",
        ])
        .stdin(std::process::Stdio::piped())
        .spawn()
        .unwrap();
    use std::io::Write;
    child
        .stdin
        .take()
        .unwrap()
        .write_all(b"((a:1,b:1):1,(c:1,d:1):1);")
        .unwrap();
    assert!(child.wait().unwrap().success());
    assert!(
        fs::read_to_string(&summary)
            .unwrap()
            .contains("sackin\t\t8")
    );

    let hommola = directory.path().join("hommola.tsv");
    let null = directory.path().join("null.tsv");
    let output = run(&[
        "cospeciation",
        "--hosts",
        fixture("host.tsv").to_str().unwrap(),
        "--associates",
        fixture("associate.tsv").to_str().unwrap(),
        "--interactions",
        fixture("interactions.tsv").to_str().unwrap(),
        "--permutations",
        "10",
        "--seed",
        "7",
        "--output",
        hommola.to_str().unwrap(),
        "--null-output",
        null.to_str().unwrap(),
    ]);
    assert_success(&output);
    assert_eq!(fs::read_to_string(&null).unwrap().lines().count(), 11);
}

#[test]
fn recursive_help_and_json_contract_are_available() {
    for command in [
        "trim",
        "distance",
        "infer",
        "compare",
        "tip-distance",
        "summarize",
        "cospeciation",
    ] {
        let output = run(&[command, "--help"]);
        assert_success(&output);
        assert!(!output.stdout.is_empty());
    }

    let directory = tempfile::tempdir().unwrap();
    let matrix = directory.path().join("distance.tsv");
    let output = run(&[
        "distance",
        "--input",
        fixture("aligned.fasta").to_str().unwrap(),
        "--output",
        matrix.to_str().unwrap(),
        "--json",
    ]);
    assert_success(&output);
    let json = String::from_utf8(output.stdout).unwrap();
    assert!(json.contains("\"schema_version\":\"1.0\""));
    assert!(json.contains("\"operation\":\"distance\""));
}

#[test]
fn incompatible_options_fail_before_output() {
    let directory = tempfile::tempdir().unwrap();
    let result = directory.path().join("result.tsv");
    let output = run(&[
        "compare",
        "--reference",
        fixture("tree1.nwk").to_str().unwrap(),
        "--comparison",
        fixture("tree2.nwk").to_str().unwrap(),
        "--metric",
        "wrf",
        "--proportion",
        "--output",
        result.to_str().unwrap(),
    ]);
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("applies only"));

    let output = run(&[
        "trim",
        "--input",
        fixture("aligned.fasta").to_str().unwrap(),
        "--column-map",
        "-",
    ]);
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("named file"));

    let output = run(&[
        "infer",
        "--input",
        fixture("distance.tsv").to_str().unwrap(),
        "--threads",
        "2",
    ]);
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("applies only"));
}
