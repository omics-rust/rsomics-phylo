import csv
import json
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import dendropy
import numpy as np
from dendropy.calculate import treemeasure
from skbio import DNA, DistanceMatrix, Sequence, TreeNode
from skbio.sequence.distance import hamming, jc69, k2p
from skbio.stats.evolve import hommola_cospeciation
from skbio.tree import nj, upgma


BINARY = Path(sys.argv[1])
TRIMAL = Path(sys.argv[2]) if len(sys.argv) > 2 else None


def execute(*arguments):
    result = subprocess.run(
        [BINARY, *map(str, arguments)], capture_output=True, text=True
    )
    if result.returncode:
        raise RuntimeError(result.stderr)


def write_matrix(path, ids, values):
    with path.open("w", newline="") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(["id", *ids])
        for identity, row in zip(ids, values):
            writer.writerow([identity, *[repr(float(value)) for value in row]])


def read_matrix(path):
    rows = list(csv.reader(path.open(), delimiter="\t"))
    return np.asarray([[float(value) for value in row[1:]] for row in rows[1:]])


def close(actual, expected, tolerance=2e-10):
    if math.isnan(expected):
        assert math.isnan(actual), (actual, expected)
    else:
        assert abs(actual - expected) <= tolerance, (actual, expected)


def distance_oracle(directory):
    rng = np.random.default_rng(4871)
    alphabet = np.asarray(list("ACGTN-."))
    sequences = ["".join(rng.choice(alphabet, 160)) for _ in range(10)]
    fasta = directory / "alignment.fa"
    fasta.write_text(
        "".join(f">s{index} description\n{sequence}\n" for index, sequence in enumerate(sequences))
    )
    maximum = 0.0
    for name, function, constructor in [
        ("hamming", hamming, Sequence),
        ("jc69", jc69, DNA),
        ("k2p", k2p, DNA),
    ]:
        output = directory / f"{name}.tsv"
        execute("distance", "-i", fasta, "-o", output, "--metric", name)
        observed = read_matrix(output)
        for row in range(len(sequences)):
            for column in range(row):
                expected = function(constructor(sequences[row]), constructor(sequences[column]))
                close(observed[row, column], expected)
                if not math.isnan(expected):
                    maximum = max(maximum, abs(observed[row, column] - expected))
    return 135, maximum


def inference_oracle(directory):
    rng = np.random.default_rng(8417)
    maximum = 0.0
    for case in range(20):
        count = int(rng.integers(4, 15))
        points = rng.normal(size=(count, 6))
        values = np.sqrt(((points[:, None, :] - points[None, :, :]) ** 2).sum(axis=2))
        ids = [f"t{index}" for index in range(count)]
        matrix = directory / "inference.tsv"
        write_matrix(matrix, ids, values)
        for method in ("nj", "upgma", "wpgma"):
            output = directory / f"{method}.nwk"
            execute("infer", "-i", matrix, "-o", output, "--method", method)
            observed = TreeNode.read([output.read_text()])
            oracle = (
                nj(DistanceMatrix(values, ids), neg_as_zero=True)
                if method == "nj"
                else upgma(DistanceMatrix(values, ids), weighted=method == "wpgma")
            )
            assert observed.compare_rfd(oracle, rooted=method != "nj") == 0
            difference = observed.compare_wrfd(
                oracle, rooted=method != "nj", metric="euclidean"
            )
            maximum = max(maximum, difference)
            assert difference <= 2e-10
    for values in [
        [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
        [[0, 1, 1, 1], [1, 0, 1, 1], [1, 1, 0, 1], [1, 1, 1, 0]],
        [[0, 2, 4, 10], [2, 0, 6, 14], [4, 6, 0, 20], [10, 14, 20, 0]],
    ]:
        ids = [f"t{index}" for index in range(len(values))]
        matrix = directory / "inference.tsv"
        output = directory / "wpgma.nwk"
        write_matrix(matrix, ids, values)
        execute("infer", "-i", matrix, "-o", output, "--method", "wpgma")
        observed = TreeNode.read([output.read_text()])
        oracle = upgma(DistanceMatrix(values, ids), weighted=True)
        assert observed.compare_rfd(oracle, rooted=True) == 0
        difference = observed.compare_wrfd(oracle, rooted=True, metric="euclidean")
        maximum = max(maximum, difference)
        assert difference <= 2e-10
    matrix = directory / "inference.tsv"
    write_matrix(matrix, ["a", "b"], [[0, 2], [2, 0]])
    output = directory / "wpgma.nwk"
    output.write_text("retained output\n")
    rejected = subprocess.run(
        [BINARY, "infer", "-i", matrix, "-o", output, "--method", "wpgma"],
        capture_output=True, text=True,
    )
    assert rejected.returncode != 0
    assert output.read_text() == "retained output\n"
    try:
        upgma(DistanceMatrix([[0, 2], [2, 0]], ["a", "b"]), weighted=True)
    except ValueError:
        pass
    else:
        raise AssertionError("source unexpectedly accepts two taxa for WPGMA")
    return 64, maximum


def tree_oracles(directory):
    first = "((a:1,b:2):1,c:4,((d:4,e:5):2,f:6):1);"
    second = "((a:3,(b:2,c:2):1):3,d:8,(e:5,f:6):2);"
    first_path = directory / "first.nwk"
    second_path = directory / "second.nwk"
    first_path.write_text(first)
    second_path.write_text(second)
    observed = TreeNode.read([first])
    comparison = TreeNode.read([second])
    expected = {
        "rf": observed.compare_rfd(comparison, rooted=False),
        "wrf": observed.compare_wrfd(comparison, rooted=False),
        "kf": observed.compare_wrfd(comparison, rooted=False, metric="euclidean"),
    }
    for metric, value in expected.items():
        output = directory / "comparison.tsv"
        execute(
            "compare",
            "--reference",
            first_path,
            "--comparison",
            second_path,
            "--metric",
            metric,
            "--output",
            output,
        )
        result = next(csv.DictReader(output.open(), delimiter="\t"))
        close(float(result["value"]), value)

    tree_source = "((a:1,b:1):1,(c:1,d:1):1);"
    tree_path = directory / "ultrametric.nwk"
    tree_path.write_text(tree_source)
    summary_path = directory / "summary.tsv"
    measures = "sackin,normalized-sackin,n-bar,colless,normalized-colless,b1,treeness,gamma"
    execute("summarize", "-i", tree_path, "-o", summary_path, "--measure", measures)
    result = {
        row["measure"]: float(row["value"])
        for row in csv.DictReader(summary_path.open(), delimiter="\t")
    }
    namespace = dendropy.TaxonNamespace()
    tree = dendropy.Tree.get(
        data=tree_source, schema="newick", taxon_namespace=namespace
    )
    expected_measures = {
        "sackin": treemeasure.sackin_index(tree, normalize=False),
        "normalized-sackin": treemeasure.sackin_index(tree, normalize=True),
        "n-bar": treemeasure.N_bar(tree),
        "colless": treemeasure.colless_tree_imbalance(tree, normalize=False),
        "normalized-colless": treemeasure.colless_tree_imbalance(tree, normalize="max"),
        "b1": treemeasure.B1(tree),
        "treeness": treemeasure.treeness(tree),
        "gamma": treemeasure.pybus_harvey_gamma(tree),
    }
    for name, value in expected_measures.items():
        close(result[name], value)
    return 11


def hommola_oracle(directory):
    hosts = np.asarray(
        [[0, 3, 8, 8, 9], [3, 0, 7, 7, 8], [8, 7, 0, 6, 7], [8, 7, 6, 0, 3], [9, 8, 7, 3, 0]],
        dtype=float,
    )
    associates = np.asarray(
        [[0, 5, 8, 8, 8], [5, 0, 7, 7, 7], [8, 7, 0, 4, 4], [8, 7, 4, 0, 2], [8, 7, 4, 2, 0]],
        dtype=float,
    )
    interaction = np.asarray(
        [[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 1, 1]],
        dtype=int,
    )
    host_ids = [f"h{index}" for index in range(5)]
    associate_ids = [f"p{index}" for index in range(5)]
    host_path = directory / "hosts.tsv"
    associate_path = directory / "associates.tsv"
    interaction_path = directory / "interactions.tsv"
    write_matrix(host_path, host_ids, hosts)
    write_matrix(associate_path, associate_ids, associates)
    with interaction_path.open("w", newline="") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(["associate", *host_ids])
        for identity, row in zip(associate_ids, interaction):
            writer.writerow([identity, *row])
    output = directory / "hommola.tsv"
    execute(
        "cospeciation",
        "--hosts",
        host_path,
        "--associates",
        associate_path,
        "--interactions",
        interaction_path,
        "--permutations",
        "0",
        "--output",
        output,
    )
    observed = float(next(csv.DictReader(output.open(), delimiter="\t"))["correlation"])
    expected = hommola_cospeciation(hosts, associates, interaction, permutations=0)[0]
    close(observed, expected)
    return abs(observed - expected)


def trimal_oracle(directory):
    if TRIMAL is None:
        return 0
    source = ">a\nACGTACGT\n>b\nAC-TACGT\n>c\nAC--ACGT\n>d\nACG-ACGT\n"
    input_path = directory / "trim.fa"
    rust_path = directory / "rust.fa"
    oracle_path = directory / "trimal.fa"
    input_path.write_text(source)
    execute(
        "trim",
        "-i",
        input_path,
        "-o",
        rust_path,
        "--maximum-gap-fraction",
        "0.25",
    )
    subprocess.run(
        [TRIMAL, "-in", input_path, "-out", oracle_path, "-gt", "0.75", "-fasta"],
        check=True,
        capture_output=True,
        text=True,
    )
    rust_sequences = [line for line in rust_path.read_text().splitlines() if not line.startswith(">")]
    oracle_sequences = [line for line in oracle_path.read_text().splitlines() if not line.startswith(">")]
    assert rust_sequences == oracle_sequences
    return 1


with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as temporary:
    directory = Path(temporary)
    distance_checks, distance_difference = distance_oracle(directory)
    inference_checks, inference_difference = inference_oracle(directory)
    tree_checks = tree_oracles(directory)
    hommola_difference = hommola_oracle(directory)
    trim_checks = trimal_oracle(directory)
    print(
        json.dumps(
            {
                "distance_checks": distance_checks,
                "distance_max_abs": distance_difference,
                "inference_checks": inference_checks,
                "inference_max_branch_score": inference_difference,
                "tree_checks": tree_checks,
                "hommola_abs": hommola_difference,
                "trim_checks": trim_checks,
            },
            sort_keys=True,
        )
    )
