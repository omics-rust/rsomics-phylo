import argparse
import csv
import gc
import json
import math
import time
from pathlib import Path

import dendropy
import numpy as np
from dendropy.calculate import treemeasure
from skbio import DNA, DistanceMatrix, Sequence, TreeNode
from skbio.sequence.distance import hamming, jc69, k2p
from skbio.stats.evolve import hommola_cospeciation
from skbio.tree import nj, upgma


def measure(operation, repetitions):
    elapsed = []
    checksum = None
    for _ in range(repetitions):
        gc.collect()
        started = time.perf_counter()
        checksum = operation()
        elapsed.append(time.perf_counter() - started)
    return {"elapsed_seconds": elapsed, "checksum": checksum}


def read_fasta(path):
    sequences = []
    current = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if current:
                sequences.append("".join(current))
                current = []
        else:
            current.append(line)
    if current:
        sequences.append("".join(current))
    return sequences


def distance(arguments):
    function, constructor = {
        "hamming": (hamming, Sequence),
        "jc69": (jc69, DNA),
        "k2p": (k2p, DNA),
    }[arguments.metric]
    sequences = [constructor(sequence) for sequence in read_fasta(arguments.input)]

    def operation():
        checksum = 0.0
        for row in range(1, len(sequences)):
            for column in range(row):
                value = function(sequences[row], sequences[column])
                if not math.isnan(value):
                    checksum += value
        return checksum

    return measure(operation, arguments.repetitions)


def euclidean_matrix(taxa):
    state = 77
    mask = (1 << 64) - 1
    points = np.empty((taxa, 8), dtype=np.float64)
    for row in range(taxa):
        for column in range(8):
            state = (state * 6364136223846793005 + 1442695040888963407) & mask
            points[row, column] = (state >> 11) / float(1 << 53)
    squares = np.sum(points * points, axis=1)
    values = np.sqrt(
        np.maximum(squares[:, None] + squares[None, :] - 2 * points @ points.T, 0.0)
    )
    np.fill_diagonal(values, 0.0)
    return DistanceMatrix(values, [f"t{index}" for index in range(taxa)], validate=False)


def inference(arguments):
    matrix = euclidean_matrix(arguments.taxa)
    return {
        "nj": measure(lambda: nj(matrix).count(tips=True), arguments.repetitions),
        "upgma": measure(lambda: upgma(matrix).count(tips=True), arguments.repetitions),
    }


def tree(arguments):
    with arguments.reference.open() as handle:
        reference = TreeNode.read(handle, convert_underscores=False)
    with arguments.comparison.open() as handle:
        comparison = TreeNode.read(handle, convert_underscores=False)
    operations = {
        "rf": lambda: reference.compare_rfd(comparison, rooted=False),
        "wrf": lambda: reference.compare_wrfd(
            comparison, metric="cityblock", rooted=False, include_tips=True
        ),
        "kf": lambda: reference.compare_wrfd(
            comparison, metric="euclidean", rooted=False, include_tips=True
        ),
        "cophenetic": lambda: reference.compare_cophenet(
            comparison, metric="unitcorr", use_length=True
        ),
        "tip-distance": lambda: float(reference.cophenet(use_length=True).data.sum()),
    }
    return {
        name: measure(operation, arguments.repetitions)
        for name, operation in operations.items()
    }


def summarize(arguments):
    phylogeny = dendropy.Tree.get(path=arguments.tree, schema="newick")

    def operation():
        values = [
            treemeasure.sackin_index(phylogeny, normalize=False),
            treemeasure.sackin_index(phylogeny, normalize=True),
            treemeasure.N_bar(phylogeny),
            treemeasure.colless_tree_imbalance(phylogeny, normalize=False),
            treemeasure.colless_tree_imbalance(phylogeny, normalize="max"),
            treemeasure.B1(phylogeny),
            treemeasure.treeness(phylogeny),
            treemeasure.pybus_harvey_gamma(phylogeny),
        ]
        return float(sum(values))

    return measure(operation, arguments.repetitions)


def labeled_matrix(path):
    rows = list(csv.reader(path.open(), delimiter="\t"))
    return np.asarray(
        [[float(value) for value in row[1:]] for row in rows[1:]], dtype=float
    )


def interaction_matrix(path):
    rows = list(csv.reader(path.open(), delimiter="\t"))
    return np.asarray(
        [[int(value) for value in row[1:]] for row in rows[1:]], dtype=int
    )


def cospeciation(arguments):
    hosts = labeled_matrix(arguments.hosts)
    associates = labeled_matrix(arguments.associates)
    interactions = interaction_matrix(arguments.interactions)

    def operation():
        result = hommola_cospeciation(
            hosts,
            associates,
            interactions,
            permutations=arguments.permutations,
            seed=7,
        )
        return float(result[0] + result[1])

    return measure(operation, arguments.repetitions)


parser = argparse.ArgumentParser()
subcommands = parser.add_subparsers(dest="command", required=True)

distance_parser = subcommands.add_parser("distance")
distance_parser.add_argument("input", type=Path)
distance_parser.add_argument("metric", choices=["hamming", "jc69", "k2p"])
distance_parser.add_argument("--repetitions", type=int, default=5)

inference_parser = subcommands.add_parser("inference")
inference_parser.add_argument("taxa", type=int)
inference_parser.add_argument("--repetitions", type=int, default=5)

tree_parser = subcommands.add_parser("tree")
tree_parser.add_argument("reference", type=Path)
tree_parser.add_argument("comparison", type=Path)
tree_parser.add_argument("--repetitions", type=int, default=5)

summary_parser = subcommands.add_parser("summarize")
summary_parser.add_argument("tree", type=Path)
summary_parser.add_argument("--repetitions", type=int, default=5)

cospeciation_parser = subcommands.add_parser("cospeciation")
cospeciation_parser.add_argument("hosts", type=Path)
cospeciation_parser.add_argument("associates", type=Path)
cospeciation_parser.add_argument("interactions", type=Path)
cospeciation_parser.add_argument("--permutations", type=int, default=99)
cospeciation_parser.add_argument("--repetitions", type=int, default=5)

arguments = parser.parse_args()
result = {
    "distance": distance,
    "inference": inference,
    "tree": tree,
    "summarize": summarize,
    "cospeciation": cospeciation,
}[arguments.command](arguments)
print(json.dumps(result, sort_keys=True))
