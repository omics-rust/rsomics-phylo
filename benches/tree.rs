use std::fs;
use std::time::Instant;

use rsomics_phylo::{
    CompareMetric, CompareOptions, CopheneticMetric, compare_trees, tip_distances,
};
use rsomics_phylo_tree::Tree;

fn main() {
    let smoke = std::env::args().any(|argument| argument == "--test")
        || std::env::var_os("RSOMICS_BENCH_TREE_A").is_none()
        || std::env::var_os("RSOMICS_BENCH_TREE_B").is_none();
    let (left, right) = if smoke {
        (
            Tree::from_newick("((a:1,b:2):1,(c:2,d:1):2);").unwrap(),
            Tree::from_newick("((a:2,c:1):1,(b:1,d:2):2);").unwrap(),
        )
    } else {
        (
            read_tree("RSOMICS_BENCH_TREE_A"),
            read_tree("RSOMICS_BENCH_TREE_B"),
        )
    };
    let repetitions = if smoke {
        1
    } else {
        environment("RSOMICS_BENCH_REPETITIONS", 5)
    };
    for metric in [
        CompareMetric::RobinsonFoulds,
        CompareMetric::WeightedRobinsonFoulds,
        CompareMetric::BranchScore,
        CompareMetric::Cophenetic,
    ] {
        let options = CompareOptions {
            cophenetic_metric: CopheneticMetric::UnitCorrelation,
            ..CompareOptions::default()
        };
        let mut elapsed = Vec::with_capacity(repetitions);
        for _ in 0..repetitions {
            let started = Instant::now();
            let result = compare_trees(&left, &right, metric, options).unwrap();
            elapsed.push(started.elapsed().as_secs_f64());
            assert!(result.value.is_finite());
        }
        println!("{metric:?}\t{}\t{elapsed:?}", left.tip_count());
    }
    let mut elapsed = Vec::with_capacity(repetitions);
    for _ in 0..repetitions {
        let started = Instant::now();
        let result = tip_distances(&left, true).unwrap();
        elapsed.push(started.elapsed().as_secs_f64());
        assert_eq!(result.size(), left.tip_count());
    }
    println!("TipDistance\t{}\t{elapsed:?}", left.tip_count());
}

fn read_tree(variable: &str) -> Tree {
    let path = std::env::var(variable).unwrap_or_else(|_| panic!("{variable} is required"));
    Tree::from_newick(&fs::read_to_string(path).unwrap()).unwrap()
}

fn environment(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(default)
}
