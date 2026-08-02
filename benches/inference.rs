use std::time::Instant;

use rsomics_phylo::{DistanceMatrix, InferenceMethod, NegativeBranchPolicy, infer_tree};

fn main() {
    let smoke = std::env::args().any(|argument| argument == "--test")
        || std::env::var_os("RSOMICS_BENCH_TAXA").is_none();
    let taxa = if smoke {
        40
    } else {
        environment("RSOMICS_BENCH_TAXA", 800)
    };
    let repetitions = if smoke {
        1
    } else {
        environment("RSOMICS_BENCH_REPETITIONS", 5)
    };
    let matrix = euclidean_matrix(taxa);
    for method in [InferenceMethod::NeighborJoining, InferenceMethod::Upgma] {
        let mut elapsed = Vec::with_capacity(repetitions);
        for _ in 0..repetitions {
            let started = Instant::now();
            let tree = infer_tree(&matrix, method, NegativeBranchPolicy::Zero).unwrap();
            elapsed.push(started.elapsed().as_secs_f64());
            assert_eq!(tree.tip_count(), taxa);
        }
        println!("{method:?}\t{taxa}\t{elapsed:?}");
    }
}

fn environment(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(default)
}

fn euclidean_matrix(taxa: usize) -> DistanceMatrix {
    let dimensions = 8;
    let mut state = 77u64;
    let mut points = vec![0.0; taxa * dimensions];
    for value in &mut points {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        *value = (state >> 11) as f64 / (1u64 << 53) as f64;
    }
    let mut values = vec![0.0; taxa * taxa];
    for row in 0..taxa {
        for column in 0..row {
            let distance = (0..dimensions)
                .map(|dimension| {
                    let difference = points[row * dimensions + dimension]
                        - points[column * dimensions + dimension];
                    difference * difference
                })
                .sum::<f64>()
                .sqrt();
            values[row * taxa + column] = distance;
            values[column * taxa + row] = distance;
        }
    }
    DistanceMatrix::new((0..taxa).map(|index| format!("t{index}")).collect(), values).unwrap()
}
