#![forbid(unsafe_code)]

mod alignment;
mod comparison;
mod cospeciation;
mod distance;
mod error;
mod inference;
mod matrix;
mod measures;

pub use alignment::{AlignedRecord, Alignment, TrimmedAlignment};
pub use comparison::{
    CompareMetric, CompareOptions, ComparisonResult, CopheneticMetric, compare_trees,
};
pub use cospeciation::{HommolaOptions, HommolaResult, InteractionMatrix, hommola_cospeciation};
pub use distance::{SequenceMetric, sequence_distances};
pub use error::{PhyloError, Result};
pub use inference::{InferenceMethod, NegativeBranchPolicy, infer_tree};
pub use matrix::DistanceMatrix;
pub use measures::{
    NodeAge, TreeMeasure, TreeSummary, node_ages, summarize_tree, tip_distances, tree_measure,
};
