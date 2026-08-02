use std::collections::{HashMap, HashSet};

use rsomics_phylo_tree::Tree;

use crate::{PhyloError, Result, tip_distances};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CompareMetric {
    RobinsonFoulds,
    WeightedRobinsonFoulds,
    BranchScore,
    FalsePositiveNegative,
    Cophenetic,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CopheneticMetric {
    UnitCorrelation,
    Euclidean,
    CityBlock,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CompareOptions {
    pub rooted: bool,
    pub include_tips: bool,
    pub proportion: bool,
    pub cophenetic_metric: CopheneticMetric,
}

impl Default for CompareOptions {
    fn default() -> Self {
        Self {
            rooted: false,
            include_tips: true,
            proportion: false,
            cophenetic_metric: CopheneticMetric::UnitCorrelation,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct ComparisonResult {
    pub metric: CompareMetric,
    pub value: f64,
    pub false_positives: Option<usize>,
    pub false_negatives: Option<usize>,
    pub taxa: usize,
    pub rooted: bool,
    pub include_tips: bool,
    pub proportion: bool,
}

pub fn compare_trees(
    reference: &Tree,
    comparison: &Tree,
    metric: CompareMetric,
    options: CompareOptions,
) -> Result<ComparisonResult> {
    let taxa = exact_taxa(reference, comparison)?;
    let (value, false_positives, false_negatives) = match metric {
        CompareMetric::RobinsonFoulds => {
            let reference_splits = topology_splits(reference, &taxa, options.rooted)?;
            let comparison_splits = topology_splits(comparison, &taxa, options.rooted)?;
            let count = reference_splits
                .symmetric_difference(&comparison_splits)
                .count();
            let value = if options.proportion {
                let denominator = reference_splits.len() + comparison_splits.len();
                if denominator == 0 {
                    0.0
                } else {
                    count as f64 / denominator as f64
                }
            } else {
                count as f64
            };
            (value, None, None)
        }
        CompareMetric::FalsePositiveNegative => {
            let reference_splits = topology_splits(reference, &taxa, options.rooted)?;
            let comparison_splits = topology_splits(comparison, &taxa, options.rooted)?;
            let false_positives = comparison_splits.difference(&reference_splits).count();
            let false_negatives = reference_splits.difference(&comparison_splits).count();
            (
                (false_positives + false_negatives) as f64,
                Some(false_positives),
                Some(false_negatives),
            )
        }
        CompareMetric::WeightedRobinsonFoulds | CompareMetric::BranchScore => {
            let reference_splits =
                weighted_splits(reference, &taxa, options.rooted, options.include_tips)?;
            let comparison_splits =
                weighted_splits(comparison, &taxa, options.rooted, options.include_tips)?;
            let keys = reference_splits
                .keys()
                .chain(comparison_splits.keys())
                .collect::<HashSet<_>>();
            let differences = keys.into_iter().map(|key| {
                reference_splits.get(key).copied().unwrap_or(0.0)
                    - comparison_splits.get(key).copied().unwrap_or(0.0)
            });
            let value = if metric == CompareMetric::WeightedRobinsonFoulds {
                differences.map(f64::abs).sum()
            } else {
                differences
                    .map(|difference| difference * difference)
                    .sum::<f64>()
                    .sqrt()
            };
            (value, None, None)
        }
        CompareMetric::Cophenetic => (
            compare_cophenetic(reference, comparison, options.cophenetic_metric)?,
            None,
            None,
        ),
    };
    Ok(ComparisonResult {
        metric,
        value,
        false_positives,
        false_negatives,
        taxa: taxa.len(),
        rooted: options.rooted,
        include_tips: options.include_tips,
        proportion: options.proportion,
    })
}

#[derive(Clone, Debug, Eq, Hash, PartialEq)]
struct Split(Vec<u64>);

fn exact_taxa(reference: &Tree, comparison: &Tree) -> Result<Vec<String>> {
    let mut reference = reference
        .tip_index()?
        .keys()
        .map(|name| (*name).to_owned())
        .collect::<Vec<_>>();
    let mut comparison = comparison
        .tip_index()?
        .keys()
        .map(|name| (*name).to_owned())
        .collect::<Vec<_>>();
    reference.sort_unstable();
    comparison.sort_unstable();
    if reference != comparison {
        let reference_set = reference.iter().collect::<HashSet<_>>();
        let comparison_set = comparison.iter().collect::<HashSet<_>>();
        let mut reference_only = reference_set
            .difference(&comparison_set)
            .map(|name| (*name).clone())
            .collect::<Vec<_>>();
        let mut comparison_only = comparison_set
            .difference(&reference_set)
            .map(|name| (*name).clone())
            .collect::<Vec<_>>();
        reference_only.sort_unstable();
        comparison_only.sort_unstable();
        return Err(PhyloError::TaxonMismatch {
            reference_only,
            comparison_only,
        });
    }
    Ok(reference)
}

fn topology_splits(tree: &Tree, taxa: &[String], rooted: bool) -> Result<HashSet<Split>> {
    Ok(split_map(tree, taxa, rooted, false, false)?
        .into_keys()
        .collect())
}

fn weighted_splits(
    tree: &Tree,
    taxa: &[String],
    rooted: bool,
    include_tips: bool,
) -> Result<HashMap<Split, f64>> {
    split_map(tree, taxa, rooted, include_tips, true)
}

fn split_map(
    tree: &Tree,
    taxa: &[String],
    rooted: bool,
    include_tips: bool,
    weighted: bool,
) -> Result<HashMap<Split, f64>> {
    let words = taxa.len().div_ceil(64);
    let positions = taxa
        .iter()
        .enumerate()
        .map(|(index, name)| (name.as_str(), index))
        .collect::<HashMap<_, _>>();
    let mut descendants = vec![None::<Vec<u64>>; tree.node_count()];
    let mut splits = HashMap::new();
    for id in tree.postorder() {
        let node = tree.node(id).expect("node exists");
        let mut bits = vec![0u64; words];
        if node.is_tip() {
            let name = node
                .name()
                .ok_or_else(|| PhyloError::Tree(format!("unnamed tip at node {}", id.index())))?;
            let position = positions[name];
            bits[position / 64] |= 1 << (position % 64);
        } else {
            for &child in node.children() {
                let child_bits = descendants[child.index()]
                    .take()
                    .expect("postorder visits children before parents");
                for (target, source) in bits.iter_mut().zip(child_bits) {
                    *target |= source;
                }
            }
        }
        if id == tree.root() {
            descendants[id.index()] = Some(bits);
            continue;
        }
        let count = bits
            .iter()
            .map(|word| word.count_ones() as usize)
            .sum::<usize>();
        let trivial = if rooted {
            count <= 1
        } else {
            count.min(taxa.len() - count) <= 1
        };
        if trivial && !include_tips {
            descendants[id.index()] = Some(bits);
            continue;
        }
        let split_bits = if rooted {
            bits.clone()
        } else {
            canonical_bipartition(&bits, taxa.len())
        };
        let length = if weighted {
            let length = node
                .branch_length()
                .ok_or(PhyloError::MissingBranchLength(id.index()))?;
            if length < 0.0 {
                return Err(PhyloError::NegativeBranchLength {
                    node: id.index(),
                    value: length,
                });
            }
            length
        } else {
            0.0
        };
        *splits.entry(Split(split_bits)).or_insert(0.0) += length;
        descendants[id.index()] = Some(bits);
    }
    Ok(splits)
}

fn canonical_bipartition(bits: &[u64], taxa: usize) -> Vec<u64> {
    let mut complement = bits.iter().map(|word| !word).collect::<Vec<_>>();
    let remainder = taxa % 64;
    if remainder != 0 {
        *complement.last_mut().expect("at least one word") &= (1u64 << remainder) - 1;
    }
    if complement.as_slice() < bits {
        complement
    } else {
        bits.to_vec()
    }
}

fn compare_cophenetic(
    reference: &Tree,
    comparison: &Tree,
    metric: CopheneticMetric,
) -> Result<f64> {
    let reference = tip_distances(reference, true)?;
    let comparison = tip_distances(comparison, true)?;
    if reference.ids() != comparison.ids() {
        return Err(PhyloError::Tree(
            "cophenetic matrices have inconsistent identities".into(),
        ));
    }
    let n = reference.size();
    let mut left = Vec::with_capacity(n * (n - 1) / 2);
    let mut right = Vec::with_capacity(left.capacity());
    for row in 1..n {
        for column in 0..row {
            left.push(reference.get(row, column));
            right.push(comparison.get(row, column));
        }
    }
    match metric {
        CopheneticMetric::Euclidean => Ok(left
            .iter()
            .zip(&right)
            .map(|(left, right)| (left - right).powi(2))
            .sum::<f64>()
            .sqrt()),
        CopheneticMetric::CityBlock => Ok(left
            .iter()
            .zip(&right)
            .map(|(left, right)| (left - right).abs())
            .sum()),
        CopheneticMetric::UnitCorrelation => {
            let mean_left = left.iter().sum::<f64>() / left.len() as f64;
            let mean_right = right.iter().sum::<f64>() / right.len() as f64;
            let mut covariance = 0.0;
            let mut variance_left = 0.0;
            let mut variance_right = 0.0;
            for (&left, &right) in left.iter().zip(&right) {
                let centered_left = left - mean_left;
                let centered_right = right - mean_right;
                covariance += centered_left * centered_right;
                variance_left += centered_left * centered_left;
                variance_right += centered_right * centered_right;
            }
            let denominator = (variance_left * variance_right).sqrt();
            if denominator == 0.0 {
                return Err(PhyloError::UndefinedMetric(
                    "cophenetic correlation requires nonconstant path distances".into(),
                ));
            }
            Ok((1.0 - covariance / denominator) / 2.0)
        }
    }
}

#[cfg(test)]
mod tests {
    use rsomics_phylo_tree::{ParseOptions, Tree};

    use super::*;

    fn parse(source: &str) -> Tree {
        Tree::from_newick_with_options(
            source,
            ParseOptions {
                convert_underscores: false,
            },
        )
        .unwrap()
    }

    #[test]
    fn scikit_bio_rf_family_example_matches() {
        let first = parse("((a:1,b:2):1,c:4,((d:4,e:5):2,f:6):1);");
        let second = parse("((a:3,(b:2,c:2):1):3,d:8,(e:5,f:6):2);");
        let options = CompareOptions::default();
        assert_eq!(
            compare_trees(&first, &second, CompareMetric::RobinsonFoulds, options)
                .unwrap()
                .value,
            4.0
        );
        assert_eq!(
            compare_trees(
                &first,
                &second,
                CompareMetric::WeightedRobinsonFoulds,
                options,
            )
            .unwrap()
            .value,
            16.0
        );
        let branch_score = compare_trees(&first, &second, CompareMetric::BranchScore, options)
            .unwrap()
            .value;
        assert!((branch_score - 6.164_414_002_968_976).abs() < 1e-12);
    }

    #[test]
    fn directional_split_errors_have_reference_semantics() {
        let reference = parse("((a,b),(c,d));");
        let comparison = parse("((a,c),(b,d));");
        let result = compare_trees(
            &reference,
            &comparison,
            CompareMetric::FalsePositiveNegative,
            CompareOptions::default(),
        )
        .unwrap();
        assert_eq!(result.false_positives, Some(1));
        assert_eq!(result.false_negatives, Some(1));
    }
}
