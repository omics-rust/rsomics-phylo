use rsomics_phylo_tree::Tree;

use crate::{DistanceMatrix, PhyloError, Result};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TreeMeasure {
    Sackin,
    NormalizedSackin,
    NBar,
    Colless,
    NormalizedColless,
    B1,
    Treeness,
    Gamma,
}

#[derive(Clone, Debug, PartialEq)]
pub struct NodeAge {
    node: usize,
    name: Option<String>,
    internal: bool,
    age: f64,
}

impl NodeAge {
    #[must_use]
    pub const fn node(&self) -> usize {
        self.node
    }

    #[must_use]
    pub fn name(&self) -> Option<&str> {
        self.name.as_deref()
    }

    #[must_use]
    pub const fn is_internal(&self) -> bool {
        self.internal
    }

    #[must_use]
    pub const fn age(&self) -> f64 {
        self.age
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct TreeSummary {
    pub tips: usize,
    pub nodes: usize,
    pub sackin: f64,
    pub normalized_sackin: f64,
    pub n_bar: f64,
    pub colless: f64,
    pub normalized_colless: f64,
    pub b1: f64,
    pub treeness: f64,
    pub gamma: f64,
}

pub fn tip_distances(tree: &Tree, use_length: bool) -> Result<DistanceMatrix> {
    let tip_index = tree.tip_index()?;
    let mut tips = tip_index
        .into_iter()
        .map(|(name, node)| (name.to_owned(), node))
        .collect::<Vec<_>>();
    tips.sort_unstable_by(|left, right| left.0.cmp(&right.0));
    if tips.len() < 2 {
        return Err(PhyloError::Tree(
            "tip distance requires at least two named tips".into(),
        ));
    }
    let n = tips.len();
    let mut values = vec![0.0; n * n];
    let mut tip_position = vec![None; tree.node_count()];
    for (position, (_, node)) in tips.iter().enumerate() {
        tip_position[node.index()] = Some(position);
    }
    let mut depths = vec![0.0; tree.node_count()];
    for id in tree.preorder() {
        let node = tree.node(id).expect("node exists");
        if let Some(parent) = node.parent() {
            let length = if use_length {
                let value = node
                    .branch_length()
                    .ok_or(PhyloError::MissingBranchLength(id.index()))?;
                if value < 0.0 {
                    return Err(PhyloError::NegativeBranchLength {
                        node: id.index(),
                        value,
                    });
                }
                value
            } else {
                1.0
            };
            depths[id.index()] = depths[parent.index()] + length;
        }
    }
    let tip_depths = tips
        .iter()
        .map(|(_, node)| depths[node.index()])
        .collect::<Vec<_>>();
    let mut descendants = vec![None::<Vec<usize>>; tree.node_count()];
    for id in tree.postorder() {
        let node = tree.node(id).expect("node exists");
        if let Some(position) = tip_position[id.index()] {
            descendants[id.index()] = Some(vec![position]);
            continue;
        }
        let mut groups = node
            .children()
            .iter()
            .map(|child| {
                descendants[child.index()]
                    .take()
                    .expect("postorder visits children first")
            })
            .collect::<Vec<_>>();
        groups.sort_unstable_by_key(Vec::len);
        let mut aggregate = groups.pop().unwrap_or_default();
        for group in groups {
            for &left in &aggregate {
                for &right in &group {
                    let distance = tip_depths[left] + tip_depths[right] - 2.0 * depths[id.index()];
                    values[left * n + right] = distance;
                    values[right * n + left] = distance;
                }
            }
            aggregate.extend(group);
        }
        descendants[id.index()] = Some(aggregate);
    }
    DistanceMatrix::new(tips.into_iter().map(|(name, _)| name).collect(), values)
}

pub fn tree_measure(tree: &Tree, measure: TreeMeasure, ultrametric_tolerance: f64) -> Result<f64> {
    match measure {
        TreeMeasure::Sackin => Ok(sackin(tree) as f64),
        TreeMeasure::NormalizedSackin | TreeMeasure::NBar => {
            Ok(sackin(tree) as f64 / tree.tip_count() as f64)
        }
        TreeMeasure::Colless => Ok(colless(tree)? as f64),
        TreeMeasure::NormalizedColless => {
            let tips = tree.tip_count();
            if tips < 3 {
                return Err(PhyloError::UndefinedMetric(
                    "normalized Colless requires at least three tips".into(),
                ));
            }
            Ok(colless(tree)? as f64 * 2.0 / ((tips - 1) * (tips - 2)) as f64)
        }
        TreeMeasure::B1 => Ok(b1(tree)),
        TreeMeasure::Treeness => treeness(tree),
        TreeMeasure::Gamma => gamma(tree, ultrametric_tolerance),
    }
}

pub fn node_ages(tree: &Tree, ultrametric_tolerance: f64) -> Result<Vec<NodeAge>> {
    if !ultrametric_tolerance.is_finite() || ultrametric_tolerance < 0.0 {
        return Err(PhyloError::UndefinedMetric(
            "ultrametric tolerance must be finite and nonnegative".into(),
        ));
    }
    let mut depths = vec![0.0; tree.node_count()];
    for id in tree.preorder() {
        let node = tree.node(id).expect("tree traversal yields existing nodes");
        if let Some(parent) = node.parent() {
            let length = node
                .branch_length()
                .ok_or(PhyloError::MissingBranchLength(id.index()))?;
            if length < 0.0 {
                return Err(PhyloError::NegativeBranchLength {
                    node: id.index(),
                    value: length,
                });
            }
            depths[id.index()] = depths[parent.index()] + length;
        }
    }
    let tip_depths = tree
        .preorder()
        .filter(|&id| tree.node(id).expect("node exists").is_tip())
        .map(|id| depths[id.index()])
        .collect::<Vec<_>>();
    let minimum = tip_depths.iter().copied().fold(f64::INFINITY, f64::min);
    let maximum = tip_depths.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if maximum - minimum > ultrametric_tolerance {
        return Err(PhyloError::NotUltrametric { minimum, maximum });
    }
    Ok(tree
        .preorder()
        .map(|id| {
            let node = tree.node(id).expect("node exists");
            NodeAge {
                node: id.index(),
                name: node.name().map(str::to_owned),
                internal: !node.is_tip(),
                age: maximum - depths[id.index()],
            }
        })
        .collect())
}

pub fn summarize_tree(tree: &Tree, ultrametric_tolerance: f64) -> Result<TreeSummary> {
    Ok(TreeSummary {
        tips: tree.tip_count(),
        nodes: tree.node_count(),
        sackin: tree_measure(tree, TreeMeasure::Sackin, ultrametric_tolerance)?,
        normalized_sackin: tree_measure(
            tree,
            TreeMeasure::NormalizedSackin,
            ultrametric_tolerance,
        )?,
        n_bar: tree_measure(tree, TreeMeasure::NBar, ultrametric_tolerance)?,
        colless: tree_measure(tree, TreeMeasure::Colless, ultrametric_tolerance)?,
        normalized_colless: tree_measure(
            tree,
            TreeMeasure::NormalizedColless,
            ultrametric_tolerance,
        )?,
        b1: tree_measure(tree, TreeMeasure::B1, ultrametric_tolerance)?,
        treeness: tree_measure(tree, TreeMeasure::Treeness, ultrametric_tolerance)?,
        gamma: tree_measure(tree, TreeMeasure::Gamma, ultrametric_tolerance)?,
    })
}

fn sackin(tree: &Tree) -> usize {
    let mut depth = vec![0usize; tree.node_count()];
    let mut total = 0usize;
    for id in tree.preorder() {
        let node = tree.node(id).expect("node exists");
        if let Some(parent) = node.parent() {
            depth[id.index()] = depth[parent.index()] + 1;
        }
        if node.is_tip() {
            total += depth[id.index()];
        }
    }
    total
}

fn colless(tree: &Tree) -> Result<usize> {
    let mut leaves = vec![0usize; tree.node_count()];
    let mut total = 0usize;
    for id in tree.postorder() {
        let node = tree.node(id).expect("node exists");
        if node.is_tip() {
            leaves[id.index()] = 1;
            continue;
        }
        if node.children().len() != 2 {
            return Err(PhyloError::NotBifurcating {
                node: id.index(),
                children: node.children().len(),
            });
        }
        let left = leaves[node.children()[0].index()];
        let right = leaves[node.children()[1].index()];
        total += left.abs_diff(right);
        leaves[id.index()] = left + right;
    }
    Ok(total)
}

fn b1(tree: &Tree) -> f64 {
    let mut maximum = vec![0usize; tree.node_count()];
    let mut total = 0.0;
    for id in tree.postorder() {
        let node = tree.node(id).expect("node exists");
        if node.is_tip() {
            continue;
        }
        let height = node
            .children()
            .iter()
            .map(|child| maximum[child.index()])
            .max()
            .unwrap_or(0)
            + 1;
        maximum[id.index()] = height;
        if id != tree.root() {
            total += 1.0 / height as f64;
        }
    }
    total
}

fn treeness(tree: &Tree) -> Result<f64> {
    let mut internal = 0.0;
    let mut external = 0.0;
    for edge in tree.edges() {
        let length = edge
            .length()
            .ok_or(PhyloError::MissingBranchLength(edge.child().index()))?;
        if length < 0.0 {
            return Err(PhyloError::NegativeBranchLength {
                node: edge.child().index(),
                value: length,
            });
        }
        if tree.node(edge.child()).expect("node exists").is_tip() {
            external += length;
        } else {
            internal += length;
        }
    }
    let total = internal + external;
    if total == 0.0 {
        return Err(PhyloError::UndefinedMetric(
            "treeness requires positive total branch length".into(),
        ));
    }
    Ok(internal / total)
}

fn gamma(tree: &Tree, tolerance: f64) -> Result<f64> {
    let tips = tree.tip_count();
    if tips <= 2 {
        return Err(PhyloError::UndefinedMetric(
            "gamma requires more than two tips".into(),
        ));
    }
    for id in tree.preorder() {
        let node = tree.node(id).expect("node exists");
        if !node.is_tip() && node.children().len() != 2 {
            return Err(PhyloError::NotBifurcating {
                node: id.index(),
                children: node.children().len(),
            });
        }
    }
    let mut ages = node_ages(tree, tolerance)?
        .into_iter()
        .filter(NodeAge::is_internal)
        .map(|age| age.age)
        .collect::<Vec<_>>();
    ages.sort_by(|left, right| right.total_cmp(left));
    let mut intervals = ages
        .windows(2)
        .map(|pair| pair[0] - pair[1])
        .collect::<Vec<_>>();
    intervals.push(*ages.last().expect("binary tree has an internal root"));
    let mut total = 0.0;
    let mut accumulated = 0.0;
    for i in 2..tips {
        total += i as f64 * intervals[i - 2];
        accumulated += total;
    }
    total += tips as f64 * intervals[tips - 2];
    let denominator_count = (tips - 2) as f64;
    let denominator = total * (1.0 / (12.0 * denominator_count)).sqrt();
    if denominator == 0.0 {
        return Err(PhyloError::UndefinedMetric(
            "gamma is undefined for zero temporal span".into(),
        ));
    }
    Ok((accumulated / denominator_count - total / 2.0) / denominator)
}

#[cfg(test)]
mod tests {
    use rsomics_phylo_tree::{ParseOptions, Tree};

    use super::*;

    fn tree() -> Tree {
        Tree::from_newick_with_options(
            "((a:1,b:1):1,(c:1,d:1):1);",
            ParseOptions {
                convert_underscores: false,
            },
        )
        .unwrap()
    }

    #[test]
    fn balanced_tree_statistics_match_closed_forms() {
        let tree = tree();
        assert_eq!(tree_measure(&tree, TreeMeasure::Sackin, 1e-5).unwrap(), 8.0);
        assert_eq!(
            tree_measure(&tree, TreeMeasure::Colless, 1e-5).unwrap(),
            0.0
        );
        assert_eq!(tree_measure(&tree, TreeMeasure::B1, 1e-5).unwrap(), 2.0);
        assert_eq!(
            tree_measure(&tree, TreeMeasure::Treeness, 1e-5).unwrap(),
            1.0 / 3.0
        );
    }

    #[test]
    fn tip_distances_are_identity_sorted() {
        let matrix = tip_distances(&tree(), true).unwrap();
        assert_eq!(matrix.ids(), &["a", "b", "c", "d"]);
        assert_eq!(matrix.get(0, 1), 2.0);
        assert_eq!(matrix.get(0, 2), 4.0);
    }

    #[test]
    fn non_ultrametric_tree_fails_timing() {
        let tree = Tree::from_newick("((a:1,b:2):1,c:2);").unwrap();
        assert!(matches!(
            node_ages(&tree, 1e-5),
            Err(PhyloError::NotUltrametric { .. })
        ));
    }
}
