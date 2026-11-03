use rsomics_phylo_tree::{Tree, TreeBuilder};

use crate::{DistanceMatrix, PhyloError, Result};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum InferenceMethod {
    NeighborJoining,
    Upgma,
    Wpgma,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum NegativeBranchPolicy {
    Zero,
    Allow,
}

pub fn infer_tree(
    matrix: &DistanceMatrix,
    method: InferenceMethod,
    negative_branches: NegativeBranchPolicy,
) -> Result<Tree> {
    matrix.require_finite_nonnegative()?;
    match method {
        InferenceMethod::NeighborJoining => neighbor_joining(matrix, negative_branches),
        InferenceMethod::Upgma => pair_group(matrix, PairGroupWeight::Taxa),
        InferenceMethod::Wpgma => {
            if matrix.size() < 3 {
                return Err(PhyloError::TooFewWpgmaTaxa);
            }
            pair_group(matrix, PairGroupWeight::Clusters)
        }
    }
}

#[derive(Debug)]
enum Shape {
    Leaf(String),
    Internal(Vec<(Shape, f64)>),
}

#[derive(Clone, Copy)]
enum PairGroupWeight {
    Taxa,
    Clusters,
}

fn pair_group(matrix: &DistanceMatrix, weight: PairGroupWeight) -> Result<Tree> {
    let n = matrix.size();
    let mut condensed = Vec::with_capacity(n * (n - 1) / 2);
    for row in 0..n {
        for column in row + 1..n {
            condensed.push(matrix.get(row, column));
        }
    }
    let merges = pair_group_linkage(&condensed, n, weight);
    let mut shapes = matrix
        .ids()
        .iter()
        .cloned()
        .map(Shape::Leaf)
        .map(Some)
        .chain(std::iter::repeat_with(|| None).take(n - 1))
        .collect::<Vec<_>>();
    let mut heights = vec![0.0; 2 * n - 1];
    for (index, merge) in merges.into_iter().enumerate() {
        let node = n + index;
        let height = merge.distance / 2.0;
        let left = shapes[merge.left].take().expect("linkage child exists");
        let right = shapes[merge.right].take().expect("linkage child exists");
        shapes[node] = Some(Shape::Internal(vec![
            (left, (height - heights[merge.left]).max(0.0)),
            (right, (height - heights[merge.right]).max(0.0)),
        ]));
        heights[node] = height;
    }
    build_tree(
        shapes[2 * n - 2]
            .take()
            .expect("final linkage cluster exists"),
    )
}

#[derive(Clone, Copy)]
struct LinkageMerge {
    left: usize,
    right: usize,
    distance: f64,
}

fn pair_group_linkage(condensed: &[f64], n: usize, weight: PairGroupWeight) -> Vec<LinkageMerge> {
    let mut distances = condensed.to_vec();
    let mut sizes = vec![1usize; n];
    let mut chain = vec![0usize; n];
    let mut chain_length = 0usize;
    let mut raw = Vec::with_capacity(n - 1);
    for _ in 0..n - 1 {
        if chain_length == 0 {
            chain[0] = sizes
                .iter()
                .position(|&size| size > 0)
                .expect("an active linkage cluster exists");
            chain_length = 1;
        }
        let (mut left, mut right, distance) = loop {
            let left = chain[chain_length - 1];
            let (mut right, mut distance) = if chain_length > 1 {
                let right = chain[chain_length - 2];
                (right, distances[condensed_index(n, left, right)])
            } else {
                (0, f64::INFINITY)
            };
            for candidate in 0..n {
                if sizes[candidate] == 0 || candidate == left {
                    continue;
                }
                let candidate_distance = distances[condensed_index(n, left, candidate)];
                if candidate_distance < distance {
                    right = candidate;
                    distance = candidate_distance;
                }
            }
            if chain_length > 1 && right == chain[chain_length - 2] {
                break (left, right, distance);
            }
            chain[chain_length] = right;
            chain_length += 1;
        };
        chain_length -= 2;
        if left > right {
            std::mem::swap(&mut left, &mut right);
        }
        let left_size = sizes[left];
        let right_size = sizes[right];
        raw.push((left, right, distance));
        sizes[left] = 0;
        sizes[right] = left_size + right_size;
        for other in 0..n {
            if sizes[other] == 0 || other == right {
                continue;
            }
            let left_distance = distances[condensed_index(n, other, left)];
            let right_distance = distances[condensed_index(n, other, right)];
            let updated = match weight {
                PairGroupWeight::Taxa => {
                    (left_distance * left_size as f64 + right_distance * right_size as f64)
                        / (left_size + right_size) as f64
                }
                PairGroupWeight::Clusters => (left_distance + right_distance) / 2.0,
            };
            let index = condensed_index(n, other, right);
            distances[index] = updated;
        }
    }
    let mut order = (0..raw.len()).collect::<Vec<_>>();
    order.sort_by(|&left, &right| raw[left].2.total_cmp(&raw[right].2));
    let mut union = LinkageUnion::new(n);
    order
        .into_iter()
        .map(|index| {
            let (left, right, distance) = raw[index];
            let left = union.find(left);
            let right = union.find(right);
            let (left, right) = if left < right {
                (left, right)
            } else {
                (right, left)
            };
            union.merge(left, right);
            LinkageMerge {
                left,
                right,
                distance,
            }
        })
        .collect()
}

fn condensed_index(n: usize, left: usize, right: usize) -> usize {
    let (left, right) = if left < right {
        (left, right)
    } else {
        (right, left)
    };
    n * left - left * (left + 1) / 2 + right - left - 1
}

fn condensed_get(distances: &[f64], capacity: usize, left: usize, right: usize) -> f64 {
    if left == right {
        0.0
    } else {
        distances[condensed_index(capacity, left, right)]
    }
}

fn condensed_set(distances: &mut [f64], capacity: usize, left: usize, right: usize, value: f64) {
    let index = condensed_index(capacity, left, right);
    distances[index] = value;
}

struct LinkageUnion {
    parent: Vec<usize>,
    next: usize,
}

impl LinkageUnion {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..2 * n - 1).collect(),
            next: n,
        }
    }

    fn find(&mut self, value: usize) -> usize {
        let mut root = value;
        while self.parent[root] != root {
            root = self.parent[root];
        }
        let mut current = value;
        while self.parent[current] != root {
            let next = self.parent[current];
            self.parent[current] = root;
            current = next;
        }
        root
    }

    fn merge(&mut self, left: usize, right: usize) {
        self.parent[left] = self.next;
        self.parent[right] = self.next;
        self.next += 1;
    }
}

fn neighbor_joining(
    matrix: &DistanceMatrix,
    negative_branches: NegativeBranchPolicy,
) -> Result<Tree> {
    let n = matrix.size();
    let capacity = 2 * n - 2;
    let mut distances = vec![0.0; capacity * (capacity - 1) / 2];
    for row in 0..n {
        for column in row + 1..n {
            condensed_set(
                &mut distances,
                capacity,
                row,
                column,
                matrix.get(row, column),
            );
        }
    }
    let mut shapes = matrix
        .ids()
        .iter()
        .cloned()
        .map(Shape::Leaf)
        .map(Some)
        .chain(std::iter::repeat_with(|| None).take(n - 2))
        .collect::<Vec<_>>();
    let mut active = vec![false; capacity];
    active[..n].fill(true);
    let mut sums = vec![0.0; capacity];
    for (row, sum) in sums.iter_mut().take(n).enumerate() {
        *sum = (0..n)
            .map(|column| condensed_get(&distances, capacity, row, column))
            .sum();
    }
    let mut sorted = vec![Vec::<u32>::new(); capacity];
    for (row, order) in sorted.iter_mut().take(n).enumerate() {
        *order = (0..n)
            .filter(|&column| column != row)
            .map(|column| column as u32)
            .collect();
        order.sort_unstable_by(|left, right| {
            condensed_get(&distances, capacity, row, *left as usize).total_cmp(&condensed_get(
                &distances,
                capacity,
                row,
                *right as usize,
            ))
        });
    }
    let mut count = n;
    let mut remaining = n;
    while remaining > 3 {
        let inverse = 1.0 / (remaining - 2) as f64;
        let maximum_u = (0..count)
            .filter(|&index| active[index])
            .map(|index| sums[index] * inverse)
            .fold(f64::NEG_INFINITY, f64::max);
        let mut best = (usize::MAX, usize::MAX);
        let mut best_q = f64::INFINITY;
        for left in 0..count {
            if !active[left] {
                continue;
            }
            let left_u = sums[left] * inverse;
            for &right in &sorted[left] {
                let right = right as usize;
                let distance = condensed_get(&distances, capacity, left, right);
                if distance - left_u - maximum_u >= best_q {
                    break;
                }
                if !active[right] {
                    continue;
                }
                let q = distance - left_u - sums[right] * inverse;
                if q < best_q {
                    best = (left, right);
                    best_q = q;
                }
            }
        }
        let (left_id, right_id) = best;
        let between = condensed_get(&distances, capacity, left_id, right_id);
        let delta = (sums[left_id] - sums[right_id]) * inverse;
        let left_length = branch((between + delta) / 2.0, negative_branches);
        let right_length = branch((between - delta) / 2.0, negative_branches);
        let left = shapes[left_id].take().expect("active shape exists");
        let right = shapes[right_id].take().expect("active shape exists");
        shapes[count] = Some(Shape::Internal(vec![
            (left, left_length),
            (right, right_length),
        ]));
        let mut new_sum = 0.0;
        for other in 0..count {
            if !active[other] || other == left_id || other == right_id {
                continue;
            }
            let left_distance = condensed_get(&distances, capacity, other, left_id);
            let right_distance = condensed_get(&distances, capacity, other, right_id);
            let value = (left_distance + right_distance - between) / 2.0;
            sums[other] = sums[other] - left_distance - right_distance + value;
            condensed_set(&mut distances, capacity, count, other, value);
            new_sum += value;
        }
        sums[count] = new_sum;
        active[left_id] = false;
        active[right_id] = false;
        active[count] = true;
        sorted[count] = (0..count)
            .filter(|&other| active[other])
            .map(|other| other as u32)
            .collect();
        sorted[count].sort_unstable_by(|left, right| {
            condensed_get(&distances, capacity, count, *left as usize).total_cmp(&condensed_get(
                &distances,
                capacity,
                count,
                *right as usize,
            ))
        });
        count += 1;
        remaining -= 1;
    }

    let rest = (0..count)
        .filter(|&index| active[index])
        .collect::<Vec<_>>();
    let children = if rest.len() == 2 {
        let distance = condensed_get(&distances, capacity, rest[0], rest[1]);
        vec![
            (
                shapes[rest[0]].take().expect("active shape exists"),
                branch(distance / 2.0, negative_branches),
            ),
            (
                shapes[rest[1]].take().expect("active shape exists"),
                branch(distance / 2.0, negative_branches),
            ),
        ]
    } else {
        let a = rest[0];
        let b = rest[1];
        let c = rest[2];
        let ab = condensed_get(&distances, capacity, a, b);
        let ac = condensed_get(&distances, capacity, a, c);
        let bc = condensed_get(&distances, capacity, b, c);
        let a_length = branch((ab + ac - bc) / 2.0, negative_branches);
        let b_length = branch((ab + bc - ac) / 2.0, negative_branches);
        let c_length = branch((ac + bc - ab) / 2.0, negative_branches);
        vec![
            (shapes[a].take().expect("active shape exists"), a_length),
            (shapes[b].take().expect("active shape exists"), b_length),
            (shapes[c].take().expect("active shape exists"), c_length),
        ]
    };
    build_tree(Shape::Internal(children))
}

fn branch(length: f64, policy: NegativeBranchPolicy) -> f64 {
    match policy {
        NegativeBranchPolicy::Zero => length.max(0.0),
        NegativeBranchPolicy::Allow => length,
    }
}

fn build_tree(root: Shape) -> Result<Tree> {
    let Shape::Internal(children) = root else {
        return Err(PhyloError::Tree(
            "inferred tree has no internal root".into(),
        ));
    };
    let mut builder = TreeBuilder::new(None, None)?;
    let root_id = builder.root();
    let mut pending = children
        .into_iter()
        .rev()
        .map(|(shape, length)| (root_id, shape, length))
        .collect::<Vec<_>>();
    while let Some((parent, shape, length)) = pending.pop() {
        match shape {
            Shape::Leaf(name) => {
                builder.add_child(parent, Some(&name), Some(length))?;
            }
            Shape::Internal(children) => {
                let node = builder.add_child(parent, None, Some(length))?;
                pending.extend(
                    children
                        .into_iter()
                        .rev()
                        .map(|(shape, length)| (node, shape, length)),
                );
            }
        }
    }
    Ok(builder.finish())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn three() -> DistanceMatrix {
        DistanceMatrix::new(
            vec!["a".into(), "b".into(), "c".into()],
            vec![0.0, 2.0, 4.0, 2.0, 0.0, 4.0, 4.0, 4.0, 0.0],
        )
        .unwrap()
    }

    #[test]
    fn upgma_is_rooted_and_ultrametric() {
        let tree =
            infer_tree(&three(), InferenceMethod::Upgma, NegativeBranchPolicy::Zero).unwrap();
        assert_eq!(tree.to_newick(), "(c:2,(a:1,b:1):1);");
    }

    #[test]
    fn nj_uses_a_three_way_unrooted_final_join() {
        let tree = infer_tree(
            &three(),
            InferenceMethod::NeighborJoining,
            NegativeBranchPolicy::Zero,
        )
        .unwrap();
        assert_eq!(tree.node(tree.root()).unwrap().children().len(), 3);
        assert_eq!(tree.tip_index().unwrap().len(), 3);
    }
}
