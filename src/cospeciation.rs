use std::collections::HashSet;

use crate::{DistanceMatrix, PhyloError, Result};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct InteractionMatrix {
    associate_ids: Vec<String>,
    host_ids: Vec<String>,
    values: Vec<bool>,
}

impl InteractionMatrix {
    pub fn new(
        associate_ids: Vec<String>,
        host_ids: Vec<String>,
        values: Vec<bool>,
    ) -> Result<Self> {
        let expected = associate_ids.len().checked_mul(host_ids.len()).ok_or(
            PhyloError::InteractionShape {
                rows: associate_ids.len(),
                columns: host_ids.len(),
                expected_rows: usize::MAX,
                expected_columns: usize::MAX,
            },
        )?;
        if values.len() != expected {
            let columns = if associate_ids.is_empty() {
                0
            } else {
                values.len() / associate_ids.len()
            };
            return Err(PhyloError::InteractionShape {
                rows: associate_ids.len(),
                columns,
                expected_rows: associate_ids.len(),
                expected_columns: host_ids.len(),
            });
        }
        validate_ids(&associate_ids)?;
        validate_ids(&host_ids)?;
        Ok(Self {
            associate_ids,
            host_ids,
            values,
        })
    }

    #[must_use]
    pub fn associate_ids(&self) -> &[String] {
        &self.associate_ids
    }

    #[must_use]
    pub fn host_ids(&self) -> &[String] {
        &self.host_ids
    }

    #[must_use]
    pub fn values(&self) -> &[bool] {
        &self.values
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct HommolaOptions {
    pub permutations: usize,
    pub seed: u64,
    pub store_null: bool,
}

impl Default for HommolaOptions {
    fn default() -> Self {
        Self {
            permutations: 999,
            seed: 0,
            store_null: false,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct HommolaResult {
    pub correlation: f64,
    pub p_value: Option<f64>,
    pub permutations: usize,
    pub seed: u64,
    pub null_distribution: Vec<f64>,
}

pub fn hommola_cospeciation(
    hosts: &DistanceMatrix,
    associates: &DistanceMatrix,
    interactions: &InteractionMatrix,
    options: HommolaOptions,
) -> Result<HommolaResult> {
    if hosts.size() < 3 || associates.size() < 3 {
        return Err(PhyloError::TooFewInteractions);
    }
    if interactions.host_ids != hosts.ids() || interactions.associate_ids != associates.ids() {
        return Err(PhyloError::InteractionIdentityMismatch);
    }
    hosts.require_finite_nonnegative()?;
    associates.require_finite_nonnegative()?;
    let edges = interactions
        .values
        .iter()
        .enumerate()
        .filter(|&(_, &present)| present)
        .map(|(index, _)| (index / hosts.size(), index % hosts.size()))
        .collect::<Vec<_>>();
    if edges.len() < 3 {
        return Err(PhyloError::TooFewInteractions);
    }
    let correlation = edge_correlation(hosts, associates, &edges, None, None)?;
    if options.permutations == 0 {
        return Ok(HommolaResult {
            correlation,
            p_value: None,
            permutations: 0,
            seed: options.seed,
            null_distribution: Vec::new(),
        });
    }
    let permutation_denominator = options
        .permutations
        .checked_add(1)
        .ok_or(PhyloError::TooManyPermutations)?;
    let mut rng = SplitMix64::new(options.seed);
    let mut associate_order = (0..associates.size()).collect::<Vec<_>>();
    let mut host_order = (0..hosts.size()).collect::<Vec<_>>();
    let mut extreme = 0usize;
    let mut null_distribution = if options.store_null {
        Vec::with_capacity(options.permutations)
    } else {
        Vec::new()
    };
    for _ in 0..options.permutations {
        rng.shuffle(&mut associate_order);
        rng.shuffle(&mut host_order);
        let statistic = edge_correlation(
            hosts,
            associates,
            &edges,
            Some(&host_order),
            Some(&associate_order),
        )?;
        extreme += usize::from(statistic >= correlation);
        if options.store_null {
            null_distribution.push(statistic);
        }
    }
    Ok(HommolaResult {
        correlation,
        p_value: Some((extreme + 1) as f64 / permutation_denominator as f64),
        permutations: options.permutations,
        seed: options.seed,
        null_distribution,
    })
}

fn validate_ids(ids: &[String]) -> Result<()> {
    let mut unique = HashSet::with_capacity(ids.len());
    for id in ids {
        if id.is_empty() {
            return Err(PhyloError::EmptyIdentity);
        }
        if !unique.insert(id.as_str()) {
            return Err(PhyloError::DuplicateIdentity(id.clone()));
        }
    }
    Ok(())
}

fn edge_correlation(
    hosts: &DistanceMatrix,
    associates: &DistanceMatrix,
    edges: &[(usize, usize)],
    host_order: Option<&[usize]>,
    associate_order: Option<&[usize]>,
) -> Result<f64> {
    let pair_count = edges.len() * (edges.len() - 1) / 2;
    let mut host_sum = 0.0;
    let mut associate_sum = 0.0;
    for left in 0..edges.len() - 1 {
        for right in left + 1..edges.len() {
            let (left_associate, left_host) = edges[left];
            let (right_associate, right_host) = edges[right];
            host_sum += hosts.get(
                permuted(left_host, host_order),
                permuted(right_host, host_order),
            );
            associate_sum += associates.get(
                permuted(left_associate, associate_order),
                permuted(right_associate, associate_order),
            );
        }
    }
    let host_mean = host_sum / pair_count as f64;
    let associate_mean = associate_sum / pair_count as f64;
    let mut covariance = 0.0;
    let mut host_variance = 0.0;
    let mut associate_variance = 0.0;
    for left in 0..edges.len() - 1 {
        for right in left + 1..edges.len() {
            let (left_associate, left_host) = edges[left];
            let (right_associate, right_host) = edges[right];
            let host = hosts.get(
                permuted(left_host, host_order),
                permuted(right_host, host_order),
            ) - host_mean;
            let associate = associates.get(
                permuted(left_associate, associate_order),
                permuted(right_associate, associate_order),
            ) - associate_mean;
            covariance += host * associate;
            host_variance += host * host;
            associate_variance += associate * associate;
        }
    }
    let denominator = (host_variance * associate_variance).sqrt();
    if denominator == 0.0 {
        return Err(PhyloError::UndefinedMetric(
            "Hommola correlation requires variable paired distances".into(),
        ));
    }
    Ok(covariance / denominator)
}

fn permuted(index: usize, order: Option<&[usize]>) -> usize {
    order.map_or(index, |order| order[index])
}

struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    const fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut value = self.state;
        value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        value ^ (value >> 31)
    }

    fn index(&mut self, bound: usize) -> usize {
        let bound = bound as u64;
        let zone = u64::MAX - u64::MAX % bound;
        loop {
            let value = self.next();
            if value < zone {
                return (value % bound) as usize;
            }
        }
    }

    fn shuffle(&mut self, values: &mut [usize]) {
        for upper in (1..values.len()).rev() {
            let lower = self.index(upper + 1);
            values.swap(upper, lower);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn matrix(values: Vec<f64>) -> DistanceMatrix {
        DistanceMatrix::new(
            vec!["a".into(), "b".into(), "c".into(), "d".into(), "e".into()],
            values,
        )
        .unwrap()
    }

    #[test]
    fn published_hommola_example_matches_observed_statistic() {
        let hosts = matrix(vec![
            0., 3., 8., 8., 9., 3., 0., 7., 7., 8., 8., 7., 0., 6., 7., 8., 7., 6., 0., 3., 9., 8.,
            7., 3., 0.,
        ]);
        let associates = matrix(vec![
            0., 5., 8., 8., 8., 5., 0., 7., 7., 7., 8., 7., 0., 4., 4., 8., 7., 4., 0., 2., 8., 7.,
            4., 2., 0.,
        ]);
        let interactions = InteractionMatrix::new(
            associates.ids().to_vec(),
            hosts.ids().to_vec(),
            vec![
                true, false, false, false, false, false, true, false, false, false, false, false,
                true, false, false, false, false, false, true, false, false, false, false, true,
                true,
            ],
        )
        .unwrap();
        let result = hommola_cospeciation(
            &hosts,
            &associates,
            &interactions,
            HommolaOptions {
                permutations: 0,
                seed: 42,
                store_null: false,
            },
        )
        .unwrap();
        assert!((result.correlation - 0.831_709_654_632_479_1).abs() < 1e-12);
        assert_eq!(result.p_value, None);
    }

    #[test]
    fn native_permutations_are_reproducible_without_storing_nulls() {
        let values = vec![0., 1., 2., 1., 0., 3., 2., 3., 0.];
        let hosts =
            DistanceMatrix::new(vec!["h1".into(), "h2".into(), "h3".into()], values.clone())
                .unwrap();
        let associates =
            DistanceMatrix::new(vec!["p1".into(), "p2".into(), "p3".into()], values).unwrap();
        let interactions = InteractionMatrix::new(
            associates.ids().to_vec(),
            hosts.ids().to_vec(),
            vec![true, false, false, false, true, false, false, false, true],
        )
        .unwrap();
        let options = HommolaOptions {
            permutations: 20,
            seed: 7,
            store_null: false,
        };
        let first = hommola_cospeciation(&hosts, &associates, &interactions, options).unwrap();
        let second = hommola_cospeciation(&hosts, &associates, &interactions, options).unwrap();
        assert_eq!(first, second);
        assert!(first.null_distribution.is_empty());
    }
}
