use rayon::prelude::*;

use crate::{Alignment, DistanceMatrix, Result};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SequenceMetric {
    Hamming,
    Jc69,
    K2p,
}

pub fn sequence_distances(alignment: &Alignment, metric: SequenceMetric) -> Result<DistanceMatrix> {
    let records = alignment.records();
    let n = records.len();
    let encoded = (metric != SequenceMetric::Hamming).then(|| {
        records
            .iter()
            .map(|record| {
                record
                    .sequence()
                    .iter()
                    .copied()
                    .map(code)
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>()
    });
    let lower = (0..n)
        .into_par_iter()
        .map(|row| {
            (0..row)
                .map(|column| match metric {
                    SequenceMetric::Hamming => {
                        hamming(records[row].sequence(), records[column].sequence())
                    }
                    SequenceMetric::Jc69 | SequenceMetric::K2p => {
                        let encoded = encoded
                            .as_ref()
                            .expect("corrected metrics have encoded sequences");
                        corrected(&encoded[row], &encoded[column], metric)
                    }
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let mut values = vec![0.0; n * n];
    for row in 1..n {
        for column in 0..row {
            let value = lower[row][column];
            values[row * n + column] = value;
            values[column * n + row] = value;
        }
    }
    DistanceMatrix::new(
        records
            .iter()
            .map(|record| record.identity().to_owned())
            .collect(),
        values,
    )
}

fn hamming(left: &[u8], right: &[u8]) -> f64 {
    let mismatches = left.iter().zip(right).filter(|&(a, b)| a != b).count();
    mismatches as f64 / left.len() as f64
}

fn corrected(left: &[u8], right: &[u8], metric: SequenceMetric) -> f64 {
    let (transitions, transversions, valid) = substitution_counts(left, right);
    if valid == 0 {
        return f64::NAN;
    }
    match metric {
        SequenceMetric::Hamming => unreachable!("Hamming does not use canonical encoding"),
        SequenceMetric::Jc69 => {
            let difference = (transitions + transversions) as f64 / valid as f64;
            let argument = 1.0 - (4.0 / 3.0) * difference;
            if argument <= 0.0 {
                f64::NAN
            } else {
                -0.75 * argument.ln() + 0.0
            }
        }
        SequenceMetric::K2p => {
            let p = transitions as f64 / valid as f64;
            let q = transversions as f64 / valid as f64;
            let first = 1.0 - 2.0 * p - q;
            let second = 1.0 - 2.0 * q;
            if first <= 0.0 || second <= 0.0 {
                f64::NAN
            } else {
                -0.5 * first.ln() - 0.25 * second.ln() + 0.0
            }
        }
    }
}

fn substitution_counts(left: &[u8], right: &[u8]) -> (usize, usize, usize) {
    let mut transitions = 0usize;
    let mut transversions = 0usize;
    let mut valid = 0usize;
    for (&left, &right) in left.iter().zip(right) {
        if left == 4 || right == 4 {
            continue;
        }
        valid += 1;
        let substitution = left ^ right;
        if substitution == 2 {
            transitions += 1;
        } else if substitution & 1 == 1 {
            transversions += 1;
        }
    }
    (transitions, transversions, valid)
}

fn code(base: u8) -> u8 {
    match base {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' | b'U' => 3,
        _ => 4,
    }
}

#[cfg(test)]
mod tests {
    use crate::AlignedRecord;

    use super::*;

    fn pair(left: &[u8], right: &[u8]) -> Alignment {
        Alignment::new(vec![
            AlignedRecord::new("left", left.to_vec()).unwrap(),
            AlignedRecord::new("right", right.to_vec()).unwrap(),
        ])
        .unwrap()
    }

    #[test]
    fn hamming_uses_all_sites_literally() {
        let matrix =
            sequence_distances(&pair(b"acgu-N", b"ACGT?N"), SequenceMetric::Hamming).unwrap();
        assert_eq!(matrix.get(0, 1), 5.0 / 6.0);
    }

    #[test]
    fn corrected_metrics_delete_noncanonical_sites_pairwise() {
        let alignment = pair(b"ACGT-N", b"AGGT?A");
        let jc = sequence_distances(&alignment, SequenceMetric::Jc69).unwrap();
        let expected_jc = -0.75 * (1.0_f64 - 4.0 / 3.0 * 0.25).ln();
        assert!((jc.get(0, 1) - expected_jc).abs() < 1e-12);

        let k2p = sequence_distances(&alignment, SequenceMetric::K2p).unwrap();
        let expected_k2p = -0.5 * (1.0_f64 - 0.25).ln() - 0.25 * (1.0_f64 - 2.0 * 0.25).ln();
        assert!((k2p.get(0, 1) - expected_k2p).abs() < 1e-12);
    }

    #[test]
    fn corrected_metrics_treat_lowercase_as_noncanonical() {
        let matrix = sequence_distances(&pair(b"acgt", b"ACGT"), SequenceMetric::Jc69).unwrap();
        assert!(matrix.get(0, 1).is_nan());
    }
}
