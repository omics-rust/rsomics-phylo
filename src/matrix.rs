use std::collections::HashSet;
use std::io::{Read, Write};

use crate::{PhyloError, Result};

#[derive(Clone, Debug, PartialEq)]
pub struct DistanceMatrix {
    ids: Vec<String>,
    values: Vec<f64>,
}

impl DistanceMatrix {
    pub fn new(ids: Vec<String>, values: Vec<f64>) -> Result<Self> {
        if ids.len() < 2 {
            return Err(PhyloError::TooFewMatrixRows);
        }
        let expected = ids
            .len()
            .checked_mul(ids.len())
            .ok_or(PhyloError::MatrixShape {
                observed: values.len(),
                expected: usize::MAX,
            })?;
        if values.len() != expected {
            return Err(PhyloError::MatrixShape {
                observed: values.len(),
                expected,
            });
        }
        let mut unique = HashSet::with_capacity(ids.len());
        for id in &ids {
            if id.is_empty() {
                return Err(PhyloError::EmptyIdentity);
            }
            if !unique.insert(id.as_str()) {
                return Err(PhyloError::DuplicateIdentity(id.clone()));
            }
        }
        let n = ids.len();
        for row in 0..n {
            if values[row * n + row] != 0.0 {
                return Err(PhyloError::NonzeroDiagonal(row));
            }
            for column in 0..n {
                let value = values[row * n + column];
                if value.is_infinite() {
                    return Err(PhyloError::InfiniteDistance { row, column });
                }
            }
            for column in 0..row {
                let left = values[row * n + column];
                let right = values[column * n + row];
                if !(left == right || (left.is_nan() && right.is_nan())) {
                    return Err(PhyloError::AsymmetricDistance {
                        left: row,
                        right: column,
                    });
                }
            }
        }
        Ok(Self { ids, values })
    }

    pub fn from_tsv(reader: impl Read) -> Result<Self> {
        let mut input = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .flexible(false)
            .from_reader(reader);
        let header = input
            .headers()
            .map_err(table_error)?
            .iter()
            .map(str::to_owned)
            .collect::<Vec<_>>();
        if header.len() < 3 || header[0] != "id" {
            return Err(PhyloError::InvalidTable(
                "distance matrix header must be id followed by at least two identities".into(),
            ));
        }
        let ids = header[1..].to_vec();
        let n = ids.len();
        let mut values = Vec::with_capacity(n * n);
        let mut rows = 0usize;
        for record in input.records() {
            let record = record.map_err(table_error)?;
            if rows >= n {
                return Err(PhyloError::InvalidTable(format!(
                    "distance matrix has more than {n} data rows"
                )));
            }
            if record.get(0) != Some(ids[rows].as_str()) {
                return Err(PhyloError::InvalidTable(format!(
                    "row {} identity must be {:?}",
                    rows + 1,
                    ids[rows]
                )));
            }
            for column in 0..n {
                let token = record.get(column + 1).ok_or_else(|| {
                    PhyloError::InvalidTable(format!(
                        "row {} is missing column {}",
                        rows + 1,
                        column + 1
                    ))
                })?;
                values.push(token.parse::<f64>().map_err(|_| {
                    PhyloError::InvalidTable(format!(
                        "row {}, column {} is not a floating-point value: {token:?}",
                        rows + 1,
                        column + 1
                    ))
                })?);
            }
            rows += 1;
        }
        if rows != n {
            return Err(PhyloError::InvalidTable(format!(
                "distance matrix has {rows} data rows, expected {n}"
            )));
        }
        Self::new(ids, values)
    }

    pub fn write_tsv(&self, writer: impl Write) -> Result<()> {
        let mut output = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_writer(writer);
        let mut header = Vec::with_capacity(self.ids.len() + 1);
        header.push("id");
        header.extend(self.ids.iter().map(String::as_str));
        output.write_record(header).map_err(table_error)?;
        let n = self.ids.len();
        for row in 0..n {
            let mut record = Vec::with_capacity(n + 1);
            record.push(self.ids[row].clone());
            record.extend(
                self.values[row * n..(row + 1) * n]
                    .iter()
                    .map(|value| value.to_string()),
            );
            output.write_record(record).map_err(table_error)?;
        }
        output
            .flush()
            .map_err(|error| PhyloError::InvalidTable(format!("writing distance matrix: {error}")))
    }

    pub fn require_finite_nonnegative(&self) -> Result<()> {
        let n = self.ids.len();
        for row in 0..n {
            for column in 0..row {
                let value = self.values[row * n + column];
                if !value.is_finite() {
                    return Err(PhyloError::NonFiniteInferenceDistance { row, column });
                }
                if value < 0.0 {
                    return Err(PhyloError::NegativeInferenceDistance { row, column, value });
                }
            }
        }
        Ok(())
    }

    #[must_use]
    pub fn ids(&self) -> &[String] {
        &self.ids
    }

    #[must_use]
    pub fn values(&self) -> &[f64] {
        &self.values
    }

    #[must_use]
    pub fn size(&self) -> usize {
        self.ids.len()
    }

    #[must_use]
    pub fn get(&self, row: usize, column: usize) -> f64 {
        self.values[row * self.ids.len() + column]
    }
}

fn table_error(error: csv::Error) -> PhyloError {
    PhyloError::InvalidTable(error.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn labelled_tsv_round_trips_nan_symmetrically() {
        let matrix = DistanceMatrix::new(
            vec!["a".into(), "b".into()],
            vec![0.0, f64::NAN, f64::NAN, 0.0],
        )
        .unwrap();
        let mut encoded = Vec::new();
        matrix.write_tsv(&mut encoded).unwrap();
        let decoded = DistanceMatrix::from_tsv(encoded.as_slice()).unwrap();
        assert_eq!(decoded.ids(), matrix.ids());
        assert!(decoded.get(0, 1).is_nan());
    }

    #[test]
    fn row_identity_order_is_checked() {
        let input = b"id\ta\tb\nb\t0\t1\na\t1\t0\n";
        assert!(matches!(
            DistanceMatrix::from_tsv(input.as_slice()),
            Err(PhyloError::InvalidTable(_))
        ));
    }
}
