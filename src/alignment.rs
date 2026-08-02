use std::collections::HashSet;

use crate::{PhyloError, Result};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AlignedRecord {
    identity: String,
    header: String,
    sequence: Vec<u8>,
}

impl AlignedRecord {
    pub fn new(header: impl Into<String>, sequence: Vec<u8>) -> Result<Self> {
        let header = header.into();
        let identity = header
            .split_ascii_whitespace()
            .next()
            .ok_or(PhyloError::EmptyIdentity)?
            .to_owned();
        if let Some(&byte) = sequence
            .iter()
            .find(|&&byte| !(b'!'..=b'~').contains(&byte))
        {
            return Err(PhyloError::InvalidAlignmentByte { id: identity, byte });
        }
        Ok(Self {
            identity,
            header,
            sequence,
        })
    }

    #[must_use]
    pub fn identity(&self) -> &str {
        &self.identity
    }

    #[must_use]
    pub fn header(&self) -> &str {
        &self.header
    }

    #[must_use]
    pub fn sequence(&self) -> &[u8] {
        &self.sequence
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Alignment {
    records: Vec<AlignedRecord>,
    width: usize,
}

impl Alignment {
    pub fn new(records: Vec<AlignedRecord>) -> Result<Self> {
        if records.len() < 2 {
            return Err(PhyloError::TooFewSequences);
        }
        let width = records[0].sequence.len();
        if width == 0 {
            return Err(PhyloError::EmptySequence);
        }
        let mut identities = HashSet::with_capacity(records.len());
        for record in &records {
            if !identities.insert(record.identity.as_str()) {
                return Err(PhyloError::DuplicateIdentity(record.identity.clone()));
            }
            if record.sequence.len() != width {
                return Err(PhyloError::RaggedAlignment {
                    id: record.identity.clone(),
                    observed: record.sequence.len(),
                    expected: width,
                });
            }
        }
        Ok(Self { records, width })
    }

    #[must_use]
    pub fn records(&self) -> &[AlignedRecord] {
        &self.records
    }

    #[must_use]
    pub const fn width(&self) -> usize {
        self.width
    }

    pub fn trim_by_gap(&self, maximum_gap_fraction: f64) -> Result<TrimmedAlignment> {
        if !maximum_gap_fraction.is_finite() || !(0.0..=1.0).contains(&maximum_gap_fraction) {
            return Err(PhyloError::InvalidGapThreshold(maximum_gap_fraction));
        }
        let count = self.records.len();
        let kept_columns = (0..self.width)
            .filter(|&column| {
                let gaps = self
                    .records
                    .iter()
                    .filter(|record| matches!(record.sequence[column], b'-' | b'.'))
                    .count();
                gaps < count && (gaps as f64) / (count as f64) <= maximum_gap_fraction
            })
            .collect::<Vec<_>>();
        if kept_columns.is_empty() {
            return Err(PhyloError::EmptyTrimResult);
        }
        let records = self
            .records
            .iter()
            .map(|record| AlignedRecord {
                identity: record.identity.clone(),
                header: record.header.clone(),
                sequence: kept_columns
                    .iter()
                    .map(|&column| record.sequence[column])
                    .collect(),
            })
            .collect();
        Ok(TrimmedAlignment {
            alignment: Self {
                records,
                width: kept_columns.len(),
            },
            kept_columns,
        })
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TrimmedAlignment {
    alignment: Alignment,
    kept_columns: Vec<usize>,
}

impl TrimmedAlignment {
    #[must_use]
    pub fn alignment(&self) -> &Alignment {
        &self.alignment
    }

    #[must_use]
    pub fn kept_columns(&self) -> &[usize] {
        &self.kept_columns
    }

    #[must_use]
    pub fn into_parts(self) -> (Alignment, Vec<usize>) {
        (self.alignment, self.kept_columns)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn alignment() -> Alignment {
        Alignment::new(vec![
            AlignedRecord::new("a first", b"A-.T".to_vec()).unwrap(),
            AlignedRecord::new("b", b"AC-T".to_vec()).unwrap(),
            AlignedRecord::new("c", b"ACGT".to_vec()).unwrap(),
        ])
        .unwrap()
    }

    #[test]
    fn trim_retains_original_headers_and_coordinates() {
        let result = alignment().trim_by_gap(0.34).unwrap();
        assert_eq!(result.kept_columns(), &[0, 1, 3]);
        assert_eq!(result.alignment().records()[0].header(), "a first");
        let strict = alignment().trim_by_gap(0.0).unwrap();
        assert_eq!(strict.kept_columns(), &[0, 3]);
        assert_eq!(strict.alignment().records()[1].sequence(), b"AT");
    }

    #[test]
    fn duplicate_first_token_is_rejected() {
        let error = Alignment::new(vec![
            AlignedRecord::new("a one", b"AC".to_vec()).unwrap(),
            AlignedRecord::new("a two", b"AC".to_vec()).unwrap(),
        ])
        .unwrap_err();
        assert_eq!(error, PhyloError::DuplicateIdentity("a".into()));
    }

    #[test]
    fn all_gap_columns_are_never_retained() {
        let alignment = Alignment::new(vec![
            AlignedRecord::new("a", b"A-".to_vec()).unwrap(),
            AlignedRecord::new("b", b"C-".to_vec()).unwrap(),
        ])
        .unwrap();
        assert_eq!(alignment.trim_by_gap(1.0).unwrap().kept_columns(), &[0]);
    }
}
