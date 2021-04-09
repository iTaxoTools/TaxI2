//! Calculating distances between sequences

use ndarray::{s, Array2, Array3, Zip};

use crate::needle::Aligner;

/// State for the distance calculation
pub struct AlignmentStats {
    total_length: usize,
    common_length: usize,
    total_gap_length: usize,
    transitions: usize,
    transversions: usize,
}

impl AlignmentStats {
    /// Zeroed state.
    pub fn new() -> Self {
        AlignmentStats {
            total_length: 0,
            common_length: 0,
            total_gap_length: 0,
            transitions: 0,
            transversions: 0,
        }
    }

    // Number of substitions calculated so far
    fn substitions(&self) -> usize {
        self.transversions + self.transitions
    }

    /// pairwise uncorrelated distance
    pub fn pdistance(&self) -> f64 {
        (self.substitions() as f64 / self.common_length as f64).abs()
    }

    /// pairwise uncorrelated distance with gaps
    pub fn pdistance_counting_gaps(&self) -> f64 {
        f64::abs(
            (self.substitions() as f64 + self.total_gap_length as f64) / self.total_length as f64,
        )
    }

    /// Jukes-Cantor distance
    pub fn jukes_cantor_distance(&self) -> f64 {
        let p = self.substitions() as f64 / self.common_length as f64;
        if p > 3.0 / 4.0 {
            f64::INFINITY
        } else {
            f64::abs(-(3.0 / 4.0) * f64::ln_1p(-(4.0 / 3.0) * p))
        }
    }

    /// Kimura's two parameter distance
    pub fn kimura2p_distance(&self) -> f64 {
        let p = self.transitions as f64 / self.common_length as f64;
        let q = self.transversions as f64 / self.common_length as f64;
        let distance =
            f64::abs(-(1.0 / 2.0) * f64::ln((1.0 - 2.0 * p - q) * f64::sqrt(1.0 - 2.0 * q)));
        if distance.is_nan() {
            f64::INFINITY
        } else {
            distance
        }
    }

    fn count_gap(&mut self) {
        self.total_length += 1;
        self.total_gap_length += 1;
    }

    fn count_match(&mut self) {
        self.total_length += 1;
        self.common_length += 1;
    }

    fn count_transition(&mut self) {
        self.total_length += 1;
        self.common_length += 1;
        self.transitions += 1;
    }
    fn count_transversion(&mut self) {
        self.total_length += 1;
        self.common_length += 1;
        self.transversions += 1;
    }

    /// Count `(x, y)` pair.
    pub fn update(&mut self, (x, y): (u8, u8)) {
        use NucleotideType::*;
        use SymbolType::*;
        match (classify(x), classify(y)) {
            (Gap, Nucleotide(_)) => self.count_gap(),
            (Nucleotide(_), Gap) => self.count_gap(),
            (Nucleotide(_), Nucleotide(_)) if x == y => self.count_match(),
            (Nucleotide(Purine), Nucleotide(Purine))
            | (Nucleotide(Pyrimidine), Nucleotide(Pyrimidine)) => self.count_transition(),
            (Nucleotide(_), Nucleotide(_)) => self.count_transversion(),
            _ => {}
        }
    }
}

enum SymbolType {
    Gap,
    Missing,
    Nucleotide(NucleotideType),
}

enum NucleotideType {
    Purine,
    Pyrimidine,
    Unknown,
}

fn classify(x: u8) -> SymbolType {
    use NucleotideType::*;
    use SymbolType::*;
    match x {
        b'-' => Gap,
        b'n' | b'N' | b'?' => Missing,
        b'a' | b'A' | b'g' | b'G' => Nucleotide(Purine),
        b'c' | b'C' | b't' | b'T' => Nucleotide(Pyrimidine),
        _ => Nucleotide(Unknown),
    }
}

/// Returns 4 distances between `target` and `query`.
///
/// Performs alignment.
pub fn seq_distances(aligner: &Aligner, target: &str, query: &str) -> [f64; 4] {
    let alignment = aligner.align(target.as_bytes(), query.as_bytes());
    let mut alignment_stats = AlignmentStats::new();
    alignment
        .common_path_iter()
        .for_each(|pair| alignment_stats.update(pair));
    [
        alignment_stats.pdistance(),
        alignment_stats.jukes_cantor_distance(),
        alignment_stats.kimura2p_distance(),
        alignment_stats.pdistance_counting_gaps(),
    ]
}

// Returns true if the character is part of a meaningful part of a sequences
fn is_nucleotide(c: char) -> bool {
    !matches!(c, '-' | 'n' | 'N' | '?')
}

// Returns the inclusive boundaries of the common non-gap part of given sequences
fn common_content(target: &str, query: &str) -> Option<(usize, usize)> {
    let target_start = target.find(is_nucleotide)?;
    let query_start = query.find(is_nucleotide)?;
    let target_end = target.rfind(is_nucleotide)?;
    let query_end = query.rfind(is_nucleotide)?;
    let start = usize::max(target_start, query_start);
    let end = usize::min(target_end, query_end);
    Some((start, end))
}

/// Returns 4 distances between `target` and `query`.
///
/// Expects aligned sequences.
pub fn seq_distances_aligned(target: &str, query: &str) -> [f64; 4] {
    let (start, end) = match common_content(target, query) {
        None => return [f64::NAN; 4],
        Some(x) => x,
    };
    let target = &target[start..=end];
    let query = &query[start..=end];
    let mut alignment_stats = AlignmentStats::new();
    target
        .bytes()
        .zip(query.bytes())
        .for_each(|pair| alignment_stats.update(pair));
    [
        alignment_stats.pdistance(),
        alignment_stats.jukes_cantor_distance(),
        alignment_stats.kimura2p_distance(),
        alignment_stats.pdistance_counting_gaps(),
    ]
}

/// Creates (n, 4) array of distances between `targets` and `queries`.
///
/// Outer iteration over `targets`.
/// Inner iteration over `queries`.
/// Performs sequence-to-sequence alignment
pub fn make_distance_array(aligner: &Aligner, targets: &[&str], queries: &[&str]) -> Array2<f64> {
    use std::mem::MaybeUninit;

    let mut distance_table = Array3::<f64>::uninit((targets.len(), queries.len(), 4));

    targets.iter().enumerate().for_each(|(i, target)| {
        queries.iter().enumerate().for_each(|(j, query)| {
            Zip::from(distance_table.slice_mut(s![i, j, ..]))
                .and(&seq_distances(aligner, target, query)[..])
                .for_each(|slot, &distance| *slot = MaybeUninit::new(distance))
        })
    });

    let distance_table = unsafe { distance_table.assume_init() };
    distance_table
        .into_shape((targets.len() * queries.len(), 4))
        .expect("Can't put into the correct shape")
}

pub fn make_distance_array_aligned(targets: &[&str], queries: &[&str]) -> Array2<f64> {
    use std::mem::MaybeUninit;

    let mut distance_table = Array3::<f64>::uninit((targets.len(), queries.len(), 4));

    targets.iter().enumerate().for_each(|(i, target)| {
        queries.iter().enumerate().for_each(|(j, query)| {
            Zip::from(distance_table.slice_mut(s![i, j, ..]))
                .and(&seq_distances_aligned(target, query)[..])
                .for_each(|slot, &distance| *slot = MaybeUninit::new(distance))
        })
    });

    let distance_table = unsafe { distance_table.assume_init() };
    distance_table
        .into_shape((targets.len() * queries.len(), 4))
        .expect("Can't put into the correct shape")
}

#[cfg(test)]
mod test_super {
    use ndarray::array;

    use super::*;

    #[test]
    fn test_distance() {
        let target = "gg-ccnccta";
        let query = "ggaccaccaa";
        let mut alignment_stats = AlignmentStats::new();
        target
            .bytes()
            .zip(query.bytes())
            .for_each(|pair| alignment_stats.update(pair));
        assert_eq!(alignment_stats.pdistance(), 1.0 / 8.0);
        assert_eq!(alignment_stats.pdistance_counting_gaps(), 2.0 / 9.0);
    }

    #[test]
    fn test_distance_table() {
        let targets = ["foo", "fao", "f-o"];
        let queries = ["foo", "bar"];

        let aligner = Aligner::default();
        let distance_table = make_distance_array(&aligner, &targets, &queries);
        let pdistances = array![0.0, 1.0, 1.0 / 3.0, 2.0 / 3.0, 0.0, 1.0];
        assert_eq!(distance_table.column(0), pdistances);
        let pdistances_gaps = array![0.0, 1.0, 1.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0, 1.0];
        assert_eq!(distance_table.column(3), pdistances_gaps);
    }
}
