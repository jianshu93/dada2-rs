//! Higher-order chimera (trimera) diagnostics over a sequence table.
//!
//! `remove-bimera-denovo` answers a binary question per sequence — is it a
//! chimera of two more-abundant parents? — and discards the coverage that led
//! to that answer. Some lower-abundance reads (notably on long amplicons such
//! as full-length 16S or `nodA`, and in low-biomass samples) survive bimera
//! removal yet look chimeric, consistent with being chimeras of *three or more*
//! parents. The bimera coverage signal is exactly what flags these: a read
//! whose best single-junction model nearly spans its full length leaves a small
//! internal gap that a third parent can fill.
//!
//! This module runs [`pooled_diagnostics`] over a [`SequenceTable`] and emits a
//! per-sequence TSV. It does not filter or modify the table.

use std::io::{self, Write};

use crate::chimera::{BimeraAlignParams, BimeraDiagnostic, pooled_diagnostics};
use crate::remove_bimera::BimeraParams;
use crate::sequence_table::SequenceTable;

/// One TSV row of diagnostics for a sequence.
pub struct DiagnosticRow<'a> {
    pub sequence_id: &'a str,
    pub length: usize,
    pub total_abundance: u32,
    pub n_samples: u32,
    pub is_bimera: bool,
    pub max_left: usize,
    pub max_right: usize,
    pub cover: usize,
    pub cover_frac: f64,
    pub gap_len: usize,
    pub gap_start: usize,
    pub gap_end: usize,
    pub best_left_parent: Option<&'a str>,
    pub best_right_parent: Option<&'a str>,
    pub third_parent: Option<&'a str>,
    pub gap_mismatches: usize,
    /// Mismatches the better *end* parent leaves across the gap (`usize::MAX` =
    /// not computed). The baseline a third parent must beat.
    pub gap_end_parent_mismatches: usize,
    /// Ends-free Hamming distance to the nearest single parent (`usize::MAX` =
    /// no candidate parent). Small = few-SNP variant; large = mosaic.
    pub nearest_parent_dist: usize,
    /// Refined screen (see [`TrimeraCriteria`]): the read is not a bimera, is
    /// far from any single parent, and a distinct third parent explains a
    /// non-trivial gap better than either end parent.
    pub trimera_suspect: bool,
}

/// Thresholds for the `trimera_suspect` call. Defaults are tuned on ~600 bp nodA
/// amplicons; see CLI `--trimera-*` flags.
#[derive(Clone, Copy, Debug)]
pub struct TrimeraCriteria {
    /// Minimum distance to the nearest single parent — rejects few-SNP variants.
    pub min_parent_dist: usize,
    /// Minimum gap length for a credible third segment — rejects one-off bimeras.
    pub min_gap_len: usize,
    /// Maximum third-parent mismatch fraction across the gap (clean fit).
    pub max_gap_error_frac: f64,
    /// Minimum length of *each* end flank (`max_left` and `max_right`). A
    /// 3-segment mosaic needs two real junctions, hence two substantial flanks;
    /// this rejects tiny-flank divergent singletons whose "gap" is nearly the
    /// whole read.
    pub min_flank: usize,
}

impl Default for TrimeraCriteria {
    fn default() -> Self {
        Self {
            min_parent_dist: 15,
            min_gap_len: 20,
            max_gap_error_frac: 0.10,
            min_flank: 30,
        }
    }
}

/// Run the pooled bimera coverage diagnostic over `table`.
///
/// `crit` sets the `trimera_suspect` thresholds. Returns one [`DiagnosticRow`]
/// per sequence, in table column order.
pub fn run_diagnostics<'a>(
    table: &'a SequenceTable,
    params: &BimeraParams,
    crit: TrimeraCriteria,
) -> Vec<DiagnosticRow<'a>> {
    let ncol = table.sequences.len();
    let nrow = table.samples.len();
    if ncol == 0 || nrow == 0 {
        return Vec::new();
    }

    let pooled: Vec<u32> = (0..ncol)
        .map(|j| (0..nrow).map(|i| table.counts[i][j] as u32).sum())
        .collect();
    let n_samples: Vec<u32> = (0..ncol)
        .map(|j| (0..nrow).filter(|&i| table.counts[i][j] > 0).count() as u32)
        .collect();

    let seq_bytes: Vec<&[u8]> = table.sequences.iter().map(|s| s.as_bytes()).collect();

    let align_params = BimeraAlignParams {
        allow_one_off: params.allow_one_off,
        min_one_off_par_dist: params.min_one_off_parent_distance,
        match_score: params.match_score,
        mismatch: params.mismatch,
        gap_p: params.gap_p,
        max_shift: params.max_shift,
        backend: params.backend,
        wfa_max_edits: params.wfa_max_edits,
    };

    let diags = pooled_diagnostics(
        &pooled,
        &seq_bytes,
        params.min_fold_parent_over_abundance,
        params.min_parent_abundance,
        &align_params,
    );

    let id = |k: Option<usize>| k.map(|k| table.sequence_ids[k].as_str());

    diags
        .iter()
        .enumerate()
        .map(|(j, d): (usize, &BimeraDiagnostic)| {
            let cover = (d.max_left + d.max_right).min(d.sqlen);
            let cover_frac = if d.sqlen > 0 {
                cover as f64 / d.sqlen as f64
            } else {
                0.0
            };
            let gap_len = d.gap_end.saturating_sub(d.gap_start);
            // Refined trimera screen. A genuine higher-order chimera: (1) is not
            // a bimera, (2) sits far from every single parent (not a SNP
            // variant), (3) has a non-trivial gap, (4) has a distinct third
            // parent that fills the gap cleanly AND strictly better than either
            // end parent does (else the gap is just variation an end parent
            // already explains).
            let gap_err_frac = if gap_len > 0 {
                d.gap_mismatches as f64 / gap_len as f64
            } else {
                1.0
            };
            let trimera_suspect = !d.is_bimera
                && d.nearest_parent_dist >= crit.min_parent_dist
                && d.max_left >= crit.min_flank
                && d.max_right >= crit.min_flank
                && gap_len >= crit.min_gap_len
                && d.third_parent.is_some()
                && gap_err_frac <= crit.max_gap_error_frac
                && d.gap_mismatches < d.gap_end_parent_mismatches;
            DiagnosticRow {
                sequence_id: table.sequence_ids[j].as_str(),
                length: d.sqlen,
                total_abundance: pooled[j],
                n_samples: n_samples[j],
                is_bimera: d.is_bimera,
                max_left: d.max_left,
                max_right: d.max_right,
                cover,
                cover_frac,
                gap_len,
                gap_start: d.gap_start,
                gap_end: d.gap_end,
                best_left_parent: id(d.best_left_parent),
                best_right_parent: id(d.best_right_parent),
                third_parent: id(d.third_parent),
                gap_mismatches: d.gap_mismatches,
                gap_end_parent_mismatches: d.gap_end_parent_mismatches,
                nearest_parent_dist: d.nearest_parent_dist,
                trimera_suspect,
            }
        })
        .collect()
}

const HEADER: &str = "sequence_id\tlength\ttotal_abundance\tn_samples\tis_bimera\t\
max_left\tmax_right\tcover\tcover_frac\tgap_len\tgap_start\tgap_end\t\
best_left_parent\tbest_right_parent\tthird_parent\tgap_mismatches\t\
gap_end_parent_mismatches\tnearest_parent_dist\ttrimera_suspect";

/// Format a `usize` count as a TSV cell, rendering the `usize::MAX` sentinel
/// ("not computed") as `NA`.
fn cell(v: usize) -> String {
    if v == usize::MAX {
        "NA".to_string()
    } else {
        v.to_string()
    }
}

/// Write diagnostics rows as TSV (header + one row per sequence) to `w`.
pub fn write_tsv<W: Write>(rows: &[DiagnosticRow<'_>], w: &mut W) -> io::Result<()> {
    writeln!(w, "{HEADER}")?;
    let na = "NA";
    for r in rows {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.sequence_id,
            r.length,
            r.total_abundance,
            r.n_samples,
            r.is_bimera,
            r.max_left,
            r.max_right,
            r.cover,
            r.cover_frac,
            r.gap_len,
            r.gap_start,
            r.gap_end,
            r.best_left_parent.unwrap_or(na),
            r.best_right_parent.unwrap_or(na),
            r.third_parent.unwrap_or(na),
            r.gap_mismatches,
            cell(r.gap_end_parent_mismatches),
            cell(r.nearest_parent_dist),
            r.trimera_suspect,
        )?;
    }
    Ok(())
}
