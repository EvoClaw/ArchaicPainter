"""
ArchaicPainter: main entry point.

Runs the full pipeline for a single query haplotype:
  1. Build modern reference panel allele frequencies
  2. Compute log-emission matrix (AMH | NEA | DEN)
  3. Run forward-backward + Viterbi
  4. Return posterior tracks and BED segments

For full-genome / multi-sample analysis, use run_chromosome() which handles
batching and parallelism.
"""
import numpy as np
import logging
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, field

from .emission import (
    compute_reference_panel_freq,
    build_log_emission_matrix,
)
from .hmm import ArchaicHMM, decode_segments, STATE_NAMES
from ..utils.genetic_map import compute_genetic_distances

logger = logging.getLogger(__name__)


@dataclass
class PainterResult:
    """Result for a single haplotype."""
    sample_id: str
    haplotype_idx: int          # 0 or 1 (paternal/maternal within diploid)
    chrom: str
    positions: np.ndarray       # (L,) bp
    posteriors: np.ndarray      # (N_STATES, L) float
    viterbi_path: np.ndarray    # (L,) int
    segments: list              # BED-format segments with state, length, posterior
    log_likelihood: float


class ArchaicPainter:
    """
    Archaic introgression detector using Li & Stephens haplotype matching.

    Parameters
    ----------
    hmm : ArchaicHMM
        Configured HMM instance (priors, tract length parameters).
    theta_amh : float
        Li & Stephens copying error rate for the AMH state.
    theta_archaic : float
        Substitution rate between query and archaic haplotype.
    min_segment_bp : int
        Minimum segment length to report (10 kb per Phase 4b plan).
    """

    def __init__(
        self,
        hmm: Optional[ArchaicHMM] = None,
        theta_amh: float = 0.01,
        theta_archaic: float = 0.01,
        min_segment_bp: int = 10_000,
    ):
        self.hmm = hmm if hmm is not None else ArchaicHMM()
        self.theta_amh = theta_amh
        self.theta_archaic = theta_archaic
        self.min_segment_bp = min_segment_bp

    def paint_haplotype(
        self,
        query_hap: np.ndarray,            # (L,) int8
        ref_panel: np.ndarray,            # (N_ref, L) int8
        archaic_genotypes: Dict[str, np.ndarray],  # "NEA"/"DEN" -> (L, 2)
        positions: np.ndarray,            # (L,) bp
        genetic_distances: np.ndarray,    # (L-1,) Morgans
        sample_id: str = "unknown",
        hap_idx: int = 0,
        chrom: str = "?",
    ) -> PainterResult:
        """
        Paint a single haplotype and return posterior + segments.
        """
        # Reference panel allele frequencies
        ref_freq = compute_reference_panel_freq(ref_panel)

        # Log-emission matrix (N_STATES, L)
        log_emit = build_log_emission_matrix(
            query_hap=query_hap,
            ref_panel_freq=ref_freq,
            archaic_genotypes=archaic_genotypes,
            state_order=STATE_NAMES,
            theta_amh=self.theta_amh,
            theta_archaic=self.theta_archaic,
        )

        # Forward-backward posteriors
        posteriors, log_lik = self.hmm.forward_backward(log_emit, genetic_distances)

        # Viterbi segmentation
        viterbi_path = self.hmm.viterbi(log_emit, genetic_distances)

        # BED segments
        segments = decode_segments(
            path=viterbi_path,
            positions=positions,
            posteriors=posteriors,
            state_names=STATE_NAMES,
            min_length_bp=self.min_segment_bp,
        )

        return PainterResult(
            sample_id=sample_id,
            haplotype_idx=hap_idx,
            chrom=chrom,
            positions=positions,
            posteriors=posteriors,
            viterbi_path=viterbi_path,
            segments=segments,
            log_likelihood=log_lik,
        )

    def paint_sample(
        self,
        sample_haplotypes: np.ndarray,       # (2, L) int8
        ref_panel: np.ndarray,               # (N_ref, L) int8
        archaic_genotypes: Dict[str, np.ndarray],
        positions: np.ndarray,
        genetic_distances: np.ndarray,
        sample_id: str = "unknown",
        chrom: str = "?",
    ) -> List[PainterResult]:
        """Paint both haplotypes of a diploid sample."""
        results = []
        for hap_idx in range(2):
            result = self.paint_haplotype(
                query_hap=sample_haplotypes[hap_idx],
                ref_panel=ref_panel,
                archaic_genotypes=archaic_genotypes,
                positions=positions,
                genetic_distances=genetic_distances,
                sample_id=sample_id,
                hap_idx=hap_idx,
                chrom=chrom,
            )
            results.append(result)
        return results


def results_to_bed(results: List[PainterResult], chrom: str) -> str:
    """Format list of PainterResults into BED-format string."""
    lines = []
    for res in results:
        for seg in res.segments:
            if seg["state"] in ("NEA", "DEN"):
                line = (
                    f"chr{chrom}\t{seg['start']}\t{seg['end']}\t"
                    f"{res.sample_id}_hap{res.haplotype_idx}\t"
                    f"{seg['state']}\t{seg['length_bp']}\t{seg['mean_posterior']:.4f}"
                )
                lines.append(line)
    return "\n".join(lines)
