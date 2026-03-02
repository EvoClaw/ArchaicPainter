"""
Emission probability models for ArchaicPainter.

Three-state HMM: AMH | NEA | DEN

Emission model:
  AMH:      allele frequency in modern reference panel (Li & Stephens 2003)
  NEA/DEN:  coalescent haplotype-match probability at each site

Archaic site emission:
  - Homozygous (a1 == a2 != -1): P(obs | archaic) = (1 - theta) if obs == a1, else theta
  - Heterozygous (a1 != a2):     P(obs | archaic) = 0.5  (uninformative)
  - Missing (-1):                 P(obs | archaic) = 1.0  (uninformative, no info)

AMH site emission (Li & Stephens):
  P(obs = 1 | AMH) = f_i + (1 - f_i) * theta_amh
  P(obs = 0 | AMH) = (1 - f_i) + f_i * theta_amh
  where f_i = ALT allele frequency in reference panel.

References:
  Li & Stephens (2003) Genetics 165:2213–2233
  Skov et al. (2020) Nat Genet 52:1124–1130 (HMMix, for comparison)
"""
import numpy as np
from typing import Dict, Optional
import logging

logger = logging.getLogger(__name__)

# Default error/mutation rate parameters
THETA_ARCHAIC = 0.01   # per-site substitution rate between query and archaic
THETA_AMH = 0.01       # Li & Stephens theta for modern panel emission

# State indices
STATE_AMH = 0
STATE_NEA = 1
STATE_DEN = 2
N_STATES = 3
STATE_NAMES = ["AMH", "NEA", "DEN"]


def compute_amh_emission(
    haplotype: np.ndarray,        # (L,) int8: 0/1/-1
    ref_panel_freq: np.ndarray,   # (L,) float: ALT frequency in modern reference
    theta: float = THETA_AMH,
) -> np.ndarray:
    """
    Li & Stephens emission probability for the AMH state.

    Returns log P(obs_i | AMH) for each site i.

    Uses the frequency of the ALT allele in the reference panel as the
    copy probability, with theta as the copying error rate.
    """
    L = len(haplotype)
    log_p = np.zeros(L, dtype=np.float64)

    obs = haplotype.astype(np.int8)
    freq = ref_panel_freq.astype(np.float64)

    # P(obs=1 | AMH) = freq * (1-theta) + (1-freq) * theta
    p_alt = freq * (1.0 - theta) + (1.0 - freq) * theta
    p_ref = (1.0 - freq) * (1.0 - theta) + freq * theta

    # Clamp for numerical safety
    p_alt = np.clip(p_alt, 1e-10, 1.0 - 1e-10)
    p_ref = np.clip(p_ref, 1e-10, 1.0 - 1e-10)

    valid = obs >= 0
    log_p[valid & (obs == 1)] = np.log(p_alt[valid & (obs == 1)])
    log_p[valid & (obs == 0)] = np.log(p_ref[valid & (obs == 0)])
    # Missing: log_p stays 0 (uninformative)
    return log_p


def compute_archaic_emission(
    haplotype: np.ndarray,         # (L,) int8: 0/1/-1  (query haplotype)
    archaic_geno: np.ndarray,      # (L, 2) int8: archaic diploid genotype
    theta: float = THETA_ARCHAIC,
    positive_only: bool = False,
) -> np.ndarray:
    """
    Emission probability for NEA or DEN state.

    Reflects the coalescent haplotype-match model:
    - If archaic is homozygous (a1 == a2): high probability for matching allele
    - If archaic is heterozygous (a1 != a2): 0.5 (no information about haplotype)
    - If archaic is missing: 1.0 (uninformative)

    Parameters
    ----------
    positive_only : bool
        If True, mismatches at homozygous archaic sites are treated as
        uninformative (log_p = 0) rather than penalized (log(theta)).
        Use for real-data inference where the reference archaic genome
        (e.g. Vindija) may differ from the true introgressing population.

    Returns log P(obs_i | archaic) for each site i.
    """
    L = len(haplotype)
    log_p = np.zeros(L, dtype=np.float64)

    obs = haplotype.astype(np.int8)
    a1 = archaic_geno[:, 0].astype(np.int8)
    a2 = archaic_geno[:, 1].astype(np.int8)

    missing_archaic = (a1 < 0) | (a2 < 0)
    missing_query = obs < 0
    homozygous = (~missing_archaic) & (a1 == a2)
    heterozygous = (~missing_archaic) & (~missing_query) & (a1 != a2)

    # Homozygous archaic sites
    match = homozygous & (~missing_query) & (obs == a1)
    mismatch = homozygous & (~missing_query) & (obs != a1)
    log_p[match] = np.log(1.0 - theta)
    # Mismatches: penalize in simulation mode; down-weight in real-data mode.
    # positive_only sets log P(mismatch|NEA) = log(0.5), so mismatches provide
    # mild negative evidence for NEA rather than strong penalisation.
    # This matches Eq.(emiss_nea_pos) in the paper: ell(mismatch) = log(1/2).
    if not positive_only:
        log_p[mismatch] = np.log(theta)
    else:
        log_p[mismatch] = np.log(0.5)

    # Heterozygous archaic: 50/50, uninformative
    log_p[heterozygous] = np.log(0.5)

    # Missing: stays 0 (log(1) = uninformative)
    return log_p


def compute_reference_panel_freq(
    ref_haplotypes: np.ndarray,    # (N_ref_hap, L) int8
) -> np.ndarray:
    """
    Compute ALT allele frequency in the modern reference panel at each site.
    Excludes missing data from the denominator.
    """
    valid_mask = ref_haplotypes >= 0
    n_valid = valid_mask.sum(axis=0).astype(np.float64)
    n_alt = (ref_haplotypes == 1).sum(axis=0).astype(np.float64)

    # Laplace smoothing to avoid 0/0
    freq = (n_alt + 0.5) / (n_valid + 1.0)
    return freq.astype(np.float32)


def build_log_emission_matrix(
    query_hap: np.ndarray,          # (L,)
    ref_panel_freq: np.ndarray,     # (L,) ALT freq in modern panel
    archaic_genotypes: Dict[str, np.ndarray],  # name -> (L, 2)
    state_order: Optional[list] = None,
    theta_amh: float = THETA_AMH,
    theta_archaic: float = THETA_ARCHAIC,
) -> np.ndarray:
    """
    Compute the full (N_states, L) log-emission matrix.

    Parameters
    ----------
    query_hap : array (L,)
        Query haplotype alleles.
    ref_panel_freq : array (L,)
        ALT allele frequency in modern reference panel.
    archaic_genotypes : dict
        Keys like "NEA" or "DEN", values are (L, 2) arrays.
    state_order : list, optional
        Order of states, default ["AMH", "NEA", "DEN"].

    Returns
    -------
    log_emit : array (N_states, L)
        Log emission probabilities.
    """
    if state_order is None:
        state_order = STATE_NAMES

    L = len(query_hap)
    n_states = len(state_order)
    log_emit = np.zeros((n_states, L), dtype=np.float64)

    for s_idx, state in enumerate(state_order):
        if state == "AMH":
            log_emit[s_idx] = compute_amh_emission(
                query_hap, ref_panel_freq, theta=theta_amh
            )
        elif state in archaic_genotypes:
            log_emit[s_idx] = compute_archaic_emission(
                query_hap, archaic_genotypes[state], theta=theta_archaic
            )
        else:
            logger.warning(f"No genotype data for state '{state}'; using zero emission")

    return log_emit
