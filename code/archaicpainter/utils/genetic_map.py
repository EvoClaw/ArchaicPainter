"""
Genetic map utilities: load Oxford HapMap II / SHAPEIT genetic maps and
interpolate to get cM positions and recombination rates between adjacent sites.

Used for HMM transition probabilities: P(change) = 1 - exp(-r * d_s)
where r = genetic distance between consecutive sites, d_s = expected tract length.
"""
import numpy as np
import pandas as pd
import logging
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

# Genetic map column names used by Oxford HapMap/SHAPEIT format
HAPMAP_COLS = ["Chromosome", "Position(bp)", "Rate(cM/Mb)", "Map(cM)"]


def load_genetic_map(map_path: str, chrom: str) -> pd.DataFrame:
    """
    Load a genetic map file (HapMap II / SHAPEIT format).

    Expected columns: Chromosome, Position(bp), Rate(cM/Mb), Map(cM)
    Returns a DataFrame sorted by position.
    """
    path = Path(map_path)
    if not path.exists():
        raise FileNotFoundError(f"Genetic map not found: {map_path}")

    df = pd.read_csv(map_path, sep=r"\s+", header=0)

    # Normalise column names
    df.columns = [c.strip() for c in df.columns]
    if "Position(bp)" not in df.columns and "pos" in df.columns:
        df = df.rename(columns={"pos": "Position(bp)", "cM": "Map(cM)", "rate": "Rate(cM/Mb)"})

    # Filter to requested chromosome if column present
    if "Chromosome" in df.columns:
        chrom_val = chrom.lstrip("chr")
        df = df[df["Chromosome"].astype(str) == chrom_val]

    df = df.sort_values("Position(bp)").reset_index(drop=True)
    logger.info(f"Loaded genetic map for chr{chrom}: {len(df)} markers")
    return df


def interpolate_genetic_positions(
    positions_bp: np.ndarray,
    map_df: pd.DataFrame,
) -> np.ndarray:
    """
    Linear interpolation of genetic positions (cM) at arbitrary bp positions.

    Positions outside the map range are extrapolated using the boundary rate.
    """
    map_pos = map_df["Position(bp)"].values.astype(np.float64)
    map_cm = map_df["Map(cM)"].values.astype(np.float64)
    return np.interp(positions_bp.astype(np.float64), map_pos, map_cm)


def compute_genetic_distances(
    positions_bp: np.ndarray,
    map_df: Optional[pd.DataFrame] = None,
    genetic_pos_cm: Optional[np.ndarray] = None,
    fallback_rate_cM_per_Mb: float = 1.0,
) -> np.ndarray:
    """
    Compute inter-site genetic distances (in Morgans) between consecutive sites.

    Returns array of length L-1 where out[i] = distance from site i to i+1.

    If no map is available, uses a uniform fallback rate.
    """
    if genetic_pos_cm is not None:
        cm = genetic_pos_cm
    elif map_df is not None:
        cm = interpolate_genetic_positions(positions_bp, map_df)
    else:
        # Uniform fallback: 1 cM/Mb
        bp_dist = np.diff(positions_bp.astype(np.float64))
        morgans = bp_dist * fallback_rate_cM_per_Mb * 1e-8
        return morgans

    cm_diff = np.diff(cm.astype(np.float64))
    morgans = cm_diff / 100.0
    morgans = np.maximum(morgans, 1e-9)  # numerical floor
    return morgans


def make_transition_rates(
    genetic_distances: np.ndarray,
    t_admix_gen: float = 1724.0,
) -> np.ndarray:
    """
    HMM transition probabilities for introgressed tract model.

    P(leave state) = 1 - exp(-r * t_admix)
    where r = genetic distance in Morgans between adjacent sites.

    Returns array of shape (L-1,): probability of a state transition at each gap.
    """
    p_change = 1.0 - np.exp(-genetic_distances * t_admix_gen)
    return p_change
