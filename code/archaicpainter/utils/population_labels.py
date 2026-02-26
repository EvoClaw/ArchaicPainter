"""
1000GP population metadata: load sample-to-population mapping,
define population groups for analysis (AMR, AFR, EUR, SAS, EAS).
"""
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional

# All 26 1000GP populations
POPULATIONS = {
    "AFR": ["YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB"],
    "AMR": ["MXL", "PUR", "CLM", "PEL"],
    "EAS": ["CHB", "JPT", "CHS", "CDX", "KHV"],
    "EUR": ["CEU", "TSI", "FIN", "GBR", "IBS"],
    "SAS": ["GIH", "PJL", "BEB", "STU", "ITU"],
}

# Populations used as modern human reference panel (exclude query populations)
DEFAULT_REFERENCE_POPS = ["YRI", "LWK", "GWD", "MSL", "ESN"]  # West African outgroup

# Primary analysis populations per method-design.md
PRIMARY_QUERY_POPS = {
    "AMR_admixed": ["PEL", "MXL", "CLM", "PUR"],   # biological discovery target
    "EUR_control": ["CEU", "FIN", "GBR"],            # known introgression benchmarks
    "EAS_control": ["CHB", "CHS", "JPT"],            # high archaic ancestry controls
    "AFR_negative": ["YRI", "LWK"],                  # low/no Neanderthal introgression
}


def load_sample_panel(panel_path: str) -> pd.DataFrame:
    """
    Load 1000GP integrated_call_samples_v3.20130502.ALL.panel or equivalent.
    Expected columns: sample, pop, super_pop, gender.
    """
    path = Path(panel_path)
    if not path.exists():
        raise FileNotFoundError(f"Panel file not found: {panel_path}")

    df = pd.read_csv(panel_path, sep=r"\s+")
    df.columns = [c.lower() for c in df.columns]
    return df


def get_samples_by_pop(panel_df: pd.DataFrame, pops: List[str]) -> List[str]:
    """Return list of sample IDs belonging to specified populations."""
    mask = panel_df["pop"].isin(pops)
    return panel_df.loc[mask, "sample"].tolist()


def get_samples_by_superpop(panel_df: pd.DataFrame, super_pops: List[str]) -> List[str]:
    """Return list of sample IDs belonging to specified super-populations."""
    col = "super_pop" if "super_pop" in panel_df.columns else "super_pop"
    mask = panel_df[col].isin(super_pops)
    return panel_df.loc[mask, "sample"].tolist()


def pop_to_superpop(panel_df: pd.DataFrame) -> Dict[str, str]:
    """Return mapping from population code to super-population code."""
    col = "super_pop" if "super_pop" in panel_df.columns else "super_pop"
    return dict(zip(panel_df["pop"], panel_df[col]))
