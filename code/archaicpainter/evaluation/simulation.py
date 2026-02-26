"""
Simulation-based ground truth generation for ArchaicPainter benchmarking.
Uses msprime + Phase 4a validated demography.
"""
import numpy as np
import msprime
import logging
from typing import List, Dict, Optional, Tuple
from .metrics import Segment

logger = logging.getLogger(__name__)

DEMOGRAPHY_PARAMS = {
    "N_CEU": 512_000,
    "N_YRI": 512_000,
    "N_NEA": 1_000,
    "N_ancestral": 7_300,
    "T_OOA": 2_000,
    "T_admix_nea": 1_724,
    "T_split_nea": 18_966,
    "T_ancient": 40_000,
    "f_nea": 0.02,
    "recomb_rate": 1e-8,
    "mut_rate": 1.4e-8,
}

def build_demography(params=None):
    p = params or DEMOGRAPHY_PARAMS
    d = msprime.Demography()
    d.add_population(name="YRI",     initial_size=p["N_YRI"])
    d.add_population(name="CEU",     initial_size=p["N_CEU"])
    d.add_population(name="NEA",     initial_size=p["N_NEA"])
    d.add_population(name="OOA",     initial_size=p["N_ancestral"])
    d.add_population(name="ANCIENT", initial_size=p["N_ancestral"])
    d.add_mass_migration(time=p["T_admix_nea"], source="CEU", dest="NEA", proportion=p["f_nea"])
    d.add_population_split(time=p["T_OOA"], derived=["CEU", "YRI"], ancestral="OOA")
    d.add_population_split(time=p["T_split_nea"], derived=["OOA", "NEA"], ancestral="ANCIENT")
    return d

def simulate_one(seed, n_query=20, n_archaic=1, seq_len=5_000_000, params=None):
    p = params or DEMOGRAPHY_PARAMS
    ts = msprime.sim_ancestry(
        samples={"CEU": n_query, "YRI": n_query, "NEA": n_archaic},
        demography=build_demography(p),
        sequence_length=seq_len,
        recombination_rate=p["recomb_rate"],
        record_migrations=True,
        random_seed=seed,
    )
    ts_mut = msprime.sim_mutations(ts, rate=p["mut_rate"], random_seed=seed+10_000, model="jc69")
    return ts, ts_mut

def extract_true_segments(ts, t_admix=1_724):
    segs = []
    for mig in ts.migrations():
        slen = mig.right - mig.left
        if slen > 0 and abs(mig.time - t_admix) < t_admix * 0.8:
            segs.append(Segment(start=int(mig.left), end=int(mig.right), source="NEA", posterior=1.0))
    return segs

def ts_to_haplotype_matrix(ts_mut, pop_name="CEU"):
    """Return (hap_matrix, positions) for sites where pop has at least one derived allele."""
    pop_id = {p.metadata.get("name", str(p.id)): p.id for p in ts_mut.populations()}
    target_pop = pop_id[pop_name]
    nodes = [s for s in ts_mut.samples() if ts_mut.node(s).population == target_pop]
    pos, rows = [], []
    for v in ts_mut.variants(samples=nodes):
        geno = v.genotypes.astype(np.int8)
        if geno.sum() > 0:  # at least one sample carries the derived allele
            pos.append(int(v.site.position))
            rows.append(geno)
    hap = np.stack(rows, axis=1) if rows else np.zeros((len(nodes), 0), dtype=np.int8)
    return hap, np.array(pos, dtype=np.int32)

def ts_to_archaic_genotypes(ts_mut, pop_name="NEA", sample_idx=0):
    """Return (geno_matrix, positions) for sites where the archaic individual carries derived allele."""
    pop_id = {p.metadata.get("name", str(p.id)): p.id for p in ts_mut.populations()}
    target_pop = pop_id[pop_name]
    nodes = [s for s in ts_mut.samples() if ts_mut.node(s).population == target_pop]
    h0, h1 = nodes[2*sample_idx], nodes[2*sample_idx+1]
    pos, rows = [], []
    for v in ts_mut.variants(samples=[h0, h1]):
        g0, g1 = int(v.genotypes[0]), int(v.genotypes[1])
        if g0 > 0 or g1 > 0:  # archaic has at least one derived allele
            pos.append(int(v.site.position))
            rows.append([g0, g1])
    geno = np.array(rows, dtype=np.int8) if rows else np.zeros((0, 2), dtype=np.int8)
    return geno, np.array(pos, dtype=np.int32)
