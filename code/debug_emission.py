"""Quick emission diagnostic for ArchaicPainter on simulated data."""
import sys, numpy as np
sys.path.insert(0, "/home/yanlin/livestock/code")
import msprime
from archaicpainter.evaluation.simulation import (
    build_demography, DEMOGRAPHY_PARAMS,
    ts_to_haplotype_matrix, ts_to_archaic_genotypes
)
from archaicpainter.core.emission import (
    compute_amh_emission, compute_archaic_emission, compute_reference_panel_freq
)

p = DEMOGRAPHY_PARAMS
ts = msprime.sim_ancestry(
    samples={"CEU": 10, "YRI": 10, "NEA": 1},
    demography=build_demography(p),
    sequence_length=5_000_000,
    recombination_rate=1e-8,
    record_migrations=True,
    random_seed=42,
)
ts_mut = msprime.sim_mutations(ts, rate=1.4e-8, random_seed=42+10000, model="jc69")

pop_id = {pp.metadata.get("name", str(pp.id)): pp.id for pp in ts.populations()}

yri_hap, yri_pos = ts_to_haplotype_matrix(ts_mut, "YRI")
ceu_hap, ceu_pos = ts_to_haplotype_matrix(ts_mut, "CEU")
arch_geno, arch_pos = ts_to_archaic_genotypes(ts_mut, "NEA")

shared = sorted(set(yri_pos.tolist()) & set(ceu_pos.tolist()) & set(arch_pos.tolist()))
shared_arr = np.array(shared, dtype=np.int32)
print(f"Shared sites: {len(shared)}")

yri_idx = {pp: i for i, pp in enumerate(yri_pos.tolist())}
ceu_idx = {pp: i for i, pp in enumerate(ceu_pos.tolist())}
arc_idx = {pp: i for i, pp in enumerate(arch_pos.tolist())}
yi = np.array([yri_idx[pp] for pp in shared])
ci = np.array([ceu_idx[pp] for pp in shared])
ai = np.array([arc_idx[pp] for pp in shared])

yri_sub = yri_hap[:, yi]
ceu_sub = ceu_hap[:, ci]
arch_sub = arch_geno[ai]

ceu_nodes = list(ts.samples(population=pop_id["CEU"]))
nea_nodes = list(ts.samples(population=pop_id["NEA"]))
archaic_node = nea_nodes[0]
t_split = p["T_split_nea"]

# Ground truth for ceu_nodes[0]
gt_intervals = []
for tree in ts.trees():
    iv = tree.interval
    mrca = tree.mrca(ceu_nodes[0], archaic_node)
    if ts.node(mrca).time < t_split * 0.5:
        gt_intervals.append((iv.left, iv.right))
print(f"Ground truth intervals for hap0: {len(gt_intervals)}")

# Mark GT sites
gt_mask = np.zeros(len(shared_arr), dtype=bool)
for start, end in gt_intervals:
    in_range = (shared_arr >= start) & (shared_arr < end)
    gt_mask |= in_range
print(f"Sites in GT intervals: {gt_mask.sum()}")

# Emissions using YRI as reference panel
query_hap = ceu_sub[0]
ref_panel = yri_sub
ref_freq = compute_reference_panel_freq(ref_panel)
log_amh = compute_amh_emission(query_hap, ref_freq)
log_nea = compute_archaic_emission(query_hap, arch_sub)

bg_mask = ~gt_mask
print(f"\nNEA emission - AMH emission (positive = NEA wins):")
nea_adv = log_nea - log_amh
print(f"  GT sites:  {nea_adv[gt_mask].mean():.4f} (n={gt_mask.sum()})")
print(f"  BG sites:  {nea_adv[bg_mask].mean():.4f} (n={bg_mask.sum()})")
print(f"  NEA>AMH at GT: {(log_nea[gt_mask] > log_amh[gt_mask]).sum()}/{gt_mask.sum()}")
print(f"  NEA>AMH at BG: {(log_nea[bg_mask] > log_amh[bg_mask]).sum()}/{bg_mask.sum()}")

# Archaic genotype stats
arch_gt = arch_sub[gt_mask]
arch_bg = arch_sub[bg_mask]
def genotype_stats(ag, label):
    n_homo_alt = ((ag[:,0]==1) & (ag[:,1]==1)).sum()
    n_homo_ref = ((ag[:,0]==0) & (ag[:,1]==0)).sum()
    n_het = ((ag[:,0]!=ag[:,1]) & (ag[:,0]>=0)).sum()
    print(f"  Archaic at {label}: hom-ALT={n_homo_alt} hom-REF={n_homo_ref} het={n_het}")
genotype_stats(arch_gt, "GT sites")
genotype_stats(arch_bg, "BG sites")

# Query allele at archaic-fixed-derived sites (expected introgression signal)
arch_fixed_alt = (arch_sub[:,0]==1) & (arch_sub[:,1]==1)
query_alt = query_hap == 1
print(f"\nArchaic-fixed-derived sites: {arch_fixed_alt.sum()}")
print(f"Query carries archaic-derived allele at GT sites: "
      f"{(arch_fixed_alt & gt_mask & query_alt).sum()}/{(arch_fixed_alt & gt_mask).sum()}")
print(f"Query carries archaic-derived allele at BG sites: "
      f"{(arch_fixed_alt & bg_mask & query_alt).sum()}/{(arch_fixed_alt & bg_mask).sum()}")
