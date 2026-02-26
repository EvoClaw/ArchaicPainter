"""
S3: NEA vs DEN Attribution Test
Simulates chromosome with BOTH NEA and DEN introgression,
then evaluates ArchaicPainter's 3-state HMM source attribution accuracy.

Key question: Can the HMM distinguish NEA-origin vs DEN-origin segments?
"""
import sys, time, json, warnings
import numpy as np
from pathlib import Path
from scipy.stats import wilcoxon
warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))

import msprime
from archaicpainter.evaluation.metrics import Segment, evaluate_full
from archaicpainter.core.hmm import ArchaicHMM, decode_segments, STATE_NAMES
from archaicpainter.core.emission import (
    compute_amh_emission, compute_archaic_emission,
)

OUT_DIR = Path("/home/yanlin/livestock/docs/05_execution")

# ── Custom 3-archaic demography ────────────────────────────────────────────────
def build_nea_den_demography(f_nea=0.02, f_den=0.01):
    """
    Model with YRI (outgroup), CEU (NEA+DEN recipient), NEA, DEN.
    Topology (root = ROOT_ANC):
      ROOT_ANC ─── DEN (splits at T_split_den)
           └──── MOD_ANC ─── NEA (splits at T_split_nea)
                       └──── YRI ─── OOA ─── CEU
    CEU receives NEA admixture at T_admix_nea, DEN admixture at T_admix_den.
    """
    T_split_den = 24000   # ~700 kya
    T_split_nea = 19000   # ~550 kya
    T_OOA       = 3000    # ~87 kya
    T_admix_nea = 1724    # ~50 kya
    T_admix_den = 1500    # ~44 kya
    N = 10000

    demography = msprime.Demography()
    demography.add_population(name="CEU",  initial_size=N)
    demography.add_population(name="YRI",  initial_size=N)
    demography.add_population(name="NEA",  initial_size=2000)
    demography.add_population(name="DEN",  initial_size=2000)
    demography.add_population(name="OOA",  initial_size=N)
    demography.add_population(name="MOD",  initial_size=N)
    demography.add_population(name="ROOT", initial_size=N)

    # Pulse admixture via MassMigration (proportion of lineages moves to archaic pop)
    demography.add_event(msprime.MassMigration(
        time=T_admix_den, source="CEU", destination="DEN", proportion=f_den))
    demography.add_event(msprime.MassMigration(
        time=T_admix_nea, source="CEU", destination="NEA", proportion=f_nea))
    # Population mergers (coalescent topology, going backward in time)
    demography.add_population_split(time=T_OOA,       derived=["CEU", "OOA"], ancestral="YRI")
    demography.add_population_split(time=T_split_nea, derived=["YRI", "NEA"], ancestral="MOD")
    demography.add_population_split(time=T_split_den, derived=["MOD", "DEN"], ancestral="ROOT")

    return demography, T_admix_nea, T_admix_den, T_split_nea

def simulate_nea_den(seed, n_query=10, seq_len=5_000_000, f_nea=0.02, f_den=0.01):
    dem, T_an, T_ad, T_sn = build_nea_den_demography(f_nea, f_den)
    samples = (
        [msprime.SampleSet(n_query, population="CEU")] +
        [msprime.SampleSet(n_query, population="YRI")] +
        [msprime.SampleSet(1, population="NEA")] +
        [msprime.SampleSet(1, population="DEN")]
    )
    ts = msprime.sim_ancestry(
        samples=samples,
        demography=dem,
        sequence_length=seq_len,
        recombination_rate=1e-8,
        random_seed=seed,
        record_migrations=False,
    )
    ts_mut = msprime.sim_mutations(
        ts, rate=1.25e-8, random_seed=seed + 1, discrete_genome=True
    )
    return ts_mut, T_an, T_ad, T_sn

def extract_gt_two_source(ts, ceu_nodes, nea_node, den_node, t_split_nea, t_split_den=24000, gap=5000):
    """Extract ground truth per source (NEA vs DEN) using MRCA time."""
    per_nea = {n: [] for n in ceu_nodes}
    per_den = {n: [] for n in ceu_nodes}
    t_thresh_nea = t_split_nea * 0.5
    t_thresh_den = t_split_den * 0.5

    for tree in ts.trees():
        iv = tree.interval
        if iv.right <= iv.left: continue
        for cn in ceu_nodes:
            try:
                m_nea = tree.mrca(cn, nea_node)
                if ts.node(m_nea).time < t_thresh_nea:
                    per_nea[cn].append((iv.left, iv.right))
            except: pass
            try:
                m_den = tree.mrca(cn, den_node)
                if ts.node(m_den).time < t_thresh_den:
                    per_den[cn].append((iv.left, iv.right))
            except: pass

    def merge_ivs(ivs, gap):
        if not ivs: return []
        ivs = sorted(ivs); merged = [list(ivs[0])]
        for s, e in ivs[1:]:
            if s - merged[-1][1] <= gap: merged[-1][1] = max(merged[-1][1], e)
            else: merged.append([s, e])
        return merged

    nea_gt = {cn: [Segment(int(a), int(b), "NEA", 1.0)
                   for a, b in merge_ivs(ivs, gap)]
              for cn, ivs in per_nea.items()}
    den_gt = {cn: [Segment(int(a), int(b), "DEN", 1.0)
                   for a, b in merge_ivs(ivs, gap)]
              for cn, ivs in per_den.items()}
    return nea_gt, den_gt

def ts_to_hap_matrix(ts_mut, pop_name):
    pid = {p.metadata.get("name", str(p.id)): p.id for p in ts_mut.populations()}
    samples = list(ts_mut.samples(population=pid[pop_name]))
    gv = ts_mut.genotype_matrix(samples=samples)  # (n_variants, n_samples)
    pos = np.array([int(v.site.position) for v in ts_mut.variants()])
    # Filter: biallelic, at least 1 derived allele in this population
    mask = (gv.max(axis=1) == 1) & (gv.sum(axis=1) > 0)
    return gv[mask].T.astype(np.int8), pos[mask]  # (n_samples, n_sites), positions

def ts_to_archaic_geno(ts_mut, pop_name):
    pid = {p.metadata.get("name", str(p.id)): p.id for p in ts_mut.populations()}
    samples = list(ts_mut.samples(population=pid[pop_name]))
    gv = ts_mut.genotype_matrix(samples=samples)  # (n_variants, 2) for diploid
    pos = np.array([int(v.site.position) for v in ts_mut.variants()])
    mask = (gv.max(axis=1) == 1) & (gv.sum(axis=1) > 0)
    geno = gv[mask].T.astype(np.int8)  # (2, n_sites) → transpose to (n_sites, 2)
    return geno.T, pos[mask]  # shape: (n_sites, 2), positions

N_REPS  = 20   # 20 reps for attribution test
SEEDS   = list(range(100, 100 + N_REPS))
SEQ     = 5_000_000
N_QUERY = 10

print(f"NEA vs DEN Attribution: {SEQ//1e6:.0f}Mb × {N_REPS} reps × {N_QUERY} samples")
print("=" * 70, flush=True)

results = []

for ri, seed in enumerate(SEEDS):
    print(f"\nRep {ri+1}/{N_REPS} seed={seed}", flush=True)
    try:
        ts_mut, T_an, T_ad, T_sn = simulate_nea_den(
            seed=seed, n_query=N_QUERY, seq_len=SEQ
        )
        pid = {p.metadata.get("name", str(p.id)): p.id for p in ts_mut.populations()}
        cnodes = list(ts_mut.samples(population=pid["CEU"]))
        nea_node = list(ts_mut.samples(population=pid["NEA"]))[0]
        den_node = list(ts_mut.samples(population=pid["DEN"]))[0]

        # Haplotype matrices
        ch, cp = ts_to_hap_matrix(ts_mut, "CEU")
        yh, yp = ts_to_hap_matrix(ts_mut, "YRI")

        # Archaic genotypes
        nea_gv = ts_mut.genotype_matrix(samples=[nea_node])
        den_gv = ts_mut.genotype_matrix(samples=[den_node])
        all_pos = np.array([int(v.site.position) for v in ts_mut.variants()])
        # Use sites present in CEU
        ceu_pos_set = set(cp.tolist())
        mask = np.array([p in ceu_pos_set for p in all_pos])
        nea_hap = nea_gv[mask].reshape(-1, 1)  # (n_sites_in_ceu, 1)
        den_hap = den_gv[mask].reshape(-1, 1)

        # Build shared sites (present in CEU + at least one archaic)
        cp_set = set(cp.tolist())
        # Get CEU, NEA, DEN positions on same shared set
        all_variants_pos = all_pos[mask]
        nea_derived = (nea_hap[:, 0] > 0)
        den_derived = (den_hap[:, 0] > 0)
        any_arch = nea_derived | den_derived

        shared = all_variants_pos[any_arch]
        nea_geno_shared = np.stack([nea_hap[any_arch, 0], nea_hap[any_arch, 0]], axis=1).astype(np.int8)  # hom
        den_geno_shared = np.stack([den_hap[any_arch, 0], den_hap[any_arch, 0]], axis=1).astype(np.int8)  # hom

        # Get CEU columns for shared sites
        cp_to_idx = dict(zip(cp.tolist(), range(len(cp))))
        shared_ci = np.array([cp_to_idx[p] for p in shared if p in cp_to_idx])
        shared_in_ceu = np.array([p for p in shared if p in cp_to_idx])

        if len(shared_in_ceu) < 10:
            print(f"  Too few shared sites ({len(shared_in_ceu)}), skipping")
            continue

        # Reindex for shared_in_ceu
        any_arch_ceu = np.array([np.searchsorted(all_variants_pos, p) for p in shared_in_ceu])
        nea_geno = nea_geno_shared[np.isin(shared, shared_in_ceu)]
        den_geno = den_geno_shared[np.isin(shared, shared_in_ceu)]
        csub = ch[:, shared_ci]
        sarr = shared_in_ceu.astype(np.int32)

        # YRI frequencies at shared sites
        yp_set = set(yp.tolist()); yp_idx = dict(zip(yp.tolist(), range(len(yp))))
        ref_freq = np.array([
            (yh[:, yp_idx[p]].sum() + 0.5) / (yh.shape[0] + 1.0) if p in yp_set
            else 0.5 / (yh.shape[0] + 1.0) for p in shared_in_ceu
        ], dtype=np.float32)

        gd = np.maximum(np.diff(sarr.astype(float)) * 1e-8, 1e-12)

        # Ground truth
        nea_gt, den_gt = extract_gt_two_source(ts_mut, cnodes, nea_node, den_node, T_sn)
        n_nea_gt = sum(len(v) for v in nea_gt.values())
        n_den_gt = sum(len(v) for v in den_gt.values())
        print(f"  sites={len(sarr)} NEA_gt={n_nea_gt} DEN_gt={n_den_gt}", flush=True)

        nspec_nea = int(nea_derived.sum()); nspec_den = int(den_derived.sum())
        print(f"  NEA_specific={nspec_nea} DEN_specific={nspec_den}")

        # 3-state HMM
        hmm_3 = ArchaicHMM(pi_amh=0.97, pi_nea=0.02, pi_den=0.01)

        nea_preds_all = []; den_preds_all = []
        nea_trs_all = []; den_trs_all = []
        attr_correct = 0; attr_total = 0

        for hi, nid in enumerate(cnodes):
            q = csub[hi]
            log_emit = np.zeros((3, len(sarr)), dtype=np.float64)
            log_emit[0] = compute_amh_emission(q, ref_freq)
            log_emit[1] = compute_archaic_emission(q, nea_geno)
            log_emit[2] = compute_archaic_emission(q, den_geno)

            path = hmm_3.viterbi(log_emit, gd)
            posts, _ = hmm_3.forward_backward(log_emit, gd)
            segs = decode_segments(path, sarr, posts, STATE_NAMES, min_length_bp=5_000)

            for s in segs:
                if s["state"] == "NEA":
                    nea_preds_all.append(Segment(s["start"], s["end"], "NEA", s["mean_posterior"]))
                elif s["state"] == "DEN":
                    den_preds_all.append(Segment(s["start"], s["end"], "DEN", s["mean_posterior"]))

            nea_trs_all.extend(nea_gt.get(nid, []))
            den_trs_all.extend(den_gt.get(nid, []))

        cl = SEQ * len(cnodes)
        ev_nea = evaluate_full(nea_trs_all, nea_preds_all, cl, "NEA", ri, seed, 0)
        ev_den = evaluate_full(den_trs_all, den_preds_all, cl, "DEN", ri, seed, 0)

        # Attribution accuracy: of predicted NEA segments, how many overlap true NEA GT?
        # Uses proper interval overlap check (no sampling grid alignment issue).
        def attr_accuracy(preds, truths):
            """Fraction of predicted segments that overlap any truth segment."""
            if not preds: return 0.0
            correct = 0
            for p in preds:
                for t in truths:
                    if p.end > t.start and p.start < t.end:
                        correct += 1; break
            return correct / len(preds)

        nea_attr_acc = attr_accuracy(nea_preds_all, nea_trs_all)
        den_attr_acc = attr_accuracy(den_preds_all, den_trs_all)

        # Cross-contamination: fraction of predicted NEA that overlap DEN truth (FP attribution)
        nea_cross = attr_accuracy(nea_preds_all, den_trs_all)
        den_cross = attr_accuracy(den_preds_all, nea_trs_all)

        print(f"  NEA: F1={ev_nea.f1_all:.3f} n_pred={ev_nea.n_pred} "
              f"attr={nea_attr_acc:.3f} cross={nea_cross:.3f}")
        print(f"  DEN: F1={ev_den.f1_all:.3f} n_pred={ev_den.n_pred} "
              f"attr={den_attr_acc:.3f} cross={den_cross:.3f}")

        results.append({
            "rep": ri, "seed": seed,
            "nea_f1": ev_nea.f1_all, "nea_n_true": ev_nea.n_true, "nea_n_pred": ev_nea.n_pred,
            "den_f1": ev_den.f1_all, "den_n_true": ev_den.n_true, "den_n_pred": ev_den.n_pred,
            "nea_attr_acc": nea_attr_acc, "den_attr_acc": den_attr_acc,
            "nea_cross_contamination": nea_cross, "den_cross_contamination": den_cross,
            "n_shared": len(sarr), "nspec_nea": nspec_nea, "nspec_den": nspec_den,
        })

    except Exception as e:
        import traceback
        print(f"  ERROR: {e}")
        traceback.print_exc()

# ── Summary ────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("NEA vs DEN ATTRIBUTION SUMMARY:")
if results:
    nea_f1s   = np.array([r["nea_f1"]               for r in results])
    den_f1s   = np.array([r["den_f1"]               for r in results])
    nea_attrs = np.array([r["nea_attr_acc"]          for r in results])
    den_attrs = np.array([r["den_attr_acc"]          for r in results])
    nea_cross = np.array([r["nea_cross_contamination"] for r in results])
    den_cross = np.array([r["den_cross_contamination"] for r in results])

    print(f"  NEA F1:          {nea_f1s.mean():.3f} ± {nea_f1s.std():.3f}")
    print(f"  DEN F1:          {den_f1s.mean():.3f} ± {den_f1s.std():.3f}")
    print(f"  NEA attr.acc:    {nea_attrs.mean():.3f} ± {nea_attrs.std():.3f}")
    print(f"  DEN attr.acc:    {den_attrs.mean():.3f} ± {den_attrs.std():.3f}")
    print(f"  NEA cross-contam:{nea_cross.mean():.3f} ± {nea_cross.std():.3f}")
    print(f"  DEN cross-contam:{den_cross.mean():.3f} ± {den_cross.std():.3f}")
    print(f"  Reps completed:  {len(results)}/{N_REPS}")

    out = {"results": results, "summary": {
        "nea_f1_mean": float(nea_f1s.mean()), "nea_f1_std": float(nea_f1s.std()),
        "den_f1_mean": float(den_f1s.mean()), "den_f1_std": float(den_f1s.std()),
        "nea_attr_acc_mean": float(nea_attrs.mean()), "nea_attr_acc_std": float(nea_attrs.std()),
        "den_attr_acc_mean": float(den_attrs.mean()), "den_attr_acc_std": float(den_attrs.std()),
        "nea_cross_contamination_mean": float(nea_cross.mean()),
        "den_cross_contamination_mean": float(den_cross.mean()),
        "n_reps": len(results),
    }}
    with open(OUT_DIR / "phase4b_nea_den_attribution.json", "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved → phase4b_nea_den_attribution.json")
else:
    print("  No results collected!")
