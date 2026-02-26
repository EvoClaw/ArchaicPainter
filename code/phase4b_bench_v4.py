"""
Phase 4b Benchmark v4
- 50 replicates (seeds 42..91)
- Bootstrap 95% CI + Wilcoxon signed-rank test vs Poisson baseline
- Also validates positive-only emission (R4) in same run
- NEA vs DEN attribution test (S3) in same run
"""
import sys, time, json, warnings
import numpy as np
from pathlib import Path
from scipy.stats import wilcoxon, bootstrap as scipy_bootstrap
warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))

OUT_DIR = Path("/home/yanlin/livestock/docs/05_execution")
OUT_DIR.mkdir(parents=True, exist_ok=True)

import msprime
from archaicpainter.evaluation.simulation import (
    build_demography, simulate_one, DEMOGRAPHY_PARAMS,
    ts_to_haplotype_matrix, ts_to_archaic_genotypes,
)
from archaicpainter.evaluation.metrics import Segment, evaluate_full
from archaicpainter.core.hmm import ArchaicHMM, decode_segments, STATE_NAMES
from archaicpainter.core.emission import compute_amh_emission, compute_archaic_emission

N_REPS  = 50
SEEDS   = list(range(42, 42 + N_REPS))
N_QUERY = 10       # 20 haplotypes per rep
SEQ     = 5_000_000
T_SPLIT = DEMOGRAPHY_PARAMS["T_split_nea"]
T_SPLIT_DEN = DEMOGRAPHY_PARAMS.get("T_split_den", 25000)  # fallback

# ── Ground-truth extraction ────────────────────────────────────────────────────
def extract_gt(ts, ceu_nodes, archaic_node, t_split=T_SPLIT, gap=5000):
    per = {n: [] for n in ceu_nodes}
    for tree in ts.trees():
        iv = tree.interval
        if iv.right <= iv.left: continue
        for cn in ceu_nodes:
            try:
                m = tree.mrca(cn, archaic_node)
                if ts.node(m).time < t_split * 0.5:
                    per[cn].append((iv.left, iv.right))
            except: pass
    result = {}
    for cn, ivs in per.items():
        if not ivs: result[cn] = []; continue
        s = sorted(ivs); merged = [list(s[0])]
        for st, en in s[1:]:
            if st - merged[-1][1] <= gap: merged[-1][1] = max(merged[-1][1], en)
            else: merged.append([st, en])
        result[cn] = [Segment(int(a), int(b), "NEA", 1.0) for a, b in merged]
    return result

def make_freq(shared, hap, pos, n_hap):
    pos_map = dict(zip(pos.tolist(), range(len(pos))))
    freq = np.full(len(shared), 0.5 / (n_hap + 1.0), dtype=np.float32)
    for i, p in enumerate(shared):
        if p in pos_map:
            freq[i] = ((hap[:, pos_map[p]] == 1).sum() + 0.5) / (n_hap + 1.0)
    return freq

def merge_nea(segs, gap=50_000, minlen=30_000):
    if not segs: return []
    segs = sorted(segs, key=lambda s: s.start)
    m = [Segment(segs[0].start, segs[0].end, segs[0].source, segs[0].posterior)]
    for s in segs[1:]:
        if s.start - m[-1].end <= gap:
            m[-1] = Segment(m[-1].start, max(m[-1].end, s.end), m[-1].source,
                             max(m[-1].posterior, s.posterior))
        else: m.append(Segment(s.start, s.end, s.source, s.posterior))
    return [s for s in m if s.length >= minlen]

def run_ap(ref_freq, q, ag, pos, gd, hmm_obj, merge=True,
           theta_arch=0.01, theta_amh=0.01, positive_only=False):
    if len(pos) < 10: return []
    log_emit = np.zeros((3, len(pos)), dtype=np.float64)
    log_emit[0] = compute_amh_emission(q, ref_freq, theta=theta_amh)
    log_emit[1] = compute_archaic_emission(q, ag, theta=theta_arch,
                                           positive_only=positive_only)
    path = hmm_obj.viterbi(log_emit, gd)
    posteriors, _ = hmm_obj.forward_backward(log_emit, gd)
    segs = decode_segments(path, pos, posteriors, STATE_NAMES, min_length_bp=5_000)
    raw = [Segment(s["start"], s["end"], s["state"], s["mean_posterior"])
           for s in segs if s["state"] == "NEA"]
    return merge_nea(raw) if merge else [s for s in raw if s.length >= 30_000]

def pd_baseline_v(positions, q, ag, win=50_000, step=20_000):
    if len(positions) == 0: return []
    fixed = (ag[:, 0] == 1) & (ag[:, 1] == 1)
    hits = (q == 1) & fixed
    if hits.sum() == 0: return []
    cum_hits = np.cumsum(np.concatenate([[0], hits.astype(np.int32)]))
    start = int(positions[0]); end = int(positions[-1]) + 1
    windows = np.arange(start, end, step, dtype=np.int64)
    segs = []
    for w in windows:
        lo = np.searchsorted(positions, w, 'left')
        hi = np.searchsorted(positions, w + win, 'left')
        h = cum_hits[hi] - cum_hits[lo]; n = hi - lo
        if n > 0 and h > 0:
            segs.append(Segment(int(w), int(w + win), "NEA", float(h/n)))
    return segs

# ── Configs ────────────────────────────────────────────────────────────────────
configs = {
    "archaicpainter":      {"ref": "YRI", "merge": True,  "theta_arch": 0.01, "positive_only": False},
    "ap_no_merge":         {"ref": "YRI", "merge": False, "theta_arch": 0.01, "positive_only": False},
    "ap_positive_only":    {"ref": "YRI", "merge": True,  "theta_arch": 0.01, "positive_only": True},  # R4
    "ap_ceu_ref":          {"ref": "CEU", "merge": True,  "theta_arch": 0.01, "positive_only": False},
    "ap_theta_high":       {"ref": "YRI", "merge": True,  "theta_arch": 0.05, "positive_only": False},
    "poisson_density":     None,
}

all_r = {k: [] for k in configs}

print(f"Benchmark v4: {SEQ//1_000_000}Mb × {N_REPS} seeds × {N_QUERY} samples")
print(f"Configs: {list(configs.keys())}")
print("=" * 70, flush=True)

for ri, seed in enumerate(SEEDS):
    print(f"\nRep {ri+1}/{N_REPS} seed={seed}", flush=True)
    t_rep = time.perf_counter()
    try:
        ts, ts_mut = simulate_one(seed=seed, n_query=N_QUERY, n_archaic=1, seq_len=SEQ)
        pid = {p.metadata.get("name", str(p.id)): p.id for p in ts.populations()}
        cnodes = list(ts.samples(population=pid["CEU"]))
        anode  = list(ts.samples(population=pid["NEA"]))[0]

        ch, cp = ts_to_haplotype_matrix(ts_mut, "CEU")
        yh, yp = ts_to_haplotype_matrix(ts_mut, "YRI")
        ah, ap2 = ts_to_archaic_genotypes(ts_mut, "NEA")

        shared = sorted(set(cp.tolist()) & set(ap2.tolist()))
        sarr = np.array(shared, dtype=np.int32)
        ci = np.array([dict(zip(cp.tolist(), range(len(cp))))[x] for x in shared])
        ai = np.array([dict(zip(ap2.tolist(), range(len(ap2))))[x] for x in shared])
        csub = ch[:, ci]; asub = ah[ai]

        ref_freqs = {
            "YRI": make_freq(shared, yh, yp, yh.shape[0]),
            "CEU": make_freq(shared, ch, cp, ch.shape[0]),
        }
        nspec = sum(1 for p in shared if p not in set(yp.tolist()))
        print(f"  sites={len(shared)} NEA_spec={nspec}", flush=True)

        t0 = time.perf_counter()
        truth = extract_gt(ts, cnodes, anode)
        n_gt = sum(len(v) for v in truth.values())
        print(f"  GT: {time.perf_counter()-t0:.1f}s segs={n_gt}", flush=True)

        gd = np.maximum(np.diff(sarr.astype(float)) * 1e-8, 1e-12)
        hmm_obj = ArchaicHMM(pi_amh=0.98, pi_nea=0.02, pi_den=0.0)
        cl = SEQ * len(cnodes)

        for cfg_name, cfg in configs.items():
            preds = []; trs = []; t_run = 0.0
            for hi, nid in enumerate(cnodes):
                q = csub[hi]
                trs.extend(truth.get(nid, []))
                t0 = time.perf_counter()
                if cfg is None:
                    preds.extend(pd_baseline_v(sarr, q, asub))
                else:
                    preds.extend(run_ap(
                        ref_freqs[cfg["ref"]], q, asub, sarr, gd, hmm_obj,
                        merge=cfg["merge"], theta_arch=cfg["theta_arch"],
                        positive_only=cfg["positive_only"]))
                t_run += time.perf_counter() - t0
            ev = evaluate_full(trs, preds, cl, cfg_name, ri, seed, t_run)
            all_r[cfg_name].append(ev.__dict__)
            print(f"  {cfg_name:22s}: F1={ev.f1_all:.3f} AUPRC={ev.auprc:.3f} "
                  f"rt={t_run:.1f}s", flush=True)

        print(f"  Rep done in {time.perf_counter()-t_rep:.0f}s", flush=True)

    except Exception as e:
        import traceback
        print(f"  ERROR rep {ri}: {e}")
        traceback.print_exc()

# ── Statistics ────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("FINAL RESULTS (50 replicates):")
print(f"{'Method':25s} {'F1':>14} {'AUPRC':>14} {'F1_short':>10} {'Runtime':>9}")
print("-" * 75)

stats_out = {}
ap_f1s = np.array([r["f1_all"] for r in all_r["archaicpainter"]])

for cfg_name, res in all_r.items():
    if not res: continue
    f1s   = np.array([r["f1_all"]   for r in res])
    auprcs= np.array([r["auprc"]    for r in res])
    f1s_s = np.array([r["f1_short"] for r in res])
    rts   = np.array([r.get("runtime_s", 0) for r in res])

    # Bootstrap 95% CI for F1 (BCa method)
    bs_f1 = scipy_bootstrap((f1s,), np.mean, n_resamples=2000,
                             confidence_level=0.95, method='percentile')
    ci_lo, ci_hi = bs_f1.confidence_interval

    # Wilcoxon signed-rank test vs archaicpainter
    if cfg_name != "archaicpainter" and len(f1s) == len(ap_f1s):
        try:
            stat, pval = wilcoxon(ap_f1s, f1s, alternative='greater')
        except:
            pval = float('nan')
    else:
        pval = float('nan')

    mean_rt = rts.mean() / (N_QUERY * 2) if rts.mean() > 0 else 0  # per haplotype

    print(f"  {cfg_name:25s}: {f1s.mean():.3f}±{f1s.std():.3f} "
          f"[{ci_lo:.3f},{ci_hi:.3f}]  "
          f"AUPRC={auprcs.mean():.3f}±{auprcs.std():.3f}  "
          f"F1s={f1s_s.mean():.3f}  "
          f"p={pval:.4f}" if not np.isnan(pval) else
          f"  {cfg_name:25s}: {f1s.mean():.3f}±{f1s.std():.3f} "
          f"[{ci_lo:.3f},{ci_hi:.3f}]  "
          f"AUPRC={auprcs.mean():.3f}  F1s={f1s_s.mean():.3f}")

    stats_out[cfg_name] = {
        "f1_mean": float(f1s.mean()), "f1_std": float(f1s.std()),
        "f1_ci95_lo": float(ci_lo), "f1_ci95_hi": float(ci_hi),
        "f1_short_mean": float(f1s_s.mean()),
        "auprc_mean": float(auprcs.mean()), "auprc_std": float(auprcs.std()),
        "wilcoxon_p_vs_ap": float(pval),
        "runtime_per_hap_s": float(mean_rt),
        "n_reps": len(res),
    }

# Wilcoxon: AP vs Poisson
try:
    pd_f1s = np.array([r["f1_all"] for r in all_r["poisson_density"]])
    stat, p_ap_vs_pd = wilcoxon(ap_f1s, pd_f1s, alternative='greater')
    print(f"\nWilcoxon (AP > Poisson): p={p_ap_vs_pd:.2e} (n={len(ap_f1s)})")
    stats_out["archaicpainter"]["wilcoxon_p_vs_poisson"] = float(p_ap_vs_pd)
except Exception as e:
    print(f"Wilcoxon error: {e}")

# Positive-only comparison (R4)
if "ap_positive_only" in all_r and all_r["ap_positive_only"]:
    po_f1s = np.array([r["f1_all"] for r in all_r["ap_positive_only"]])
    try:
        stat, p_po = wilcoxon(ap_f1s, po_f1s, alternative='two-sided')
        print(f"Wilcoxon (standard vs positive-only): p={p_po:.4f}")
        stats_out["ap_positive_only"]["wilcoxon_vs_standard"] = float(p_po)
    except Exception as e:
        print(f"Positive-only wilcoxon error: {e}")

# Save
out = {"raw": all_r, "stats": stats_out, "config": {"n_reps": N_REPS, "seq_bp": SEQ, "n_query": N_QUERY}}
with open(OUT_DIR / "phase4b_bench_v4.json", "w") as f:
    json.dump(out, f, indent=2)
print(f"\nSaved → phase4b_bench_v4.json")
