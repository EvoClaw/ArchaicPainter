"""
Phase 4b Large-Scale Simulation Benchmark
- 50 Mb x 20 seeds x 20 query haplotypes
- Ablation: YRI vs CEU reference, merging on/off, theta ablation
"""
import sys, time, json, warnings, numpy as np
from pathlib import Path
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
from archaicpainter.core.painter import ArchaicPainter
from archaicpainter.core.hmm import ArchaicHMM, decode_segments, STATE_NAMES
from archaicpainter.core.emission import compute_amh_emission, compute_archaic_emission

# ── Helpers ───────────────────────────────────────────────────────────────────
def extract_gt(ts, ceu_nodes, archaic_node, t_split=18966, gap=5000):
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
    pim = {v: i for i, v in enumerate(pos.tolist())}
    freq = np.full(len(shared), 0.5 / (n_hap + 1.0), dtype=np.float32)
    for i, p in enumerate(shared):
        if p in pim:
            j = pim[p]; freq[i] = ((hap[:, j] == 1).sum() + 0.5) / (n_hap + 1.0)
    return freq

def merge_nea(segs, gap=50_000, minlen=30_000):
    if not segs: return []
    segs = sorted(segs, key=lambda s: s.start)
    m = [Segment(segs[0].start, segs[0].end, segs[0].source, segs[0].posterior)]
    for s in segs[1:]:
        if s.start - m[-1].end <= gap:
            m[-1] = Segment(m[-1].start, max(m[-1].end, s.end), m[-1].source, max(m[-1].posterior, s.posterior))
        else: m.append(Segment(s.start, s.end, s.source, s.posterior))
    return [s for s in m if s.length >= minlen]

def run_ap(ref_freq, q, ag, pos, theta_arch=0.01, theta_amh=0.01):
    if len(pos) < 10: return []
    gd = np.maximum(np.diff(pos.astype(float)) * 1e-8, 1e-12)
    hmm = ArchaicHMM(pi_amh=0.98, pi_nea=0.02, pi_den=0.0)
    log_emit = np.zeros((3, len(pos)), dtype=np.float64)
    log_emit[0] = compute_amh_emission(q, ref_freq, theta=theta_amh)
    log_emit[1] = compute_archaic_emission(q, ag, theta=theta_arch)
    path = hmm.viterbi(log_emit, gd)
    posteriors, _ = hmm.forward_backward(log_emit, gd)
    segs = decode_segments(path, pos, posteriors, STATE_NAMES, min_length_bp=5_000)
    raw = [Segment(s["start"], s["end"], s["state"], s["mean_posterior"])
           for s in segs if s["state"] == "NEA"]
    return merge_nea(raw)

def pd_baseline(positions, q, ag, win=50000, step=10000):
    if len(positions) == 0: return []
    fixed = (ag[:,0]==1) & (ag[:,1]==1)
    hits = (q==1) & fixed
    clen = int(positions[-1]) + 1
    segs = []
    w = 0
    while w < clen:
        mask = (positions >= w) & (positions < w + win)
        n = mask.sum()
        if n > 0:
            f = hits[mask].sum() / n
            if f > 0: segs.append(Segment(w, w+win, "NEA", float(f)))
        w += step
    return segs

# ── Main Benchmark ────────────────────────────────────────────────────────────
SEEDS = list(range(42, 62))   # 20 seeds
N_REPS = 20; N_QUERY = 20; SEQ = 50_000_000
T_SPLIT = DEMOGRAPHY_PARAMS["T_split_nea"]

configs = {
    "archaicpainter":    {"ref": "YRI", "theta_arch": 0.01, "theta_amh": 0.01},
    "ap_ceu_ref":        {"ref": "CEU", "theta_arch": 0.01, "theta_amh": 0.01},  # ablation
    "ap_theta_high":     {"ref": "YRI", "theta_arch": 0.05, "theta_amh": 0.01},  # ablation
    "poisson_density":   None,
}

all_r = {k: [] for k in configs}

print(f"Large benchmark: {SEQ//1_000_000} Mb x {N_REPS} seeds x {N_QUERY} samples")
print("=" * 70)

for ri, seed in enumerate(SEEDS[:N_REPS]):
    print(f"\nRep {ri+1}/{N_REPS} seed={seed}", flush=True)
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
        nspec_count = sum(1 for p in shared if p not in set(yp.tolist()))
        print(f"  sites={len(shared)} NEA_spec={nspec_count}", flush=True)

        t0 = time.perf_counter()
        truth = extract_gt(ts, cnodes, anode, T_SPLIT)
        print(f"  GT:{time.perf_counter()-t0:.1f}s segs={sum(len(v) for v in truth.values())}", flush=True)

        cl = SEQ * len(cnodes)
        rep_results = {}

        for cfg_name, cfg in configs.items():
            preds = []; trs = []; t_run = 0.0
            for hi, nid in enumerate(cnodes):
                q = csub[hi]
                trs.extend(truth.get(nid, []))
                t0 = time.perf_counter()
                if cfg is None:
                    preds.extend(pd_baseline(sarr, q, asub))
                else:
                    preds.extend(run_ap(ref_freqs[cfg["ref"]], q, asub, sarr,
                                        theta_arch=cfg["theta_arch"],
                                        theta_amh=cfg["theta_amh"]))
                t_run += time.perf_counter() - t0
            ev = evaluate_full(trs, preds, cl, cfg_name, ri, seed, t_run)
            all_r[cfg_name].append(ev.__dict__)
            rep_results[cfg_name] = ev

        # Print summary line
        ap = rep_results["archaicpainter"]
        pd_ = rep_results["poisson_density"]
        print(f"  AP: F1={ap.f1_all:.3f} AUPRC={ap.auprc:.3f} FDR={ap.fdr_5pct:.3f} npred={ap.n_pred}")
        print(f"  PD: F1={pd_.f1_all:.3f} AUPRC={pd_.auprc:.3f}", flush=True)

    except Exception as e:
        import traceback
        print(f"  ERROR: {e}")
        traceback.print_exc()

print("\n" + "=" * 70)
print("FINAL RESULTS:")
for cfg_name, res in all_r.items():
    if not res: continue
    f1   = np.mean([r["f1_all"]   for r in res])
    f1s  = np.mean([r["f1_short"] for r in res])
    au   = np.mean([r["auprc"]    for r in res])
    fd   = np.mean([r["fdr_5pct"] for r in res])
    n_p  = np.mean([r["n_pred"]   for r in res])
    print(f"  {cfg_name:30s}: F1={f1:.3f} F1s={f1s:.3f} AUPRC={au:.3f} FDR={fd:.3f} npred={n_p:.0f}")

with open(OUT_DIR / "phase4b_benchmark_large.json", "w") as f:
    json.dump(all_r, f, indent=2)
print("Saved → phase4b_benchmark_large.json")
