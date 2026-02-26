"""
Phase 4b Benchmark v3 - Optimized: 20Mb × 12 seeds × 10 samples
  - Vectorized pd_baseline
  - Optimized evaluate_full
  - Ablation: YRI vs CEU reference, segment merging, theta
"""
import sys, time, json, warnings
import numpy as np
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
from archaicpainter.core.hmm import ArchaicHMM, decode_segments, STATE_NAMES
from archaicpainter.core.emission import compute_amh_emission, compute_archaic_emission

SEEDS    = list(range(42, 54))  # 12 seeds
N_REPS   = 12
N_QUERY  = 10   # 10 samples = 20 haplotypes
SEQ      = 5_000_000
T_SPLIT  = DEMOGRAPHY_PARAMS["T_split_nea"]

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

# ── Frequency computation ──────────────────────────────────────────────────────
def make_freq(shared, hap, pos, n_hap):
    pos_map = dict(zip(pos.tolist(), range(len(pos))))
    freq = np.full(len(shared), 0.5 / (n_hap + 1.0), dtype=np.float32)
    for i, p in enumerate(shared):
        if p in pos_map:
            j = pos_map[p]
            freq[i] = ((hap[:, j] == 1).sum() + 0.5) / (n_hap + 1.0)
    return freq

# ── Segment merging ────────────────────────────────────────────────────────────
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

# ── ArchaicPainter ─────────────────────────────────────────────────────────────
def run_ap(ref_freq, q, ag, pos, gd, hmm_obj, merge=True, theta_arch=0.01, theta_amh=0.01):
    if len(pos) < 10: return []
    log_emit = np.zeros((3, len(pos)), dtype=np.float64)
    log_emit[0] = compute_amh_emission(q, ref_freq, theta=theta_amh)
    log_emit[1] = compute_archaic_emission(q, ag, theta=theta_arch)
    path = hmm_obj.viterbi(log_emit, gd)
    posteriors, _ = hmm_obj.forward_backward(log_emit, gd)
    segs = decode_segments(path, pos, posteriors, STATE_NAMES, min_length_bp=5_000)
    raw = [Segment(s["start"], s["end"], s["state"], s["mean_posterior"])
           for s in segs if s["state"] == "NEA"]
    return merge_nea(raw) if merge else [s for s in raw if s.length >= 30_000]

# ── Poisson density baseline (vectorized) ─────────────────────────────────────
def pd_baseline_v(positions, q, ag, win=50_000, step=20_000):
    if len(positions) == 0: return []
    fixed = (ag[:, 0] == 1) & (ag[:, 1] == 1)
    hits = (q == 1) & fixed
    if hits.sum() == 0: return []

    # Sort-based approach: for each window [w, w+win), count hits
    pos_f = positions.astype(np.float64)
    start = int(positions[0])
    end   = int(positions[-1]) + 1
    windows = np.arange(start, end, step, dtype=np.int64)

    # Use cumsum for vectorized window counting
    cum_hits = np.cumsum(np.concatenate([[0], hits.astype(np.int32)]))
    cum_tot  = np.cumsum(np.concatenate([[0], np.ones(len(positions), dtype=np.int32)]))

    def count_in_window(w_start, w_end):
        i_lo = np.searchsorted(positions, w_start, 'left')
        i_hi = np.searchsorted(positions, w_end, 'left')
        return cum_hits[i_hi] - cum_hits[i_lo], cum_tot[i_hi] - cum_tot[i_lo]

    segs = []
    for w in windows:
        h, n = count_in_window(w, w + win)
        if n > 0 and h > 0:
            score = h / n
            segs.append(Segment(int(w), int(w + win), "NEA", float(score)))
    return segs

# ── Configs ────────────────────────────────────────────────────────────────────
configs = {
    "archaicpainter":      {"ref": "YRI", "merge": True,  "theta_arch": 0.01},
    "ap_no_merge":         {"ref": "YRI", "merge": False, "theta_arch": 0.01},   # ablation
    "ap_ceu_ref":          {"ref": "CEU", "merge": True,  "theta_arch": 0.01},   # ablation
    "ap_theta_high":       {"ref": "YRI", "merge": True,  "theta_arch": 0.05},   # ablation
    "poisson_density":     None,
}

all_r = {k: [] for k in configs}

print(f"Benchmark v3: {SEQ//1_000_000} Mb × {N_REPS} seeds × {N_QUERY} samples/rep")
print(f"Configs: {list(configs.keys())}")
print("=" * 70, flush=True)

for ri, seed in enumerate(SEEDS[:N_REPS]):
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
        gt_time = time.perf_counter() - t0
        n_gt = sum(len(v) for v in truth.values())
        print(f"  GT:{gt_time:.1f}s segs={n_gt}", flush=True)

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
                    preds.extend(run_ap(ref_freqs[cfg["ref"]], q, asub, sarr, gd,
                                        hmm_obj, merge=cfg["merge"],
                                        theta_arch=cfg["theta_arch"]))
                t_run += time.perf_counter() - t0
            ev = evaluate_full(trs, preds, cl, cfg_name, ri, seed, t_run)
            all_r[cfg_name].append(ev.__dict__)
            print(f"  {cfg_name:20s}: F1={ev.f1_all:.3f} AUPRC={ev.auprc:.3f} "
                  f"FDR={ev.fdr_5pct:.3f} npred={ev.n_pred} rt={t_run:.1f}s",
                  flush=True)

        print(f"  Rep done in {time.perf_counter()-t_rep:.0f}s", flush=True)

    except Exception as e:
        import traceback
        print(f"  ERROR: {e}")
        traceback.print_exc()

print("\n" + "=" * 70)
print("FINAL RESULTS:")
for cfg_name, res in all_r.items():
    if not res: continue
    vals = {k: np.mean([r[k] for r in res]) for k in ("f1_all","f1_short","auprc","fdr_5pct","n_pred")}
    stds = {k: np.std([r[k] for r in res])  for k in ("f1_all","auprc")}
    print(f"  {cfg_name:22s}: F1={vals['f1_all']:.3f}±{stds['f1_all']:.3f} "
          f"F1s={vals['f1_short']:.3f} AUPRC={vals['auprc']:.3f}±{stds['auprc']:.3f} "
          f"FDR={vals['fdr_5pct']:.3f} npred={vals['n_pred']:.0f}")

with open(OUT_DIR / "phase4b_bench_v3.json", "w") as f:
    json.dump(all_r, f, indent=2)
print("\nSaved → phase4b_bench_v3.json")
