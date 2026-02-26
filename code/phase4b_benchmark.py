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
from archaicpainter.core.hmm import ArchaicHMM

def extract_gt(ts, ceu_nodes, archaic_node, t_split=18966, gap=5000):
    per = {n: [] for n in ceu_nodes}
    for tree in ts.trees():
        iv = tree.interval
        if iv.right <= iv.left:
            continue
        for cn in ceu_nodes:
            try:
                m = tree.mrca(cn, archaic_node)
                if ts.node(m).time < t_split * 0.5:
                    per[cn].append((iv.left, iv.right))
            except:
                pass
    result = {}
    for cn, ivs in per.items():
        if not ivs:
            result[cn] = []
            continue
        s = sorted(ivs)
        merged = [list(s[0])]
        for st, en in s[1:]:
            if st - merged[-1][1] <= gap:
                merged[-1][1] = max(merged[-1][1], en)
            else:
                merged.append([st, en])
        result[cn] = [Segment(int(a), int(b), "NEA", 1.0) for a, b in merged]
    return result

def make_yri_freq(shared, yh, yp, n_yri):
    yim = {v: i for i, v in enumerate(yp.tolist())}
    freq = np.full(len(shared), 0.5 / (n_yri + 1.0), dtype=np.float32)
    for i, pos in enumerate(shared):
        if pos in yim:
            j = yim[pos]
            n_alt = (yh[:, j] == 1).sum()
            freq[i] = (n_alt + 0.5) / (n_yri + 1.0)
    return freq

def pd_baseline(positions, q, ag, win=50000, step=10000):
    if len(positions) == 0:
        return []
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
            if f > 0:
                segs.append(Segment(w, w+win, "NEA", float(f)))
        w += step
    return segs

def merge_nea_segments(segs, merge_gap_bp=50_000, min_len_bp=30_000):
    """Merge NEA fragments within merge_gap_bp, then filter by min_len_bp."""
    if not segs:
        return []
    segs = sorted(segs, key=lambda s: s.start)
    merged = [Segment(segs[0].start, segs[0].end, segs[0].source, segs[0].posterior)]
    for s in segs[1:]:
        if s.start - merged[-1].end <= merge_gap_bp:
            merged[-1] = Segment(merged[-1].start, max(merged[-1].end, s.end),
                                  merged[-1].source, max(merged[-1].posterior, s.posterior))
        else:
            merged.append(Segment(s.start, s.end, s.source, s.posterior))
    return [s for s in merged if s.length >= min_len_bp]

def run_ap(ref_freq, q, ag, pos):
    if len(pos) < 10:
        return []
    from archaicpainter.core.emission import compute_amh_emission, compute_archaic_emission
    from archaicpainter.core.hmm import ArchaicHMM, decode_segments, STATE_NAMES
    gd = np.maximum(np.diff(pos.astype(float)) * 1e-8, 1e-12)
    hmm = ArchaicHMM(pi_amh=0.98, pi_nea=0.02, pi_den=0.0)
    log_emit = np.zeros((3, len(pos)), dtype=np.float64)
    log_emit[0] = compute_amh_emission(q, ref_freq)
    log_emit[1] = compute_archaic_emission(q, ag)
    posteriors, _ = hmm.forward_backward(log_emit, gd)
    path = hmm.viterbi(log_emit, gd)
    segs = decode_segments(path, pos, posteriors, STATE_NAMES, min_length_bp=5_000)
    raw = [Segment(s["start"], s["end"], s["state"], s["mean_posterior"])
           for s in segs if s["state"] == "NEA"]
    # Merge fragments within 50 kb gap, then keep >= 30 kb
    return merge_nea_segments(raw, merge_gap_bp=50_000, min_len_bp=30_000)

SEEDS = [42, 123, 456, 789, 1024]
N_REPS = 5; N_QUERY = 10; SEQ = 5_000_000
T_SPLIT = DEMOGRAPHY_PARAMS["T_split_nea"]
all_r = {"archaicpainter": [], "poisson_density": []}

print("Phase 4b Benchmark: 2-way intersection + YRI freq AMH emission")
print("=" * 70)

for ri, seed in enumerate(SEEDS[:N_REPS]):
    print(f"Rep {ri+1}/{N_REPS} seed={seed}", flush=True)
    try:
        ts, ts_mut = simulate_one(seed=seed, n_query=N_QUERY, n_archaic=1, seq_len=SEQ)
        pid = {p.metadata.get("name", str(p.id)): p.id for p in ts.populations()}
        cnodes = list(ts.samples(population=pid["CEU"]))
        anode = list(ts.samples(population=pid["NEA"]))[0]
        ch, cp = ts_to_haplotype_matrix(ts_mut, "CEU")
        yh, yp = ts_to_haplotype_matrix(ts_mut, "YRI")
        ah, ap2 = ts_to_archaic_genotypes(ts_mut, "NEA")
        shared = sorted(set(cp.tolist()) & set(ap2.tolist()))
        sarr = np.array(shared, dtype=np.int32)
        ci = np.array([dict(zip(cp.tolist(), range(len(cp))))[x] for x in shared])
        ai = np.array([dict(zip(ap2.tolist(), range(len(ap2))))[x] for x in shared])
        csub = ch[:, ci]; asub = ah[ai]
        ref_freq = make_yri_freq(shared, yh, yp, yh.shape[0])
        in_yri = sum(1 for x in shared if x in set(yp.tolist()))
        nea_spec = len(shared) - in_yri
        print(f"  sites={len(shared)} in_YRI={in_yri} NEA_spec={nea_spec}", end=" ", flush=True)
        if ri == 0:
            from archaicpainter.core.emission import compute_amh_emission, compute_archaic_emission
            q0 = csub[0]
            la = compute_amh_emission(q0, ref_freq)
            ln = compute_archaic_emission(q0, asub)
            hom = (asub[:,0]==1) & (asub[:,1]==1)
            nspec = np.array([x not in set(yp.tolist()) for x in shared])
            print()
            print(f"  DIAG: arch_hom={hom.sum()} nea_spec={nspec.sum()}")
            print(f"        NEA-AMH@all={((ln-la).mean()):.3f}")
            print(f"        NEA-AMH@nea_spec={((ln-la)[nspec].mean() if nspec.sum()>0 else 0):.3f}")
            print(f"        NEA-AMH@hom_nspec={((ln-la)[hom&nspec].mean() if (hom&nspec).sum()>0 else 0):.3f} n={(hom&nspec).sum()}")
        else:
            print()
        t0 = time.perf_counter()
        truth = extract_gt(ts, cnodes, anode, T_SPLIT)
        print(f"  GT:{time.perf_counter()-t0:.1f}s segs={sum(len(v) for v in truth.values())}", flush=True)
        aps = []; pds = []; trs = []
        at = pt = 0.0
        for hi, nid in enumerate(cnodes):
            q = csub[hi]
            trs.extend(truth.get(nid, []))
            t0 = time.perf_counter(); aps.extend(run_ap(ref_freq, q, asub, sarr)); at += time.perf_counter()-t0
            t0 = time.perf_counter(); pds.extend(pd_baseline(sarr, q, asub)); pt += time.perf_counter()-t0
        cl = SEQ * len(cnodes)
        ae = evaluate_full(trs, aps, cl, "archaicpainter", ri, seed, at)
        pe = evaluate_full(trs, pds, cl, "poisson_density", ri, seed, pt)
        all_r["archaicpainter"].append(ae.__dict__)
        all_r["poisson_density"].append(pe.__dict__)
        print(f"  AP: F1={ae.f1_all:.3f} F1s={ae.f1_short:.3f} AUPRC={ae.auprc:.3f} npred={ae.n_pred} t={at:.1f}s")
        print(f"  PD: F1={pe.f1_all:.3f} F1s={pe.f1_short:.3f} npred={pe.n_pred} ntrue={pe.n_true}")
    except Exception as e:
        import traceback; print(f"  ERROR: {e}"); traceback.print_exc()

print("=" * 70)
for meth, res in all_r.items():
    if not res: continue
    f1 = np.mean([r["f1_all"] for r in res])
    f1s = np.mean([r["f1_short"] for r in res])
    au = np.mean([r["auprc"] for r in res])
    fd = np.mean([r["fdr_5pct"] for r in res])
    print(f"{meth}: F1={f1:.3f} F1s={f1s:.3f} AUPRC={au:.3f} FDR={fd:.3f}")
with open(OUT_DIR/"phase4b_benchmark_results.json","w") as f:
    json.dump(all_r, f, indent=2)
print("Done")
