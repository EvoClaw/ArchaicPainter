"""
S4: Scalability Test - Runtime vs n_haplotypes
Tests ArchaicPainter runtime scaling: 10 to 1000 haplotypes
on a fixed 5Mb region (seed=42).
"""
import sys, time, json, warnings
import numpy as np
from pathlib import Path
warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))

import msprime
from archaicpainter.evaluation.simulation import (
    build_demography, simulate_one, DEMOGRAPHY_PARAMS,
    ts_to_haplotype_matrix, ts_to_archaic_genotypes,
)
from archaicpainter.core.hmm import ArchaicHMM, decode_segments, STATE_NAMES
from archaicpainter.core.emission import compute_amh_emission, compute_archaic_emission

OUT_DIR = Path("/home/yanlin/livestock/docs/05_execution")
OUT_DIR.mkdir(parents=True, exist_ok=True)
SEQ = 5_000_000
SEED = 42
N_HAP_LIST = [10, 20, 50, 100, 200, 500, 1000]


def make_freq(shared, hap, pos, n_hap):
    pos_map = dict(zip(pos.tolist(), range(len(pos))))
    freq = np.full(len(shared), 0.5 / (n_hap + 1.0), dtype=np.float32)
    for i, p in enumerate(shared):
        if p in pos_map:
            freq[i] = ((hap[:, pos_map[p]] == 1).sum() + 0.5) / (n_hap + 1.0)
    return freq


print("Scalability test: {}Mb, seed={}".format(SEQ // 1000000, SEED))
print("n_hap: {}".format(N_HAP_LIST))
print("=" * 60, flush=True)

results = []

for n_hap in N_HAP_LIST:
    n_query = n_hap // 2
    print("\nn_hap={} (n_query={})".format(n_hap, n_query), flush=True)
    try:
        ts, ts_mut = simulate_one(seed=SEED, n_query=n_query, n_archaic=1, seq_len=SEQ)

        ch, cp = ts_to_haplotype_matrix(ts_mut, "CEU")
        yh, yp = ts_to_haplotype_matrix(ts_mut, "YRI")
        ah, ap2 = ts_to_archaic_genotypes(ts_mut, "NEA")

        shared = sorted(set(cp.tolist()) & set(ap2.tolist()))
        sarr = np.array(shared, dtype=np.int32)
        cp_idx = dict(zip(cp.tolist(), range(len(cp))))
        ap_idx = dict(zip(ap2.tolist(), range(len(ap2))))
        ci = np.array([cp_idx[x] for x in shared])
        ai = np.array([ap_idx[x] for x in shared])
        csub = ch[:, ci]
        asub = ah[ai]
        ref_freq = make_freq(shared, yh, yp, yh.shape[0])
        gd = np.maximum(np.diff(sarr.astype(float)) * 1e-8, 1e-12)
        hmm_obj = ArchaicHMM(pi_amh=0.98, pi_nea=0.02, pi_den=0.0)

        t0 = time.perf_counter()
        for hi in range(ch.shape[0]):
            q = csub[hi]
            log_emit = np.zeros((3, len(sarr)), dtype=np.float64)
            log_emit[0] = compute_amh_emission(q, ref_freq)
            log_emit[1] = compute_archaic_emission(q, asub)
            path = hmm_obj.viterbi(log_emit, gd)
            _, _ = hmm_obj.forward_backward(log_emit, gd)
        t_total = time.perf_counter() - t0
        t_per_hap = t_total / ch.shape[0]

        print("  n_hap={}: total={:.2f}s  per_hap={:.4f}s  n_sites={}".format(
            n_hap, t_total, t_per_hap, len(sarr)), flush=True)

        results.append({
            "n_haplotypes": n_hap,
            "n_sites": len(sarr),
            "total_time_s": t_total,
            "per_hap_s": t_per_hap,
            "seq_mb": SEQ // 1000000,
        })

    except Exception as e:
        import traceback
        print("  ERROR: {}".format(e))
        traceback.print_exc()

print("\n" + "=" * 60)
print("SCALABILITY SUMMARY:")
print("{:>8} {:>12} {:>10} {:>8}".format("n_hap", "t_total(s)", "t/hap(s)", "n_sites"))
for r in results:
    print("  {:>6}  {:>10.2f}  {:>8.4f}  {:>8}".format(
        r["n_haplotypes"], r["total_time_s"], r["per_hap_s"], r["n_sites"]))

if results:
    r_last = results[-1]
    genome_scale = 3000000000 / SEQ
    t_genome_est = r_last["per_hap_s"] * r_last["n_haplotypes"] * genome_scale / 3600
    print("\nExtrapolated: {} hap x full genome (~3Gb): {:.1f} cpu-hours".format(
        r_last["n_haplotypes"], t_genome_est))

    nhaps = np.array([r["n_haplotypes"] for r in results], dtype=float)
    times = np.array([r["total_time_s"] for r in results], dtype=float)
    slope = np.polyfit(np.log(nhaps), np.log(times), 1)[0]
    print("Scaling exponent (log-log slope): {:.3f}  (1.0 = linear O(n))".format(slope))

with open(str(OUT_DIR / "phase4b_scalability.json"), "w") as f:
    json.dump({"results": results, "seq_bp": SEQ, "seed": SEED}, f, indent=2)
print("\nSaved -> phase4b_scalability.json")
