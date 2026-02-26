"""
Phase 4b Real Data Pipeline v2 - with hg19→hg38 liftover for Vindija
"""
import sys, time, json, warnings, subprocess, io
import numpy as np
from pathlib import Path
warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))

from archaicpainter.core.emission import compute_amh_emission, compute_archaic_emission
from archaicpainter.core.hmm import ArchaicHMM, decode_segments, STATE_NAMES
from archaicpainter.evaluation.metrics import Segment

VCF_1KGP = "/home/yanlin/livestock/1000GP/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
VCF_VINDIJA = "/home/yanlin/livestock/data/archaic/vindija33.19_chr21.vcf.gz"
CHAIN_FILE = "/home/yanlin/livestock/data/archaic/hg19ToHg38.over.chain.gz"
SAMPLES_INFO = "/home/yanlin/livestock/1000GP/samples.info"
CROSSMAP = "/home/yanlin/miniconda3/bin/CrossMap"
OUT_DIR = Path("/home/yanlin/livestock/docs/05_execution")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# use FULL chr21 in hg38 (after telomere gap)
CHR21_REGION = "chr21:5010000-46709983"  # hg38 chr21

N_CEU = 20    # 20 individuals = 40 haplotypes
N_YRI = 50    # 50 individuals = 100 haplotypes

# ── Sample loading ─────────────────────────────────────────────────────────────
def load_samples(info_file, pop, n=None):
    samples = []
    with open(info_file) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2 and parts[1] == pop:
                samples.append(parts[0])
    return samples[:n] if n else samples

print("Loading samples...", flush=True)
ceu_samples = load_samples(SAMPLES_INFO, "CEU", N_CEU)
yri_samples = load_samples(SAMPLES_INFO, "YRI", N_YRI)
print(f"  CEU: {len(ceu_samples)} x YRI: {len(yri_samples)}", flush=True)

Path("/tmp/ceu_s.txt").write_text("\n".join(ceu_samples))
Path("/tmp/yri_s.txt").write_text("\n".join(yri_samples))

# ── 1000GP haplotype matrix ────────────────────────────────────────────────────
def read_1kgp_hap(samples_file, n_samples):
    print(f"  Parsing 1kGP VCF for {n_samples} samples...", flush=True)
    cmd = ["bcftools", "view", "-v", "snps", "-m2", "-M2",
           "-S", samples_file, "-r", CHR21_REGION, VCF_1KGP]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    positions, rows = [], []
    for line in io.TextIOWrapper(proc.stdout):
        if line.startswith("#"): continue
        parts = line.rstrip().split("\t")
        pos = int(parts[1])
        fmt = parts[8]; gt_idx = fmt.split(":").index("GT")
        gts = []
        for sf in parts[9:]:
            gt_raw = sf.split(":")[gt_idx]
            a = gt_raw.replace("|", "/").split("/")
            try: gts.extend([int(x) for x in a])
            except: gts.extend([0, 0])
        positions.append(pos)
        rows.append(gts)
    proc.wait()
    if not rows:
        return np.zeros((2*n_samples, 0), dtype=np.int8), np.array([], dtype=np.int32)
    hap = np.array(rows, dtype=np.int8).T
    pos = np.array(positions, dtype=np.int32)
    print(f"    → {hap.shape[1]} SNPs, {hap.shape[0]} haplotypes", flush=True)
    return hap, pos

# ── Vindija liftover ───────────────────────────────────────────────────────────
def load_vindija_lifted():
    """Extract Vindija variants, liftover hg19→hg38, return arrays."""
    print("  Extracting Vindija variants...", flush=True)
    # Step 1: extract to BED
    cmd = ["bcftools", "query", "-r", "21",
           "-f", "21\t%POS0\t%POS\t%REF\t%ALT\t[%GT]\n", VCF_VINDIJA]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    lines = [l for l in r.stdout.strip().split("\n")
             if l and l.split("\t")[4] not in ("", ".")]
    Path("/tmp/vindija_hg19.bed").write_text("\n".join(lines))
    print(f"    → {len(lines)} Vindija variants extracted", flush=True)

    # Step 2: CrossMap liftover
    print("  CrossMap hg19 → hg38...", flush=True)
    subprocess.run([CROSSMAP, "bed", CHAIN_FILE,
                    "/tmp/vindija_hg19.bed", "/tmp/vindija_hg38.bed"],
                   capture_output=True, timeout=120)

    # Step 3: Parse lifted BED
    positions, genos = [], []
    with open("/tmp/vindija_hg38.bed") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 6: continue
            pos = int(parts[2])  # 1-based
            gt_raw = parts[5]
            alleles = gt_raw.replace("|", "/").split("/")
            try:
                g0, g1 = int(alleles[0]), int(alleles[1])
            except:
                continue
            if g0 > 0 or g1 > 0:
                positions.append(pos)
                genos.append([g0, g1])
    print(f"    → {len(positions)} Vindija sites with derived allele (hg38 coords)", flush=True)
    return np.array(genos, dtype=np.int8), np.array(positions, dtype=np.int32)

# ── Main ──────────────────────────────────────────────────────────────────────
print("\n1. Reading 1000GP haplotypes...", flush=True)
t0 = time.perf_counter()
ceu_hap, ceu_pos = read_1kgp_hap("/tmp/ceu_s.txt", N_CEU)
yri_hap, yri_pos = read_1kgp_hap("/tmp/yri_s.txt", N_YRI)
print(f"  Done in {time.perf_counter()-t0:.1f}s", flush=True)

print("\n2. Loading Vindija genotypes (with hg19→hg38 liftover)...", flush=True)
t0 = time.perf_counter()
vind_geno, vind_pos = load_vindija_lifted()
print(f"  Done in {time.perf_counter()-t0:.1f}s", flush=True)

print("\n3. Building intersection...", flush=True)
shared = sorted(set(ceu_pos.tolist()) & set(vind_pos.tolist()))
sarr = np.array(shared, dtype=np.int32)
print(f"  Shared sites (CEU ∩ Vindija_hg38): {len(shared)}", flush=True)

if len(shared) < 100:
    print("ERROR: Too few shared sites. Dumping first 10 Vindija positions...")
    print("  Vindija:", vind_pos[:10].tolist())
    print("  CEU:    ", ceu_pos[:10].tolist())
    sys.exit(1)

cp_map = dict(zip(ceu_pos.tolist(), range(len(ceu_pos))))
vp_map = dict(zip(vind_pos.tolist(), range(len(vind_pos))))
yp_map = dict(zip(yri_pos.tolist(), range(len(yri_pos))))

ci = np.array([cp_map[x] for x in shared])
vi = np.array([vp_map[x] for x in shared])
ceu_sub = ceu_hap[:, ci]
vind_sub = vind_geno[vi]

# YRI reference frequencies
n_yri_hap = yri_hap.shape[0]
ref_freq = np.full(len(shared), 0.5 / (n_yri_hap + 1.0), dtype=np.float32)
for k, pos in enumerate(shared):
    if pos in yp_map:
        j = yp_map[pos]
        ref_freq[k] = ((yri_hap[:, j] == 1).sum() + 0.5) / (n_yri_hap + 1.0)

nea_spec = sum(1 for p in shared if p not in yp_map)
mean_vf = ref_freq.mean()
print(f"  Vindija-specific (absent from YRI): {nea_spec} ({100*nea_spec/len(shared):.1f}%)")
print(f"  Mean YRI allele freq at shared sites: {mean_vf:.4f}")

# ── ArchaicPainter ────────────────────────────────────────────────────────────
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

print("\n4. Running ArchaicPainter on CEU chr21...", flush=True)
gd = np.maximum(np.diff(sarr.astype(float)) * 1e-8, 1e-12)
hmm = ArchaicHMM(pi_amh=0.98, pi_nea=0.02, pi_den=0.0)

all_segments = []
ancestry_fracs = []
t_total = 0.0

for hi in range(ceu_sub.shape[0]):
    q = ceu_sub[hi]
    t0 = time.perf_counter()
    log_emit = np.zeros((3, len(sarr)), dtype=np.float64)
    log_emit[0] = compute_amh_emission(q, ref_freq)
    log_emit[1] = compute_archaic_emission(q, vind_sub)
    path = hmm.viterbi(log_emit, gd)
    posteriors, _ = hmm.forward_backward(log_emit, gd)
    t_total += time.perf_counter() - t0

    raw = decode_segments(path, sarr, posteriors, STATE_NAMES, min_length_bp=5_000)
    nea_raw = [Segment(s["start"], s["end"], "NEA", s["mean_posterior"])
               for s in raw if s["state"] == "NEA"]
    nea = merge_nea(nea_raw)
    total_bp = sum(s.length for s in nea)
    af = total_bp / (sarr[-1] - sarr[0]) if len(sarr) > 1 else 0
    ancestry_fracs.append(af)
    all_segments.extend([(ceu_samples[hi//2], hi%2, s) for s in nea])
    if hi % 8 == 0:
        print(f"  hap {hi}/{ceu_sub.shape[0]}: segs={len(nea)} AF={af:.4f}", flush=True)

# ── Summary ───────────────────────────────────────────────────────────────────
mean_af = np.mean(ancestry_fracs)
seg_lens = [s.length for _, _, s in all_segments]
print(f"\n{'='*60}")
print(f"  Haplotypes: {ceu_sub.shape[0]}")
print(f"  Total NEA segments: {len(all_segments)}")
print(f"  Mean NEA ancestry fraction: {mean_af:.4f} ({100*mean_af:.2f}%)")
if seg_lens:
    print(f"  Segment length: mean={np.mean(seg_lens)/1e3:.1f}kb "
          f"median={np.median(seg_lens)/1e3:.1f}kb "
          f"max={np.max(seg_lens)/1e3:.1f}kb")
    print(f"  n segs per haplotype: {len(all_segments)/ceu_sub.shape[0]:.1f}")
print(f"  Runtime: {t_total:.1f}s ({t_total/ceu_sub.shape[0]:.2f}s/hap)", flush=True)

# ── Save outputs ──────────────────────────────────────────────────────────────
bed_path = OUT_DIR / "archaicpainter_chr21_CEU_v2.bed"
with open(bed_path, "w") as f:
    f.write("#chrom\tstart\tend\tsample\thap\tstate\tposterior\tlength_kb\n")
    for samp, hi, seg in sorted(all_segments, key=lambda x: (x[0], x[1], x[2].start)):
        f.write(f"chr21\t{seg.start}\t{seg.end}\t{samp}\t{hi}\tNEA\t"
                f"{seg.posterior:.4f}\t{seg.length/1000:.2f}\n")
print(f"Saved BED → {bed_path}")

summary = {
    "n_individuals": N_CEU, "n_haplotypes": ceu_sub.shape[0],
    "shared_sites": len(shared), "nea_specific_sites": nea_spec,
    "n_segments": len(all_segments),
    "mean_ancestry_fraction": float(mean_af),
    "mean_segment_length_kb": float(np.mean(seg_lens)/1e3) if seg_lens else 0,
    "median_segment_length_kb": float(np.median(seg_lens)/1e3) if seg_lens else 0,
    "per_haplotype_ancestry": [float(x) for x in ancestry_fracs],
    "runtime_total_s": t_total,
}
with open(OUT_DIR / "realdata_summary_v2.json", "w") as f:
    json.dump(summary, f, indent=2)
print(f"Saved JSON → realdata_summary_v2.json")
