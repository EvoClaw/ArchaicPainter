"""
Phase 4b Real Data Pipeline
Runs ArchaicPainter on 1000GP chr21 CEU samples against Vindija33.19 Neanderthal.

Data:
  - Query:    1000GP chr21 CEU haplotypes (phased, high-coverage)
  - Reference: 1000GP chr21 YRI samples (for AMH emission)
  - Archaic:  Vindija33.19 chr21 VCF (Prufer et al. 2017)

Output:
  - Per-sample predicted introgressed segments BED file
  - Summary statistics (ancestry fraction, segment length distribution)
"""
import sys, time, json, warnings, subprocess, io
import numpy as np
from pathlib import Path
warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))

from archaicpainter.core.emission import compute_amh_emission, compute_archaic_emission
from archaicpainter.core.hmm import ArchaicHMM, decode_segments, STATE_NAMES
from archaicpainter.evaluation.metrics import Segment

# ── Config ─────────────────────────────────────────────────────────────────────
VCF_1KGP = "/home/yanlin/livestock/1000GP/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
VCF_VINDIJA = "/home/yanlin/livestock/data/archaic/vindija33.19_chr21.vcf.gz"
SAMPLES_INFO = "/home/yanlin/livestock/1000GP/samples.info"
OUT_DIR = Path("/home/yanlin/livestock/docs/05_execution")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Chr21 region: skip telomeric + centromeric gaps
CHR21_REGION = "chr21:9000000-46709983"
VINDIJA_REGION = "21:9000000-46709983"  # Vindija uses "21" naming

N_CEU = 20    # number of CEU samples (= 40 haplotypes)
N_YRI = 50    # number of YRI samples for reference panel

# ── Load samples ──────────────────────────────────────────────────────────────
def load_samples(info_file, pop):
    samples = []
    with open(info_file) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2 and parts[1] == pop:
                samples.append(parts[0])
    return samples

print("Loading sample lists...", flush=True)
ceu_samples = load_samples(SAMPLES_INFO, "CEU")[:N_CEU]
yri_samples = load_samples(SAMPLES_INFO, "YRI")[:N_YRI]
print(f"  CEU: {len(ceu_samples)} samples = {2*len(ceu_samples)} haplotypes")
print(f"  YRI: {len(yri_samples)} samples = {2*len(yri_samples)} haplotypes")

# Write sample files
Path("/tmp/ceu_samples.txt").write_text("\n".join(ceu_samples))
Path("/tmp/yri_samples.txt").write_text("\n".join(yri_samples))

# ── Extract VCF data ──────────────────────────────────────────────────────────
def bcftools_query(vcf, region, samples_file, fields):
    """Run bcftools query and return lines."""
    cmd = ["bcftools", "query",
           "-r", region,
           "-s", open(samples_file).read().strip().replace("\n", ","),
           "-f", fields,
           vcf]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if result.returncode != 0:
        print(f"bcftools error: {result.stderr[:200]}", flush=True)
    return result.stdout.strip().split("\n")

def read_haplotype_matrix(vcf, region, samples_file, n_samples):
    """
    Returns:
      hap: (2*n_samples, n_sites) int8
      pos: (n_sites,) int32 (1-based positions)
    """
    print(f"  Reading VCF {Path(vcf).name} region={region} samples={n_samples}", flush=True)
    # bcftools: filter biallelic SNPs, get GT for target samples
    cmd = ["bcftools", "view", "-v", "snps", "-m2", "-M2",
           "-S", samples_file, "-r", region, vcf]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)

    positions = []
    rows = []
    for line in io.TextIOWrapper(proc.stdout):
        if line.startswith("#"): continue
        parts = line.rstrip().split("\t")
        chrom, pos, _, ref, alt = parts[:5]
        fmt = parts[8]
        gt_idx = fmt.split(":").index("GT")
        gts = []
        for sample_field in parts[9:]:
            gt_raw = sample_field.split(":")[gt_idx]
            alleles = gt_raw.replace("|", "/").split("/")
            try:
                gts.extend([int(a) for a in alleles])
            except:
                gts.extend([0, 0])
        positions.append(int(pos))
        rows.append(gts)

    proc.wait()
    if not rows:
        return np.zeros((2*n_samples, 0), dtype=np.int8), np.array([], dtype=np.int32)
    hap = np.array(rows, dtype=np.int8).T  # (2*n_samples, n_sites)
    pos = np.array(positions, dtype=np.int32)
    print(f"    → {hap.shape[1]} biallelic SNPs, {hap.shape[0]} haplotypes", flush=True)
    return hap, pos

def read_archaic_genotypes(vcf, region):
    """
    Returns:
      geno: (n_sites, 2) int8 (0/1 derived alleles)
      pos: (n_sites,) int32
    """
    print(f"  Reading archaic VCF {Path(vcf).name}", flush=True)
    cmd = ["bcftools", "view", "-v", "snps", "-m2", "-M2", "-r", region, vcf]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)

    positions = []
    genos = []
    for line in io.TextIOWrapper(proc.stdout):
        if line.startswith("#"): continue
        parts = line.rstrip().split("\t")
        fmt = parts[8]
        gt_idx = fmt.split(":").index("GT")
        gt_raw = parts[9].split(":")[gt_idx]
        alleles = gt_raw.replace("|", "/").split("/")
        try:
            g0, g1 = int(alleles[0]), int(alleles[1])
        except:
            continue
        if g0 > 0 or g1 > 0:  # only sites where archaic has derived allele
            positions.append(int(parts[1]))
            genos.append([g0, g1])

    proc.wait()
    if not genos:
        return np.zeros((0, 2), dtype=np.int8), np.array([], dtype=np.int32)
    geno = np.array(genos, dtype=np.int8)
    pos = np.array(positions, dtype=np.int32)
    print(f"    → {len(pos)} sites with archaic derived allele", flush=True)
    return geno, pos

# ── Run ────────────────────────────────────────────────────────────────────────
print("\n1. Extracting 1000GP chr21 haplotypes...", flush=True)
t0 = time.perf_counter()
ceu_hap, ceu_pos = read_haplotype_matrix(VCF_1KGP, CHR21_REGION, "/tmp/ceu_samples.txt", N_CEU)
yri_hap, yri_pos = read_haplotype_matrix(VCF_1KGP, CHR21_REGION, "/tmp/yri_samples.txt", N_YRI)
print(f"  Done in {time.perf_counter()-t0:.1f}s", flush=True)

print("\n2. Extracting Vindija chr21 derived alleles...", flush=True)
t0 = time.perf_counter()
vind_geno, vind_pos = read_archaic_genotypes(VCF_VINDIJA, VINDIJA_REGION)
print(f"  Done in {time.perf_counter()-t0:.1f}s", flush=True)

print("\n3. Building intersection...", flush=True)
shared = sorted(set(ceu_pos.tolist()) & set(vind_pos.tolist()))
sarr = np.array(shared, dtype=np.int32)
print(f"  Shared sites (CEU ∩ Vindija): {len(shared)}", flush=True)

if len(shared) < 100:
    print("ERROR: Too few shared sites. Check chromosome naming."); sys.exit(1)

# Index into matrices
cp_map = dict(zip(ceu_pos.tolist(), range(len(ceu_pos))))
vp_map = dict(zip(vind_pos.tolist(), range(len(vind_pos))))
yp_map = dict(zip(yri_pos.tolist(), range(len(yri_pos))))

ci = np.array([cp_map[x] for x in shared])
vi = np.array([vp_map[x] for x in shared])
ceu_sub = ceu_hap[:, ci]   # (n_ceu_hap, n_shared)
vind_sub = vind_geno[vi]   # (n_shared, 2)

# YRI frequencies at shared sites
n_yri_hap = yri_hap.shape[0]
ref_freq = np.full(len(shared), 0.5 / (n_yri_hap + 1.0), dtype=np.float32)
for k, pos in enumerate(shared):
    if pos in yp_map:
        j = yp_map[pos]
        ref_freq[k] = ((yri_hap[:, j] == 1).sum() + 0.5) / (n_yri_hap + 1.0)

yri_pos_set = set(yri_pos.tolist())
nea_spec = sum(1 for p in shared if p not in yri_pos_set)
print(f"  Vindija-specific (not in YRI): {nea_spec} ({100*nea_spec/len(shared):.1f}%)")

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

print("\n4. Running ArchaicPainter on CEU samples...", flush=True)
gd = np.maximum(np.diff(sarr.astype(float)) * 1e-8, 1e-12)
hmm = ArchaicHMM(pi_amh=0.98, pi_nea=0.02, pi_den=0.0)

all_segments = []
ancestry_fracs = []
t_total = 0.0

for hi in range(ceu_sub.shape[0]):
    sample_name = ceu_samples[hi // 2]
    hap_idx = hi % 2
    q = ceu_sub[hi]

    t0 = time.perf_counter()
    log_emit = np.zeros((3, len(sarr)), dtype=np.float64)
    log_emit[0] = compute_amh_emission(q, ref_freq)
    log_emit[1] = compute_archaic_emission(q, vind_sub)

    path = hmm.viterbi(log_emit, gd)
    posteriors, _ = hmm.forward_backward(log_emit, gd)
    t_total += time.perf_counter() - t0

    raw_segs = decode_segments(path, sarr, posteriors, STATE_NAMES, min_length_bp=5_000)
    nea_raw = [Segment(s["start"], s["end"], "NEA", s["mean_posterior"])
               for s in raw_segs if s["state"] == "NEA"]
    nea_merged = merge_nea(nea_raw)

    total_nea_bp = sum(s.length for s in nea_merged)
    chr21_len = 46709983 - 9000000
    af = total_nea_bp / chr21_len
    ancestry_fracs.append(af)
    all_segments.extend([(sample_name, hap_idx, s) for s in nea_merged])

    if hi % 10 == 0:
        print(f"  hi={hi}: NEA segs={len(nea_merged)} frac={af:.4f}", flush=True)

# ── Summary ────────────────────────────────────────────────────────────────────
n_hap = ceu_sub.shape[0]
mean_af = np.mean(ancestry_fracs)
seg_lens = [s.length for _, _, s in all_segments]

print(f"\n{'='*60}")
print(f"Results: {N_CEU} CEU individuals ({n_hap} haplotypes)")
print(f"  Total predicted NEA segments: {len(all_segments)}")
print(f"  Mean NEA ancestry fraction:   {mean_af:.4f} ({100*mean_af:.2f}%)")
if seg_lens:
    print(f"  Segment length: mean={np.mean(seg_lens)/1000:.1f}kb "
          f"median={np.median(seg_lens)/1000:.1f}kb "
          f"max={np.max(seg_lens)/1000:.1f}kb")
print(f"  Total runtime: {t_total:.1f}s ({t_total/n_hap:.2f}s/haplotype)")

# ── Save BED ──────────────────────────────────────────────────────────────────
bed_path = OUT_DIR / "archaicpainter_chr21_CEU.bed"
with open(bed_path, "w") as f:
    f.write("#chrom\tstart\tend\tsample\thap\tstate\tposterior\tlength_kb\n")
    for samp, hi, seg in sorted(all_segments, key=lambda x: (x[0], x[1], x[2].start)):
        f.write(f"chr21\t{seg.start}\t{seg.end}\t{samp}\t{hi}\tNEA\t"
                f"{seg.posterior:.4f}\t{seg.length/1000:.2f}\n")
print(f"\nSaved BED → {bed_path}")

# Save JSON summary
summary = {
    "n_individuals": N_CEU,
    "n_haplotypes": n_hap,
    "n_segments": len(all_segments),
    "mean_ancestry_fraction": float(mean_af),
    "mean_segment_length_kb": float(np.mean(seg_lens)/1000) if seg_lens else 0,
    "median_segment_length_kb": float(np.median(seg_lens)/1000) if seg_lens else 0,
    "per_haplotype_ancestry": [float(x) for x in ancestry_fracs],
}
with open(OUT_DIR / "realdata_summary.json", "w") as f:
    json.dump(summary, f, indent=2)
print(f"Saved JSON → realdata_summary.json")
