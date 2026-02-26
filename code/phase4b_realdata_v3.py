"""
Phase 4b Real Data Pipeline v3 - Allele-Aware, Positive-Only Emission
Uses chr21:10-20Mb (euchromatic, good coverage in both VCFs).
Fixes allele coding mismatch by checking REF/ALT at each shared position.
Uses positive-only archaic emission: mismatches are uninformative (not penalized)
because Vindija may differ from the actual introgressing Neanderthal population.
"""
import sys, time, json, warnings, subprocess
import numpy as np
from pathlib import Path
warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))

from archaicpainter.core.emission import compute_amh_emission, compute_archaic_emission
from archaicpainter.core.hmm import ArchaicHMM, decode_segments, STATE_NAMES
from archaicpainter.evaluation.metrics import Segment

VCF_1KGP  = "/home/yanlin/livestock/1000GP/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
VCF_VIND  = "/home/yanlin/livestock/data/archaic/vindija33.19_chr21.vcf.gz"
CHAIN     = "/home/yanlin/livestock/data/archaic/hg19ToHg38.over.chain.gz"
CROSSMAP  = "/home/yanlin/miniconda3/bin/CrossMap"
SAMPLES   = "/home/yanlin/livestock/1000GP/samples.info"
OUT_DIR   = Path("/home/yanlin/livestock/docs/05_execution")

REGION_HG38 = "chr21:10000000-20000000"
REGION_HG19 = "21:10000000-21000000"  # approximate hg19 range for 10-20Mb hg38 region

N_CEU = 20
N_YRI = 50

COMP = {"A":"T","T":"A","C":"G","G":"C"}

def load_samples(f, pop, n=None):
    s = [l.split()[0] for l in open(f) if len(l.split())>=2 and l.split()[1]==pop]
    return s[:n] if n else s

print("=== Phase 4b Real Data v3 ===")
ceu = load_samples(SAMPLES, "CEU", N_CEU)
yri = load_samples(SAMPLES, "YRI", N_YRI)
Path("/tmp/ceu_v3.txt").write_text("\n".join(ceu))
Path("/tmp/yri_v3.txt").write_text("\n".join(yri))
print(f"CEU: {len(ceu)}, YRI: {len(yri)}")

# ── Step 1: Extract Vindija variants in hg19 + liftover ──────────────────────
print("\n1. Extracting + lifting Vindija variants...", flush=True)
cmd = (
    f"bcftools view -v snps -m2 -M2 -r 21 {VCF_VIND} | "
    f"bcftools query -f '21\\t%POS0\\t%POS\\t%REF\\t%ALT\\t[%GT]\\n'"
)
r = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=180)
lines = [l for l in r.stdout.strip().split("\n")
         if l and len(l.split("\t")) >= 6 and l.split("\t")[4] not in ("", ".")]
Path("/tmp/vind_hg19_v3.bed").write_text("\n".join(lines))
print(f"  {len(lines)} Vindija biallelic SNPs extracted")

subprocess.run([CROSSMAP, "bed", CHAIN, "/tmp/vind_hg19_v3.bed", "/tmp/vind_hg38_v3.bed"],
               capture_output=True, timeout=120)

# Load lifted Vindija: key = pos_hg38, val = (ref_hg19, alt_hg19, g0, g1)
vind_map = {}
with open("/tmp/vind_hg38_v3.bed") as f:
    for line in f:
        p = line.strip().split("\t")
        if len(p) < 6: continue
        pos = int(p[2])
        if not (10_000_000 <= pos <= 20_000_000): continue
        try: g0, g1 = int(p[5][0]), int(p[5][2])
        except: continue
        if g0 > 0 or g1 > 0:
            vind_map[pos] = (p[3], p[4], g0, g1)  # ref_hg19, alt_hg19, g0, g1

print(f"  {len(vind_map)} Vindija sites with derived allele in 10-20Mb (hg38)")

# ── Step 2: Load 1000GP extraction ───────────────────────────────────────────
print("\n2. Loading 1000GP CEU+YRI allele table...", flush=True)

def parse_1kgp_table(fname, n_samples):
    """Returns pos_map: {pos: (ref, alt, hap_array)} and positions array."""
    pos_map = {}
    with open(fname) as f:
        for line in f:
            parts = line.rstrip().split("\t")
            if len(parts) < 4: continue
            pos, ref, alt = int(parts[0]), parts[1], parts[2]
            gts = []
            for gt_field in parts[3:]:
                gt_field = gt_field.strip()
                if not gt_field: continue
                a0 = gt_field[0]; a1 = gt_field[2] if len(gt_field) > 2 else a0
                try: gts.extend([int(a0), int(a1)])
                except: gts.extend([0, 0])
            if len(gts) == 2*n_samples:
                pos_map[pos] = (ref, alt, np.array(gts, dtype=np.int8))
    return pos_map

t0 = time.perf_counter()
ceu_map = parse_1kgp_table("/tmp/1kgp_chr21_10_20mb_ceu.txt", N_CEU)
yri_map = parse_1kgp_table("/tmp/1kgp_chr21_10_20mb_yri.txt", N_YRI)
print(f"  CEU: {len(ceu_map)} SNPs, YRI: {len(yri_map)} SNPs, {time.perf_counter()-t0:.1f}s")

# ── Step 3: Allele-aware intersection ────────────────────────────────────────
print("\n3. Building allele-aware intersection...", flush=True)
shared_pos = []
vind_genos = []   # (g0, g1) relative to 1000GP convention
ceu_haps = []     # (n_ceu_hap, n_sites) — will be columnar
yri_freqs = []

n_compatible = n_swapped = n_strand = n_incompatible = 0

for pos in sorted(vind_map.keys()):
    if pos not in ceu_map: continue
    ref_hg19, alt_hg19, g0, g1 = vind_map[pos]
    ref_1kgp, alt_1kgp, hap_arr = ceu_map[pos]

    # Check allele compatibility
    ref_h19 = ref_hg19.upper(); alt_h19 = alt_hg19.upper()
    ref_1k = ref_1kgp.upper(); alt_1k = alt_1kgp.upper()

    if ref_h19 == ref_1k and alt_h19 == alt_1k:
        # Direct match: Vindija's derived (ALT) = 1kGP's derived (ALT)
        vg0, vg1 = g0, g1
        n_compatible += 1
    elif ref_h19 == alt_1k and alt_h19 == ref_1k:
        # Allele swap: Vindija's ALT = 1kGP's REF → flip
        vg0, vg1 = 1-g0, 1-g1
        n_swapped += 1
    elif (COMP.get(ref_h19,"") == ref_1k and COMP.get(alt_h19,"") == alt_1k):
        # Complement match (strand flip, same orientation)
        vg0, vg1 = g0, g1
        n_strand += 1
    elif (COMP.get(ref_h19,"") == alt_1k and COMP.get(alt_h19,"") == ref_1k):
        # Complement + swap
        vg0, vg1 = 1-g0, 1-g1
        n_strand += 1
    else:
        n_incompatible += 1
        continue

    if vg0 == 0 and vg1 == 0: continue  # no derived allele in Vindija after correction
    # Only keep homozygous Vindija sites (1/1): heterozygous sites are uninformative
    if vg0 != vg1: continue

    # YRI frequency
    if pos in yri_map:
        _, _, yri_arr = yri_map[pos]
        yri_f = float(yri_arr.sum() + 0.5) / (len(yri_arr) + 1.0)
    else:
        yri_f = 0.5 / (2*N_YRI + 1.0)

    # Only keep archaic-specific sites (YRI freq < 1%) for positive-only model
    if yri_f >= 0.01: continue

    shared_pos.append(pos)
    vind_genos.append([vg0, vg1])
    ceu_haps.append(hap_arr)
    yri_freqs.append(yri_f)

print(f"  Shared (compatible): {len(shared_pos)}")
print(f"  Direct match: {n_compatible}, Swapped: {n_swapped}, Strand: {n_strand}, Incompatible: {n_incompatible}")
nspec = sum(1 for f in yri_freqs if f < 0.01)
print(f"  Vindija-specific (YRI<1%): {nspec} ({100*nspec/max(len(shared_pos),1):.1f}%)")
print(f"  Mean YRI freq: {np.mean(yri_freqs):.4f}")

if len(shared_pos) < 50:
    print("ERROR: Too few shared sites. Exiting.")
    sys.exit(1)

sarr = np.array(shared_pos, dtype=np.int32)
asub = np.array(vind_genos, dtype=np.int8)
ceu_matrix = np.stack(ceu_haps, axis=1)  # (n_hap, n_sites)
ref_freq = np.array(yri_freqs, dtype=np.float32)

print(f"  ceu_matrix: {ceu_matrix.shape}, asub: {asub.shape}")

# ── Step 4: Run ArchaicPainter ─────────────────────────────────────────────
def merge_nea(segs, gap=10_000, minlen=10_000):
    if not segs: return []
    segs = sorted(segs, key=lambda s: s.start)
    m = [Segment(segs[0].start, segs[0].end, segs[0].source, segs[0].posterior)]
    for s in segs[1:]:
        if s.start - m[-1].end <= gap:
            m[-1] = Segment(m[-1].start, max(m[-1].end, s.end), m[-1].source,
                             max(m[-1].posterior, s.posterior))
        else: m.append(Segment(s.start, s.end, s.source, s.posterior))
    return [s for s in m if s.length >= minlen]

print("\n4. Running ArchaicPainter on CEU chr21:10-20Mb...", flush=True)
gd = np.maximum(np.diff(sarr.astype(float)) * 1e-8, 1e-12)
hmm = ArchaicHMM(pi_amh=0.98, pi_nea=0.02, pi_den=0.0)

segments = []
afs = []
t_total = 0.0

for hi in range(ceu_matrix.shape[0]):
    q = ceu_matrix[hi]
    t0 = time.perf_counter()
    log_emit = np.zeros((3, len(sarr)), dtype=np.float64)
    log_emit[0] = compute_amh_emission(q, ref_freq)
    log_emit[1] = compute_archaic_emission(q, asub, positive_only=True)
    path = hmm.viterbi(log_emit, gd)
    posts, _ = hmm.forward_backward(log_emit, gd)
    t_total += time.perf_counter() - t0

    raw = decode_segments(path, sarr, posts, STATE_NAMES, min_length_bp=5_000)
    nea_raw = [Segment(s["start"], s["end"], "NEA", s["mean_posterior"])
               for s in raw if s["state"] == "NEA"]
    nea = merge_nea(nea_raw)
    af = sum(s.length for s in nea) / 10_000_000  # 10 Mb region
    afs.append(af)
    segments.extend([(ceu[hi//2], hi%2, s) for s in nea])
    if hi % 8 == 0:
        print(f"  hi={hi}: segs={len(nea)} AF={af:.4f}", flush=True)

# ── Summary ─────────────────────────────────────────────────────────────────
n_hap = ceu_matrix.shape[0]
mean_af = np.mean(afs)
lens = [s.length for _,_,s in segments]
print(f"\n{'='*60}")
print(f"  Haplotypes: {n_hap}")
print(f"  NEA segments: {len(segments)} ({len(segments)/n_hap:.1f}/hap)")
print(f"  Mean AF: {mean_af:.4f} ({100*mean_af:.2f}%)")
if lens:
    print(f"  Segment lengths: mean={np.mean(lens)/1e3:.1f}kb median={np.median(lens)/1e3:.1f}kb max={np.max(lens)/1e3:.1f}kb")
print(f"  Runtime: {t_total:.1f}s ({t_total/n_hap:.2f}s/hap)")

# Save
bed_path = OUT_DIR / "archaicpainter_chr21_10_20mb_CEU_v3.bed"
with open(bed_path, "w") as f:
    f.write("#chrom\tstart\tend\tsample\thap\tstate\tposterior\tlength_kb\n")
    for samp, hi, seg in sorted(segments, key=lambda x:(x[0],x[1],x[2].start)):
        f.write(f"chr21\t{seg.start}\t{seg.end}\t{samp}\t{hi}\tNEA\t{seg.posterior:.4f}\t{seg.length/1e3:.2f}\n")

summary = {
    "region": REGION_HG38, "n_individuals": N_CEU, "n_haplotypes": n_hap,
    "shared_sites": len(shared_pos), "nea_specific_sites": nspec,
    "allele_direct": n_compatible, "allele_swapped": n_swapped,
    "allele_strand": n_strand, "allele_incompatible": n_incompatible,
    "n_segments": len(segments),
    "mean_ancestry_fraction": float(mean_af),
    "mean_segment_length_kb": float(np.mean(lens)/1e3) if lens else 0,
    "per_haplotype_ancestry": [float(x) for x in afs],
}
with open(OUT_DIR / "realdata_summary_v3.json", "w") as f:
    json.dump(summary, f, indent=2)
print(f"\nSaved → {bed_path.name}, realdata_summary_v3.json")
