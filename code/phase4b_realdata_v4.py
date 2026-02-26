"""
Phase 4b Real Data v4 - Full chr21 (10-48Mb) x CEU + CHB
Extends v3 to full euchromatic chr21 and adds CHB population.
"""
import sys, time, json, warnings, subprocess
import numpy as np
from pathlib import Path
warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))

from archaicpainter.core.emission import compute_amh_emission, compute_archaic_emission
from archaicpainter.core.hmm import ArchaicHMM, decode_segments, STATE_NAMES
from archaicpainter.evaluation.metrics import Segment

VCF_1KGP = "/home/yanlin/livestock/1000GP/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
VCF_VIND = "/home/yanlin/livestock/data/archaic/vindija33.19_chr21.vcf.gz"
CHAIN    = "/home/yanlin/livestock/data/archaic/hg19ToHg38.over.chain.gz"
CROSSMAP = "/home/yanlin/miniconda3/bin/CrossMap"
SAMPLES  = "/home/yanlin/livestock/1000GP/samples.info"
OUT_DIR  = Path("/home/yanlin/livestock/docs/05_execution")

REGION_HG38 = "chr21:10000000-48000000"
REGION_LEN  = 38_000_000
N_CEU = 20
N_CHB = 20
N_YRI = 50
COMP = {"A":"T","T":"A","C":"G","G":"C"}


def load_samples(fname, pop, n=None):
    s = [l.split()[0] for l in open(fname) if len(l.split()) >= 2 and l.split()[1] == pop]
    return s[:n] if n else s


print("=== Phase 4b Real Data v4: full chr21 CEU+CHB ===")
ceu = load_samples(SAMPLES, "CEU", N_CEU)
chb = load_samples(SAMPLES, "CHB", N_CHB)
yri = load_samples(SAMPLES, "YRI", N_YRI)
Path("/tmp/ceu_v4.txt").write_text("\n".join(ceu))
Path("/tmp/chb_v4.txt").write_text("\n".join(chb))
Path("/tmp/yri_v4.txt").write_text("\n".join(yri))
print("  CEU:{} CHB:{} YRI:{}".format(len(ceu), len(chb), len(yri)))

print("\n1. Extracting Vindija variants...", flush=True)
cmd = ("bcftools view -v snps -m2 -M2 -r 21 " + VCF_VIND +
       " | bcftools query -f '21\\t%POS0\\t%POS\\t%REF\\t%ALT\\t[%GT]\\n'")
r = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=300)
lines = [l for l in r.stdout.strip().split("\n")
         if l and len(l.split("\t")) >= 6 and l.split("\t")[4] not in ("", ".")]
Path("/tmp/vind_hg19_v4.bed").write_text("\n".join(lines))
print("  {} Vindija SNPs extracted".format(len(lines)))
subprocess.run([CROSSMAP, "bed", CHAIN, "/tmp/vind_hg19_v4.bed", "/tmp/vind_hg38_v4.bed"],
               capture_output=True, timeout=180)

vind_map = {}
with open("/tmp/vind_hg38_v4.bed") as f:
    for line in f:
        p = line.strip().split("\t")
        if len(p) < 6: continue
        pos = int(p[2])
        if not (10_000_000 <= pos <= 48_000_000): continue
        try: g0, g1 = int(p[5][0]), int(p[5][2])
        except: continue
        if g0 > 0 or g1 > 0:
            vind_map[pos] = (p[3], p[4], g0, g1)
print("  {} Vindija derived sites 10-48Mb".format(len(vind_map)))


print("\n2. Querying 1000GP chr21:10-48Mb...", flush=True)


def query_1kgp(pop_file, label, n_samples):
    cmd = ("bcftools view -v snps -m2 -M2 -r chr21:10000000-48000000 -S " + pop_file +
           " " + VCF_1KGP +
           " | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%GT\\t]\\n'")
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=900)
    pos_map = {}
    for line in r.stdout.strip().split("\n"):
        if not line: continue
        parts = line.split("\t")
        if len(parts) < 5: continue
        try:
            pos = int(parts[1]); ref = parts[2]; alt = parts[3]; gts = parts[4:]
            haps = []
            for g in gts:
                g = g.strip().rstrip("\t")
                sep = "|" if "|" in g else "/"
                if sep in g:
                    a, b = g.split(sep)
                    try: haps.extend([int(a), int(b)])
                    except: haps.extend([0, 0])
                else:
                    haps.extend([0, 0])
            arr = np.array(haps[:n_samples * 2], dtype=np.int8)
            if arr.sum() > 0:
                pos_map[pos] = (ref, alt, arr)
        except: continue
    print("  {}: {} sites".format(label, len(pos_map)))
    return pos_map


ceu_map = query_1kgp("/tmp/ceu_v4.txt", "CEU", N_CEU)
chb_map = query_1kgp("/tmp/chb_v4.txt", "CHB", N_CHB)
yri_map = query_1kgp("/tmp/yri_v4.txt", "YRI", N_YRI)


print("\n3. Allele-aware intersections...", flush=True)


def build_intersection(pop_map, n_pop, yri_m, vind_m, label):
    sp = []; vg = []; ph = []; yf = []
    nc = nsw = nst = ni = 0
    for pos in sorted(vind_m.keys()):
        if pos not in pop_map: continue
        rh19, ah19, g0, g1 = vind_m[pos]
        r1k, a1k, hap = pop_map[pos]
        rh = rh19.upper(); ah = ah19.upper(); r1 = r1k.upper(); a1 = a1k.upper()
        if rh == r1 and ah == a1:
            vg0, vg1 = g0, g1; nc += 1
        elif rh == a1 and ah == r1:
            vg0, vg1 = 1 - g0, 1 - g1; nsw += 1
        elif COMP.get(rh, "") == r1 and COMP.get(ah, "") == a1:
            vg0, vg1 = g0, g1; nst += 1
        elif COMP.get(rh, "") == a1 and COMP.get(ah, "") == r1:
            vg0, vg1 = 1 - g0, 1 - g1; nst += 1
        else:
            ni += 1; continue
        if vg0 == 0 and vg1 == 0: continue
        if vg0 != vg1: continue
        if pos in yri_m:
            _, _, ya = yri_m[pos]; yfrac = float(ya.sum() + 0.5) / (len(ya) + 1.0)
        else:
            yfrac = 0.5 / (2 * N_YRI + 1.0)
        if yfrac >= 0.01: continue
        sp.append(pos); vg.append([vg0, vg1])
        ph.append(hap[:n_pop * 2]); yf.append(yfrac)
    print("  {} shared={} compat={} swap={} strand={} incompat={}".format(
        label, len(sp), nc, nsw, nst, ni))
    return sp, vg, ph, yf


ceu_sh, ceu_vg, ceu_hp, ceu_yf = build_intersection(ceu_map, N_CEU, yri_map, vind_map, "CEU")
chb_sh, chb_vg, chb_hp, chb_yf = build_intersection(chb_map, N_CHB, yri_map, vind_map, "CHB")


def merge_nea(segs, gap=10_000, minlen=10_000):
    if not segs: return []
    segs = sorted(segs, key=lambda s: s.start)
    m = [Segment(segs[0].start, segs[0].end, segs[0].source, segs[0].posterior)]
    for s in segs[1:]:
        if s.start - m[-1].end <= gap:
            m[-1] = Segment(m[-1].start, max(m[-1].end, s.end),
                             m[-1].source, max(m[-1].posterior, s.posterior))
        else:
            m.append(Segment(s.start, s.end, s.source, s.posterior))
    return [s for s in m if s.length >= minlen]


def run_pop(sh, vgn, hps, yfreqs, pop_name, snames, n_pop):
    if len(sh) < 50:
        print("  {}: too few sites ({})".format(pop_name, len(sh)))
        return [], []
    sarr = np.array(sh, dtype=np.int32)
    asub = np.array(vgn, dtype=np.int8)
    mat = np.stack(hps, axis=1)
    ref_freq = np.array(yfreqs, dtype=np.float32)
    gd = np.maximum(np.diff(sarr.astype(float)) * 1e-8, 1e-12)
    hmm = ArchaicHMM(pi_amh=0.98, pi_nea=0.02, pi_den=0.0)
    print("\n{}: {} haplotypes, {} sites...".format(pop_name, mat.shape[0], len(sarr)), flush=True)
    segments = []; afs = []; t_total = 0.0
    for hi in range(mat.shape[0]):
        q = mat[hi]; t0 = time.perf_counter()
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
        af = sum(s.length for s in nea) / REGION_LEN
        afs.append(af)
        sn = snames[hi // 2] if hi // 2 < len(snames) else "s{}".format(hi // 2)
        segments.extend([(sn, hi % 2, s) for s in nea])
        if hi % 10 == 0:
            print("  {} hi={}: segs={} AF={:.4f}".format(pop_name, hi, len(nea), af), flush=True)
    lens = [s.length for _, _, s in segments]
    print("  {} done: {} segs AF={:.4f} ({:.2f}%) rt={:.1f}s".format(
        pop_name, len(segments), np.mean(afs) if afs else 0,
        100 * (np.mean(afs) if afs else 0), t_total))
    if lens:
        print("  Lengths: mean={:.1f}kb median={:.1f}kb max={:.1f}kb".format(
            np.mean(lens) / 1e3, np.median(lens) / 1e3, np.max(lens) / 1e3))
    return segments, afs


print("\n4. Running ArchaicPainter...", flush=True)
ceu_segs, ceu_afs = run_pop(ceu_sh, ceu_vg, ceu_hp, ceu_yf, "CEU", ceu, N_CEU)
chb_segs, chb_afs = run_pop(chb_sh, chb_vg, chb_hp, chb_yf, "CHB", chb, N_CHB)

print("\n" + "=" * 60)
print("POPULATION COMPARISON (chr21 10-48Mb):")
for pop, segs, afs in [("CEU", ceu_segs, ceu_afs), ("CHB", chb_segs, chb_afs)]:
    if not segs:
        print("  {}: no segments".format(pop)); continue
    lens = [s.length for _, _, s in segs]
    print("  {:4}: segs={:>5}  AF={:.3f}%  med={:.1f}kb  mean={:.1f}kb".format(
        pop, len(segs), 100 * np.mean(afs), np.median(lens) / 1e3, np.mean(lens) / 1e3))


def save_bed(segs, path):
    with open(str(path), "w") as f:
        f.write("#chrom\tstart\tend\tsample\thap\tstate\tposterior\tlength_kb\n")
        for sn, hi, seg in sorted(segs, key=lambda x: (x[0], x[1], x[2].start)):
            f.write("chr21\t{}\t{}\t{}\t{}\tNEA\t{:.4f}\t{:.2f}\n".format(
                seg.start, seg.end, sn, hi, seg.posterior, seg.length / 1e3))


save_bed(ceu_segs, OUT_DIR / "archaicpainter_chr21_full_CEU_v4.bed")
save_bed(chb_segs, OUT_DIR / "archaicpainter_chr21_full_CHB_v4.bed")


def pop_sum(segs, afs, sh, n_pop, pop):
    lens = [s.length for _, _, s in segs] if segs else []
    return {
        "population": pop, "region": REGION_HG38,
        "n_individuals": n_pop, "n_haplotypes": n_pop * 2,
        "shared_vindija_sites": len(sh),
        "n_segments": len(segs),
        "mean_ancestry_fraction": float(np.mean(afs)) if afs else 0.0,
        "std_ancestry_fraction": float(np.std(afs)) if afs else 0.0,
        "mean_segment_length_kb": float(np.mean(lens) / 1e3) if lens else 0,
        "median_segment_length_kb": float(np.median(lens) / 1e3) if lens else 0,
        "per_haplotype_ancestry": [float(x) for x in afs],
    }


summary = {
    "CEU": pop_sum(ceu_segs, ceu_afs, ceu_sh, N_CEU, "CEU"),
    "CHB": pop_sum(chb_segs, chb_afs, chb_sh, N_CHB, "CHB"),
}
with open(str(OUT_DIR / "realdata_summary_v4.json"), "w") as f:
    json.dump(summary, f, indent=2)
print("\nSaved -> realdata_summary_v4.json")
