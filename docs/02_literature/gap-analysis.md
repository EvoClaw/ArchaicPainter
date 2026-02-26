# Gap Analysis — Archaic Introgression Research
# Generated: 2026-02-26

## Gap 1: No Scalable Haplotype-Level Multi-Source Joint Inference
**Status**: CONFIRMED (explicitly stated in Peede et al. 2025)
**Confidence**: high (full-text review + DAIseg paper)

- IBDmix: scales to 10,000+, but single-source and individual-level (not haplotype)
- ARGweaver-D: haplotype-level, multi-source, but computational infeasibility limits to ~10 individuals
- DAIseg: joint multi-source HMM, haplotype-level, but window-based (no long-range LD), limited scale testing
- **Gap**: A method that is SIMULTANEOUSLY: (a) haplotype-level, (b) multi-source (Altai/Vindija/Chagyrskaya Neanderthal + Denisovan + ghost), (c) uses full sequence context (not just variant density windows), (d) scales to 2,500+ individuals

**Why not solved?**:
- Computational cost: ARG inference is O(n² per site)
- Demographic misspecification: DL methods fail when trained-on demography ≠ real demography
- Multiple Neanderthal genomes only recently all available at high coverage

**Our approach**: Deep learning (transformer/CNN) on haplotype sequences with domain adaptation for demographic robustness

---

## Gap 2: Lack of Principled Donor Attribution (Multi-Neanderthal)
**Status**: CONFIRMED (explicitly stated in Peede et al. 2025)
**Confidence**: high

- We have 3 high-coverage Neanderthal genomes (Altai ~500ka, Vindija 33.19 ~50ka, Chagyrskaya ~80ka) representing different populations and times
- Current practice: compare to ONE Neanderthal reference, or use ad-hoc match rate thresholds
- No method systematically asks which specific Neanderthal population contributed to which modern human population
- "the field lacks a standardized, principled approach" (Peede et al. 2025)

**Biological question unlocked**: Did East Asians receive Neanderthal DNA from a different Neanderthal population than Europeans? Is there evidence that some modern populations have Chagyrskaya-specific (vs. Vindija-specific) ancestry?

---

## Gap 3: Archaic Introgression in Recently Admixed Populations (Americas)
**Status**: CONFIRMED (explicitly stated in Peede et al. 2025 + Villanea & Witt 2022)
**Confidence**: high

- 1000GP has 7 AMR populations: MXL, CLM, PUR, PEL, IBS, PEL, etc.
- These have mixed African + European + Indigenous American ancestry
- European ancestry carries Neanderthal; African ancestry may carry ghost archaic
- Existing methods can't properly disentangle: which Neanderthal alleles came via European ancestry vs. which might have come via a distinct introgression into Indigenous American ancestors?
- "There is a need for methods that can distinguish introgressed segments originating from different genetic backgrounds" (Peede et al. 2025)
- DAIseg demo'd Demography 2 (modern + archaic admixture) but not applied to 1000GP AMR at genome-wide scale

**Biological question**: Are there archaic alleles in AMR populations that preferentially segregate in the Indigenous American ancestry background? Do introgressed alleles show different selection patterns depending on the ancestry background in which they sit?

---

## Gap 4: Pangenome-Based Archaic Introgression Detection
**Status**: IDENTIFIED (suggested as future direction in Liang et al. 2025)
**Confidence**: medium

- T2T paper itself says: "Incorporating a pangenome graph approach... may provide a promising solution"
- HPRC (Human Pangenome Reference Consortium) now has 47 diverse haplotype assemblies
- Current approaches: align archaic reads to linear reference → miss structural variations → lose complex regions
- Pangenome graph can represent all known human variation, allowing archaic reads to align to the path most similar to archaic sequences
- **Gap**: No study has aligned archaic hominin reads to human pangenome graph to detect SV-mediated introgression

**Challenge**: Archaic reads are short-read Illumina; pangenome graph alignment of short reads in complex regions is still imperfect

---

## Gap 5: Temporal Selection Dynamics on Specific Introgressed Haplotypes
**Status**: IDENTIFIED (partially addressed by recent work)
**Confidence**: medium

- We know introgressed alleles were negatively selected overall (introgression deserts)
- We know some were positively selected (adaptive introgression)
- But we don't know: for the positive ones, WHEN did selection occur? Was it immediate? Or delayed by 10,000-20,000 years?
- SINGER provides ARG with uncertainty quantification — this enables extracting local tree topologies that can be combined with selection inference
- Bayesian time-series allele frequency inference now possible (arXiv 2502.12279)
- **Gap**: Systematic analysis of selection timing for ALL known adaptive introgressed haplotypes using ARG + WF diffusion

---

## Gap 6: Ghost Archaic Introgression in ALL 26 1000GP Populations
**Status**: PARTIALLY ADDRESSED but incomplete
**Confidence**: medium

- Ghost introgression (from unknown archaic populations) identified in Africa (Durvasula & Sankararaman 2019, Genome Biology 2019)
- But systematic analysis of all 26 1000GP populations for ghost archaic signals is lacking
- In principle, after accounting for known Neanderthal + Denisovan ancestry, residual archaic signal = ghost archaic
- Requires: reliable multi-source subtraction of known archaic → characterize residual signal

**Biological implication**: Are there additional archaic populations that contributed to specific non-African populations beyond the known Denisovans?
