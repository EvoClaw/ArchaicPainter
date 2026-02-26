# Literature Review — Deep Reading Table
# Field: Archaic Introgression / Paleogenomics
# Generated: 2026-02-26

## Paper 1: Liang et al. 2025 (Genome Biology) — T2T Neanderthal Analysis
**Access**: full_text | **Venue**: Genome Biology 26:32

| Field | Content |
|-------|---------|
| **Problem** | Prior introgression analyses relied on GRCh37/GRCh38; T2T-CHM13 (complete human reference) had not been applied |
| **Method** | Remapped Altai Neanderthal + Denisovan reads to GRCh37/GRCh38/T2T-CHM13; applied IBDmix (reference-free); 2504 individuals from 26 populations (1000GP) |
| **Key results** | T2T identifies ~51Mb of Neanderthal sequences unique vs. GRCh38; 80% of these overlap with SVs distinguishing T2T from GRCh38; pre-phasing filtering strategy has 15–40% effect on introgression calls; 10 novel population-specific adaptive haplotypes (AFR, EAS, EUR-specific); created ASH database |
| **Limitations** | Only uses Altai Neanderthal (one reference); only IBDmix method; no Denisovan segment analysis beyond remapping; no multi-source disentanglement; MUSK gene example but limited functional follow-up; no AMR populations highlighted; no uncertainty quantification |
| **Open invitations** | "Incorporating a pangenome graph approach into archaic genome research in the future may provide a promising solution"; calls for HMM or S* validation of IBDmix T2T findings; long-read sequencing needed for centromeric/complex regions |

---

## Paper 2: Peede et al. 2025 (Genome Research) — Methods Review
**Access**: full_text | **Venue**: Genome Research 36:239-256

| Field | Content |
|-------|---------|
| **Problem** | Comprehensive review of state-of-art introgression methods; identifies outstanding challenges |
| **Method categories** | Site pattern / f-statistics (D, f4, fd); ARG-based (ARGweaver, ARGweaver-D, SINGER); HMM segment-based (HMMix, admixfrog, IBDmix); ML-based (IntroUNET, genomatnn, ArchIE, MaLAdapt) |
| **Key gaps explicitly stated** | (1) **Donor disambiguation**: "the field lacks a standardized, principled approach to annotate most likely archaic donor"; (2) **Multi-source at haplotype level**: "Future methods: haplotype-level, multiple archaic sources, scalable"; (3) **Admixed populations**: need to jointly infer local ancestry + introgressed tracts; (4) **ML demographic misspecification**: domain adaptation as solution; (5) **Low-coverage aDNA**: bias correction for imputed segment lengths |
| **ML methods** | IntroUNET (semantic segmentation CNN), genomatnn (CNN, adaptive introgression), ArchIE (random forest), MaLAdapt (random forest) — all trained on simulations with pre-specified demographics → sensitive to demographic misspecification |
| **ARG methods** | SINGER (2025) is breakthrough: Bayesian ARG sampling 100x faster than ARGweaver, applied to 1000GP British + African individuals; but SINGER applied narrowly; not yet comprehensive |
| **Critical observation** | ARGweaver-D handles multiple archaic sources at haplotype level but is computationally infeasible at 1000GP scale; IBDmix scales but only one source and individual-level |

---

## Paper 3: DAIseg — Planche et al. 2025 (bioRxiv)
**Access**: full_text (HTML) | **Venue**: bioRxiv March 2025

| Field | Content |
|-------|---------|
| **Problem** | Most methods can't jointly distinguish multiple archaic sources natively; require complex post-processing |
| **Method** | HMM extension of HMMix; models multiple archaic ancestries simultaneously; parameterizes emission probabilities as Poisson with mean = mutation_rate × coalescent_time; transition probabilities extended for multi-event admixture |
| **Key results** | 24% more Neanderthal + 17% more Denisovan segments detected vs. HMMix in Papuans; cleaner separation of Denisovan from Neanderthal; evidence for 2 Denisovan introgression events without post-processing; works on unphased data |
| **Limitations** | HMM architecture = window-based, no full haplotype context; limited to populations in paper (Papuans); not applied to 1000GP-scale; parameter estimation may be sensitive to demographic model choice; short segments may be missed |
| **Open invitations** | Extension to all 1000GP populations including AMR; integration with pangenome; incorporating more archaic genomes |

---

## Paper 4: Ongaro & Huerta-Sanchez 2024 (Nature Genetics)
**Access**: abstract_only | **Confidence**: high (multiple secondary sources)

| Field | Content |
|-------|---------|
| **Problem** | How many separate Denisovan populations contributed to modern humans? |
| **Key results** | Evidence for ≥3 separate Denisovan introgression events; different Denisovan populations had distinct environments and relatedness to sequenced Altai Denisovan; complex geographic distribution |
| **Implication** | The "Denisovan" in modern genomes is not a single source — creating a method challenge for any reference-based approach using only the Altai Denisovan genome |

---

## Paper 5: SINGER — Li et al. 2025 (Nature Genetics)
**Access**: abstract_only | **Confidence**: high

| Field | Content |
|-------|---------|
| **Problem** | Previous ARG inference either slow (ARGweaver) or inaccurate under model misspecification (tsinfer+tsdate, Relate) |
| **Method** | Bayesian MCMC ARG sampling; 100x faster than ARGweaver; samples posterior distribution |
| **Key results** | Applied to 1000GP British + African populations; detected archaic introgression signals + uncertainty quantification; strong support for ancient polymorphism at HLA locus |
| **Limitations** | Applied to small subset of 1000GP; not comprehensive; computationally still heavy for 2504 samples |

---

## Summary: Method Landscape (as of Feb 2026)

| Method | Multi-source | Haplotype-level | Scale (n) | Ref-free | DL-based | Year |
|--------|------------|-----------------|-----------|----------|----------|------|
| IBDmix | No (single) | No (individual) | 10,000+ | Yes | No | 2020 |
| HMMix | No (post-hoc) | Yes | 1,000+ | Yes | No | 2020 |
| DAIseg | Yes (joint) | Yes | 100s | Yes | No | 2025 |
| ARGweaver-D | Yes (joint) | Yes | ~10 | No | No | 2020 |
| SINGER | Yes (signals) | Yes (ARG) | Hundreds | No | No | 2025 |
| ArchaicSeeker 2.0 | Yes (2 sources) | No | 1,000+ | No | No | 2021 |
| IntroUNET | No | Limited | Moderate | No | Yes (CNN) | 2024 |
| **MISSING** | **Yes (3+ sources)** | **Yes** | **2,500+** | **Yes** | **Yes** | — |
