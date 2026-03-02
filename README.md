# ArchaicPainter

<div align="center">

**A Haplotype Matching Framework for Per-Haplotype Archaic Introgression Detection**

[![Paper PDF](https://img.shields.io/badge/📄_Paper-PDF-red?style=for-the-badge)](https://github.com/EvoClaw/ArchaicPainter/raw/main/paper/main.pdf)
[![Built with Amplify](https://img.shields.io/badge/🔬_Built_with-Amplify-blueviolet?style=for-the-badge)](https://evoclaw.github.io/amplify/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow?style=for-the-badge)](LICENSE)

> Detecting archaic hominin ancestry in modern human genomes using coalescent-theoretic haplotype match-length signals

</div>

---

## 📄 Paper

> **ArchaicPainter: A Haplotype Matching Framework for Per-Haplotype Archaic Introgression Detection**
>
> *[Authors]*

**[→ Download the full paper (PDF)](https://github.com/EvoClaw/ArchaicPainter/raw/main/paper/main.pdf)**

---

## What Is This?

Modern humans carry 1–4% Neanderthal DNA — a legacy of interbreeding that occurred ~50,000 years ago. Identifying *exactly which genomic segments* are of Neanderthal or Denisovan origin is critical for understanding the functional and adaptive consequences of archaic introgression.

**ArchaicPainter** is a probabilistic method that detects and delineates archaic introgressed segments from phased haplotype data. Unlike existing methods that count archaic variant *density*, ArchaicPainter exploits **haplotype match length** — a signal that coalescent theory predicts to be ~11× stronger in introgressed segments than in non-introgressed ones.

### Core idea

When a modern human haplotype is derived from an archaic ancestor, it will exhibit long stretches of exact matching (or near-exact matching) to the archaic reference genome, because the coalescence time is short (~50 kya). In contrast, non-introgressed segments coalesce with the archaic reference much deeper (~550 kya), producing only short, fragmented matches. ArchaicPainter captures this difference explicitly.

```
Introgressed segment:  ██████████████████████  ~58 kb match  ← easy to detect
Non-introgressed:      █  █   █  █             ~5 kb match   ← 11× shorter
```

### Method

ArchaicPainter implements a **three-state Hidden Markov Model** with states for:
- **AMH** — anatomically modern human ancestry
- **NEA** — Neanderthal ancestry  
- **DEN** — Denisovan ancestry

It extends the Li & Stephens (2003) haplotype copying framework by deriving **emission probabilities from coalescent theory** (match-length distributions under each ancestry state) and **transition probabilities** as `exp(−d·t)` where `d` is the genetic distance and `t` is admixture time. This yields a principled, parameter-efficient model with interpretable components.

Two emission modes are provided:
- **Bidirectional** (default for simulation/validation): exploits both matches *and* mismatches
- **Positive-only** (recommended for real data): appropriate when the available archaic reference diverges from the true introgressing population

---

## Key Results

| Metric | ArchaicPainter |
|--------|---------------|
| Segment F1 (50 replicates) | **0.480 ± 0.095** |
| AUPRC | **0.577** |

**Real data (Chromosome 21, 1000 Genomes Phase 3):**
- CEU (Europeans): **1.10% ± 0.60%** Neanderthal ancestry per haplotype (chr21)
- CHB (East Asians): **0.98% ± 0.52%** Neanderthal ancestry per haplotype (chr21)
- Consistent with published genome-wide estimates of 1–2%

**Scalability:**
- Near-linear empirical scaling (exponent **1.21**)
- **0.11 s/haplotype** at n = 1,000 on a single CPU core

**Functional annotation (chr21):**
- 75 named genes overlap detected introgressed regions
- Complete **interferon receptor gene cluster** (*IFNAR1*, *IFNAR2*, *IFNGR2*, *IL10RB*) recovered in Europeans — consistent with adaptive retention of archaic antiviral immunity
- 7 of 16 published adaptive introgression candidates recovered

---

## Repository Structure

```
ArchaicPainter/
├── paper/                    # Full manuscript (LaTeX source + compiled PDF)
│   ├── main.pdf              ← compiled paper
│   ├── main.tex
│   ├── preamble.tex
│   ├── sections/             # Modular section files
│   │   ├── abstract.tex
│   │   ├── introduction.tex
│   │   ├── related-work.tex
│   │   ├── method.tex
│   │   ├── results.tex
│   │   ├── discussion.tex
│   │   ├── conclusion.tex
│   │   └── figures.tex
│   ├── figures/              # All paper figures (PDF + PNG)
│   │   ├── fig1_hmm_schematic.pdf
│   │   ├── fig2_benchmark.pdf
│   │   ├── fig3_realdata.pdf
│   │   ├── fig4_scalability.pdf
│   │   └── fig5_functional_annotation.pdf
│   └── references.bib
├── docs/                     # Research documentation (Amplify phases)
│   ├── 01_intake/            # Research anchor, domain definition
│   ├── 02_literature/        # Literature review, gap analysis
│   ├── 03_plan/              # Method design, evaluation protocol
│   ├── 05_execution/         # Experiment logs, raw results (JSON)
│   └── 06_integration/       # Argument blueprint, integration notes
└── README.md
```

---

## Figures

| Figure | Description |
|--------|-------------|
| **Fig 1** | HMM architecture and 11× coalescent match-length contrast |
| **Fig 2** | Simulation benchmark: F1 box plots (n=50) and AUPRC |
| **Fig 3** | Real-data ancestry fractions (CEU/CHB) and IBDmix comparison |
| **Fig 4** | Near-linear computational scalability (log-log) |
| **Fig 5** | Functional annotation of chr21 introgressed regions |

---

## Data

- **Reference genomes**: Altai Neanderthal (Prüfer et al. 2014) and Altai Denisovan (Meyer et al. 2012), VCF lifted to GRCh38
- **Query panel**: 1000 Genomes Project Phase 3 — 40 CEU + 40 CHB haplotypes, chromosome 21 (10–48 Mb)
- **Informative sites**: filtered for archaic/modern bi-allelic sites with archaic minor allele frequency ≥ 5%

All data accession numbers and preprocessing steps are documented in `docs/03_plan/method-design.md` and Methods section of the paper.

---

## Comparison with Existing Methods

| Method | Unit | Signal used | Per-haplotype | Neanderthal + Denisovan |
|--------|------|-------------|:---:|:---:|
| HMMix | Diploid | Archaic variant density | ✗ | ✗ |
| DAIseg | Diploid | Archaic allele HMM | ✗ | ✗ |
| IBDmix | Diploid | IBD-based scoring | ✗ | ✗ |
| S'PRIME | Haplotype | LD / haplotype | ✓ | ✗ |
| ArchIE | Haplotype | Logistic regression | ✓ | ✗ |
| **ArchaicPainter** | **Haplotype** | **Li & Stephens match length** | **✓** | **✓** |

---

## Citation

If you use ArchaicPainter in your work, please cite:

```bibtex
@article{archaicpainter2026,
  title   = {{ArchaicPainter}: A Haplotype Matching Framework for
             Per-Haplotype Archaic Introgression Detection},
  author  = {[Authors]},
  year    = {2026},
  note    = {Preprint}
}
```

---

## Built with Amplify

This research — from literature review to experiment design, execution, results integration, and paper writing — was conducted autonomously using **[Amplify](https://evoclaw.github.io/amplify/)**, an open-source agentic research automation framework built on top of Cursor IDE.

> Amplify turns your AI coding assistant into an autonomous research agent that conducts the full scientific workflow with built-in rigor and human oversight: literature review → problem validation → method design → experiment execution → results integration → paper writing.

[![Built with Amplify](https://img.shields.io/badge/🔬_Automated_Research_with-Amplify-blueviolet?style=flat-square)](https://evoclaw.github.io/amplify/)

The entire research pipeline (7 phases, 4 mandatory gates, 34 verified citations, multi-agent deliberation) was executed by Claude claude-4.6-sonnet-medium-thinking via Amplify Skillsets. **No human editing of the paper text.** For demonstration purposes — please treat results with appropriate scientific caution.

---

## License

MIT License — see [LICENSE](LICENSE) for details.
