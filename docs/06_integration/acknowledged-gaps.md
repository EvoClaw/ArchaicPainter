# Acknowledged Gaps — ArchaicPainter

These are experiments or analyses that were identified as important but explicitly skipped or deferred by author decision. All must appear in the paper's Limitations and/or Future Work sections.

## User-Overridden (Explicitly Skipped)

| Gap | Overridden In | Paper Treatment |
|-----|---------------|-----------------|
| Direct comparison to HMMix on same simulation | Phase 4b (tool installation failed; user said skip) | Limitations: "Direct comparison to HMMix and DAIseg was not performed; the Poisson density baseline serves as a simplified proxy for density-based methods" |
| Direct comparison to DAIseg on same simulation | Same as above | Same as above |
| Quantitative F1/AUPRC comparison to IBDmix (both on simulations with ground truth) | Phase 4b | Limitations: "IBDmix comparison is qualitative (segment count, ancestry fraction, overlap) as no common simulation ground truth was used" |

## Deferred to Future Work

| Gap | Reason for Deferral | Paper Treatment |
|-----|---------------------|-----------------|
| Full 1000GP genome-wide atlas (2,504 individuals, 26 populations) | Scope too large for current phase | Future Work: "Full genome-wide application to the 1000 Genomes Project is planned using PBWT-based acceleration" |
| Real Denisova genome integration for DEN state | Denisova VCF not in current data pipeline | Limitations: "DEN attribution results use a simulated Denisova reference; real Denisova genome integration is planned for v2" |
| Multi-chromosome generalization (beyond chr21) | Time constraints | Limitations: "Analysis is limited to chr21; genome-wide behavior is extrapolated from scalability benchmarks" |
| AMR (Admixed American) population analysis | Background admixture confounds require additional methodology | Future Work: "Application to admixed populations (AMR) requires handling of complex background admixture and is deferred" |

## Soft Limitations (Inherent to the Field)

| Gap | Nature |
|-----|--------|
| Real-data ground truth | No reliable gold standard for archaic introgression in real human genomes; validation must be indirect |
| Complex demography in simulation | Single admixture pulse used; multiple pulses, selection, background admixture not modeled |
| Reference-source divergence | Vindija 33.19 ≠ actual introgressing Neanderthal; handled by positive-only emission but creates information loss |

## Tracking

- Phase 5 deliberation: 2 rounds completed, CONDITIONAL PASS
- All REQUIRED supplements overridden by user
- Venue target adjusted to PLOS Genetics (primary) given tool comparison gap
