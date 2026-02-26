# ArchaicPainter Technical Evaluation

**Evaluator perspective:** Senior ML/computational biology engineer; HMM-based genomic tools; methods that work in practice.  
**Target venue:** Genome Biology. **Research type:** Type H.

---

## 1. IMPLEMENTABILITY

### Design clarity: **Mostly clean, with one critical ambiguity**

**Unambiguous:**
- HMM structure (3 states, biallelic SNPs) is standard
- Transition model P(change) = 1 - exp(-r·d_s) is well-defined; Oxford map + d_s from literature
- Archaic emission for homozygous sites: P = (1-θ) if match, θ if mismatch — clear
- Heterozygous archaic: P = 0.5 marginalization — correct and avoids phase assumption
- Forward-backward for posteriors, Viterbi for segments — textbook

**Ambiguous / underspecified:**
- **AMH emission:** Spec says "derived allele frequency in modern reference panel (PBWT-based)". This is ambiguous:
  - **Interpretation A (allele frequency):** P(obs|AMH) = f_derived at that SNP. Simple, but **loses LD** — equivalent to independent-SNP model. Contradicts "Li & Stephens" framing.
  - **Interpretation B (L&S via PBWT):** Emission = L&S copy probability: P(obs|copy h) = (1-θ) if match, θ if mismatch; marginalize over panel. PBWT finds *which* haplotypes match; emission is computed from sparse set (SparsePainter-style). **This preserves LD.**
- **Recommendation:** Lock Interpretation B. Document explicitly: "AMH emission = Li & Stephens emission with sparse PBWT approximation (Q longest matches per site)." If Interpretation A is intended, the core novelty claim (haplotype LD > Poisson) is undermined.

### Main engineering challenges

| Challenge | Severity | Notes |
|-----------|----------|-------|
| **PBWT → emission bridge** | HIGH | PBWT outputs set-maximal/long matches, not emission probabilities. Need explicit algorithm: ReportQLongestMatches → which haplotypes → L&S emission formula. SparsePainter does this; reuse or reimplement. |
| **Coordinate systems** | MED | Query (phased target), modern panel (5008), archaic (6 diploid, unphased) must share same SNP set. Archaic genomes have different coverage; missing sites need handling (mask or marginalize). |
| **Recombination map** | MED | Oxford map is per-bp; need r_i for each SNP pair. Map may not cover all SNPs; interpolation strategy required. |
| **Parameter grid** | LOW | θ, d_s, π(AMH) — literature values exist; EM on chr21 is optional. Baum-Welch for 3-state HMM is standard. |

---

## 2. PBWT SPECIFICS

### Is PBWT used correctly?

**Yes, for its intended role:** PBWT efficiently finds long haplotype matches between query and panel in O(n·L). For 5008 haplotypes × 5M SNPs, naive O(n²·L) is infeasible; PBWT is the right tool.

### How does PBWT help Li & Stephens?

- **Classical L&S:** At each site, emission = Σ_h [P(copy h) × P(obs | copy h)]. Requires considering all n haplotypes.
- **PBWT contribution:** Identifies which haplotypes have long matches covering the current site. Only those contribute meaningfully; others have negligible weight under L&S (mismatch-dominated).
- **Sparse approximation (SparsePainter):** Keep only Q longest matches per site. Forward-backward runs over Q "effective" templates instead of n. Complexity drops from O(n·L) to O(Q·L).

### "PBWT-based emission" — theory vs implementation gap

**Gap:** PBWT does **not** output emission probabilities. It outputs *match intervals* (start, end, haplotype id).

**Required bridge:**
1. At site k, run ReportQLongestMatches (or similar) → get haplotypes {h_1, …, h_Q} with match lengths.
2. For each h_j: P(obs_k | copy h_j) = (1-θ) if query matches h_j at k, else θ.
3. AMH emission: P(obs_k | AMH) = Σ_j w_j · P(obs_k | copy h_j), with weights w_j from L&S (or simplified uniform over Q).
4. SparsePainter uses a HashMap over these Q haplotypes and "vanishing mutation rate" for transitions.

**Conclusion:** The design is sound **if** the implementation explicitly implements this bridge. The phrase "PBWT-based emission" is shorthand; the actual emission is *derived from PBWT match output* via L&S formula. Document this in the method section to avoid reviewer confusion.

---

## 3. HMM STABILITY

### Numerical issues

| Issue | Risk | Mitigation |
|-------|------|------------|
| **Underflow in forward-backward** | HIGH | 5M SNPs × 3 states: product of 5M probabilities underflows in float64. **Mandatory:** Log-space forward-backward (log-sum-exp for normalization). |
| **Emission zeros** | MED | θ = 1e-4 ensures no exact zeros. Archaic het P = 0.5 is fine. |
| **Transition matrix** | LOW | 3×3; well-conditioned. |
| **θ too small** | LOW | θ = 1e-4 is standard. If θ → 0, log(θ) = -9.2; manageable. |

**Recommendation:** Implement log-space forward-backward from day one. Use `scipy.special.logsumexp` or equivalent. Checkpointing (O(log L) memory) is available if full forward/backward matrices are too large, though for 3 states × 5M × 8 bytes × 2 ≈ 240 MB per chromosome, it may be acceptable to store.

---

## 4. RESOURCE FIT

### Memory (order-of-magnitude)

| Component | Size |
|-----------|------|
| Reference panel (5008 hap × 5M SNPs, binary) | ~3 GB |
| PBWT structures | ~3–6 GB (run-length arrays) |
| Archaic genomes (6 × 5M × 2 alleles) | ~60 MB |
| Forward-backward per chr (3 × 5M × 8 × 2) | ~240 MB |
| Per-individual peak (22 chr, streaming) | ~5–10 GB |
| **Total for pipeline** | **~20–30 GB** |

**Verdict:** 1 TB RAM is more than sufficient. 367 GB disk for intermediate files (PBWT, precomputed emissions) is tight but workable.

### Runtime

- **PBWT build:** O(n·L) ≈ 5008 × 5e6 × 22 ≈ 5e11 operations. Single-threaded ~hours; parallelizable.
- **Per-haplotype forward-backward:** O(L·K²) = 5e6 × 9 ≈ 45e6 ops. 5008 haplotypes × 22 chr ≈ 110k runs.
- **Rough estimate:** 110k × 0.1 sec ≈ 3 hours per chromosome for FB; 22 chr × 3 h ≈ 66 hours single-threaded. With 128 cores: **~30–60 min per individual** for HMM. PBWT + I/O adds overhead.
- **2504 individuals:** ~20–50 CPU-days total. With 128 cores: **~4–10 hours wall-clock** for full cohort.

**Verdict:** 8–12 month timeline is dominated by method development, validation, and paper writing — not compute. Resources are adequate.

### GPU utilization

The design is CPU-bound (PBWT, forward-backward). GPUs are not naturally used. Options: (a) accept CPU-only; (b) explore GPU-accelerated PBWT (research project); (c) use GPUs for simulation/validation (e.g., msprime, SLiM). For core method: **CPU is fine.**

---

## 5. FAILURE MODE ANALYSIS

### Top 3 ways this breaks in practice

**1. Archaic heterozygosity washes out NEA/DEN signal**

- **Mechanism:** At heterozygous archaic sites, P(obs|NEA) = P(obs|DEN) = 0.5. Archaic genomes (Altai, Vindija, Denisovan) have non-trivial het rates from low coverage and damage. If 20–30% of archaic sites are het, the emission is uninformative at those sites.
- **Impact:** Reduced power to call short archaic segments; inflated false negatives in archaic-rich regions.
- **Mitigation:** (a) Use high-quality archaic calls (GQ filter); (b) consider phased archaic where available (e.g., from pedigree); (c) sensitivity analysis: mask archaic hets and compare segment calls.

**2. Reference panel composition biases AMH emission**

- **Mechanism:** 1000GP has strong population structure. A haplotype from an under-represented population may have few long matches in the panel. Sparse PBWT returns Q matches that might all be poor. AMH emission becomes artificially low → bias toward NEA/DEN.
- **Impact:** False positive archaic calls in populations poorly represented in the panel (e.g., Indigenous Americans, some African populations).
- **Mitigation:** (a) Use population-specific or diverse panels where possible; (b) validate on simulated admixed individuals with known panel composition; (c) report population-specific precision/recall in evaluation.

**3. Recombination map error distorts segment boundaries**

- **Mechanism:** Transition rate uses r from Oxford map. Map errors (especially in low-LD or poorly studied regions) cause wrong P(change). Overestimated r → too many state switches → fragmented segments; underestimated r → too few switches → bloated segments.
- **Impact:** Degraded boundary precision; affects segment length statistics and downstream analyses.
- **Mitigation:** (a) Compare Oxford vs. other maps (e.g., HapMap, 1000GP); (b) ablation: uniform map vs. Oxford (in your plan); (c) use map from same reference build as genotypes.

---

## 6. SUGGESTED IMPROVEMENTS

1. **Clarify AMH emission in writing:** State explicitly: "AMH emission = Li & Stephens P(obs|copy h) marginalized over Q longest PBWT matches; weights from L&S copy process." Avoid "allele frequency" unless you mean the simpler model (and accept LD loss).

2. **Precompute emission tables:** For each (query_haplotype, chromosome), compute emission at each SNP for AMH (from PBWT) and NEA/DEN (from archaic comparison) once. Store in compressed form. Forward-backward then reads precomputed emissions — cleaner and debuggable.

3. **Log-space from the start:** Do not implement linear-space forward-backward first. Use log-space from the initial prototype.

4. **Missing data in archaic:** Define policy: if archaic has no call at site i, use P(obs_i | archaic) = 0.5 (uninformative) or marginalize over possible alleles. Document and keep consistent.

5. **Segment post-processing:** Raw Viterbi can produce single-SNP flicker. Add a minimum segment length filter (e.g., 10 kb) or run-length encoding merge for segments shorter than expected tract length.

6. **Checkpointing for very long chromosomes:** Chr1 has ~6M SNPs. If memory becomes an issue, implement checkpointing (e.g., store every 100k SNPs) to bound memory at O(log L).

---

## 7. VERDICT

### **CONDITIONAL PASS**

**Conditions:**
1. **Lock AMH emission:** Implement full Li & Stephens emission with sparse PBWT (Interpretation B), not raw allele frequency. Document clearly.
2. **Log-space HMM:** Implement log-space forward-backward in the first version.
3. **Validate failure modes:** Include in evaluation protocol: (a) archaic het rate sensitivity; (b) panel composition (e.g., AMR with Indigenous reference); (c) map comparison.

**Rationale:**
- Design is theoretically sound and novel (L&S for archaic; coalescent tract-length signal).
- PBWT usage is appropriate; the bridge to emission must be explicit.
- Resources and timeline are adequate.
- Main risks (archaic het, panel bias, map error) are identifiable and mitigatable.

**If conditions are met:** Strong candidate for Genome Biology. Method is implementable, scalable, and addresses a clear gap (haplotype-based archaic introgression with multi-source attribution).
