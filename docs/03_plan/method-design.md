# ArchaicPainter Method Design (Phase 3 Draft)

## Innovation Point

Core novelty: Extend Li and Stephens (2003) haplotype matching to archaic genomes
(Neanderthal, Denisovan) as reference panel members, with coalescent-theoretic
emission probabilities. This signal is absent in ALL existing introgression methods.

Existing tools:
- HMMix (Skov 2020): Poisson variant density emission; ignores LD
- DAIseg (Zhou 2025): same Poisson emission; multi-source but no haplotype structure
- IBDmix (Chen 2020 Cell): reference-free IBD; no archaic ref; no source attribution
- SparsePainter (2024): Li and Stephens for modern LAI, NEVER archaic

Theoretical core:
Under archaic introgression (t_admix ~ 1724 gen, 50 kya):
  E[match length to archaic] = 1/(r * t_admix) ~ 58 kb
Under AMH background (t_split ~ 19000 gen, 550 kya):
  E[match length to archaic] = 1/(r * t_split) ~ 5 kb
11x signal difference. Unused by all existing tools.

## Value Proposition

1. Theoretically grounded: every parameter = biological quantity
2. Higher sensitivity for short segments (<50 kb): haplotype LD vs independent SNPs
3. Multi-source HMM: states {AMH, NEA, DEN}; posterior probabilities natural output
4. Expert knowledge directly encoded in model structure
5. Scalable: PBWT for 5008 modern haplotypes; direct comparison for 6 archaic genomes

## HMM Structure

States: AMH, NEA, DEN (ghost archaic deferred to v2)
Observations: biallelic SNPs along chromosome

Emission P(o_i | state s):
  AMH: allele frequency in reference modern panel (Li and Stephens painting)
  NEA/DEN: coalescent haplotype match probability (see below)

Transition P(s_i -> s_i+1):
  P(change) = 1 - exp(-r_i_i+1 * d_s)
  r = genetic distance from Oxford map (GRCh38)
  d_s = 1/(r_mean * t_admix_s) = expected tract length

Prior: pi(AMH) ~ 0.97-0.99; pi(NEA), pi(DEN) from literature per population.

## Emission for Unphased Archaic Genomes

Archaic genomes are diploid and unphased.
At archaic position i with alleles {a1, a2}:
  Homozygous (a1==a2): P(o_i | archaic) = (1-theta) if o_i==a1, else theta
  Heterozygous (a1!=a2): P(o_i | archaic) = 0.5
    (uninformative; prevents heterozygote inflation without masking)

AMH state: standard Li and Stephens; frequency in reference panel.

## Scalability

Modern panel (5008 haplotypes): PBWT-based O(n*L) sparse matching
Archaic panel (6 genomes): direct comparison (trivially cheap)
HMM: forward-backward for posteriors; Viterbi for segments; O(L * 9) per haplotype

## Output

BED file per individual: chr, start, end, state, posterior_mean, length_bp
Posterior track: P(AMH), P(NEA), P(DEN) at each SNP position
Summary: total archaic bp, segment count per state

## Baselines (M-Step 3)

IBDmix: Chen et al. 2020 Cell; github.com/PrincetonUniversity/ibdmix
HMMix: Skov et al. 2020 Nat Genet; github.com/LauritsSkov/Introgression-detection
DAIseg: Zhou et al. 2025 bioRxiv; TBD
Sprime: Browning et al. 2018 PLOS Genet; browning/sprime.html
Simple density window: in-house 50kb sliding window SNP count (trivial floor)

## Ablation Plan (M-Step 4)

A1 (HIGH): Replace coalescent emission with Poisson density -> DEGRADE F1 on short segments (tests core claim)
A2 (HIGH): Single archaic state vs NEA+DEN -> DEGRADE source attribution
A3 (MED): EM vs fixed parameters from literature -> individual variation
A4 (MED): Oxford recombination map vs uniform -> DEGRADE boundary precision
A5 (LOW): Naive het assignment vs marginalization -> increase FP at archaic het sites
A6 (LOW): 500 vs 5008 reference haplotypes -> weaker AMH emission

## Theoretical Analysis Plan (M-Step 4b)

T1: Proposition: L&S sensitivity > Poisson for segment length L < L_crossover
    where L_crossover = f(t_admix, t_split, r, Ne)
T2: Derivation: E[match_length] = 1/(r * t_admix_gen) from coalescent theory
T3: Complexity: PBWT O(n*L) vs naive O(n^2 * L)
