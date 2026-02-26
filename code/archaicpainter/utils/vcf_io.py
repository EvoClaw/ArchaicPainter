"""
VCF I/O utilities: parse 1000GP and archaic VCF files into numpy arrays.
Handles GRCh37/38, biallelic SNP filtering, and phased/unphased genotype extraction.
"""
import numpy as np
import pysam
import logging
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional, List, Tuple, Dict

logger = logging.getLogger(__name__)


@dataclass
class VariantData:
    """Container for a single chromosome's variant data."""
    chrom: str
    positions: np.ndarray           # (L,) physical positions (bp)
    genetic_pos: np.ndarray         # (L,) genetic positions (cM)
    ref_alleles: List[str]          # (L,) reference alleles
    alt_alleles: List[str]          # (L,) alternate alleles
    haplotypes: np.ndarray          # (N_hap, L) int8 {0, 1, -1=missing}
    sample_ids: List[str]           # (N_samples,)
    n_samples: int = field(init=False)
    n_hap: int = field(init=False)
    n_sites: int = field(init=False)

    def __post_init__(self):
        self.n_samples = len(self.sample_ids)
        self.n_hap = self.haplotypes.shape[0]
        self.n_sites = self.haplotypes.shape[1]


@dataclass
class ArchaicVariantData:
    """Container for archaic genome variants (unphased diploid)."""
    chrom: str
    positions: np.ndarray           # (L,) physical positions
    ref_alleles: List[str]
    alt_alleles: List[str]
    genotypes: np.ndarray           # (L, 2) int8 alleles {0, 1, -1=missing}
    sample_id: str

    @property
    def n_sites(self) -> int:
        return len(self.positions)

    def allele_pair(self, i: int) -> Tuple[int, int]:
        """Return (a1, a2) at site i; -1 for missing."""
        return int(self.genotypes[i, 0]), int(self.genotypes[i, 1])


def read_1000gp_vcf(
    vcf_path: str,
    chrom: str,
    samples: Optional[List[str]] = None,
    min_af: float = 0.0,
    max_af: float = 1.0,
) -> VariantData:
    """
    Read phased 1000GP VCF into a VariantData object.

    Parameters
    ----------
    vcf_path : str
        Path to bgzipped+tabix-indexed VCF.
    chrom : str
        Chromosome name (e.g. '21' or 'chr21').
    samples : list, optional
        Subset of sample IDs. If None, load all.
    min_af, max_af : float
        Allele frequency filters (applied to ALT allele).
    """
    vcf_path = Path(vcf_path)
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF not found: {vcf_path}")

    vcf = pysam.VariantFile(str(vcf_path))

    # Subset samples
    if samples is not None:
        available = set(vcf.header.samples)
        missing = set(samples) - available
        if missing:
            logger.warning(f"Samples not found in VCF: {missing}")
        samples = [s for s in samples if s in available]
        vcf.subset_samples(samples)
    else:
        samples = list(vcf.header.samples)

    n_samples = len(samples)

    positions = []
    ref_list = []
    alt_list = []
    haplotype_rows = []   # list of (2*N,) arrays

    for rec in vcf.fetch(chrom):
        # Biallelic SNPs only
        if len(rec.alts) != 1:
            continue
        ref, alt = rec.ref, rec.alts[0]
        if len(ref) != 1 or len(alt) != 1:
            continue  # skip indels

        # Compute AF for filter
        ac = rec.info.get("AC", None)
        an = rec.info.get("AN", None)
        if ac is not None and an is not None and an > 0:
            af = ac[0] / an
            if af < min_af or af > max_af:
                continue

        row = np.full(2 * n_samples, -1, dtype=np.int8)
        for s_idx, sample in enumerate(samples):
            gt = rec.samples[sample]["GT"]
            if gt[0] is not None:
                row[2 * s_idx] = int(gt[0])
            if gt[1] is not None:
                row[2 * s_idx + 1] = int(gt[1])

        positions.append(rec.pos)
        ref_list.append(ref)
        alt_list.append(alt)
        haplotype_rows.append(row)

    vcf.close()

    if len(positions) == 0:
        raise ValueError(f"No biallelic SNPs found on {chrom} in {vcf_path}")

    positions_arr = np.array(positions, dtype=np.int32)
    hap_matrix = np.stack(haplotype_rows, axis=1)  # (2*N, L)

    logger.info(
        f"Loaded {n_samples} samples ({2*n_samples} haplotypes), "
        f"{len(positions)} SNPs from {chrom}"
    )

    return VariantData(
        chrom=chrom,
        positions=positions_arr,
        genetic_pos=np.zeros(len(positions_arr), dtype=np.float32),  # filled later
        ref_alleles=ref_list,
        alt_alleles=alt_list,
        haplotypes=hap_matrix,
        sample_ids=samples,
    )


def read_archaic_vcf(
    vcf_path: str,
    chrom: str,
    sample_id: Optional[str] = None,
) -> ArchaicVariantData:
    """
    Read an archaic genome VCF (Neanderthal/Denisovan) into ArchaicVariantData.

    These genomes are diploid but unphased. Heterozygous sites contribute 0.5
    emission probability (uninformative) per the method design.
    """
    vcf_path = Path(vcf_path)
    if not vcf_path.exists():
        raise FileNotFoundError(f"Archaic VCF not found: {vcf_path}")

    vcf = pysam.VariantFile(str(vcf_path))
    if sample_id is None:
        sample_id = list(vcf.header.samples)[0]

    positions = []
    ref_list = []
    alt_list = []
    geno_rows = []

    for rec in vcf.fetch(chrom):
        if len(rec.alts) != 1:
            continue
        ref, alt = rec.ref, rec.alts[0]
        if len(ref) != 1 or len(alt) != 1:
            continue

        gt = rec.samples[sample_id]["GT"]
        a1 = int(gt[0]) if gt[0] is not None else -1
        a2 = int(gt[1]) if gt[1] is not None else -1

        positions.append(rec.pos)
        ref_list.append(ref)
        alt_list.append(alt)
        geno_rows.append([a1, a2])

    vcf.close()

    logger.info(
        f"Loaded archaic {sample_id}: {len(positions)} SNPs on {chrom}"
    )

    return ArchaicVariantData(
        chrom=chrom,
        positions=np.array(positions, dtype=np.int32),
        ref_alleles=ref_list,
        alt_alleles=alt_list,
        genotypes=np.array(geno_rows, dtype=np.int8),
        sample_id=sample_id,
    )


def intersect_sites(
    modern: VariantData,
    archaic: ArchaicVariantData,
) -> Tuple[VariantData, ArchaicVariantData]:
    """
    Restrict to biallelic SNP positions shared between modern and archaic VCFs,
    requiring REF/ALT allele consistency.
    """
    archaic_pos_to_idx = {p: i for i, p in enumerate(archaic.positions)}
    keep_modern = []
    keep_archaic = []

    for i, pos in enumerate(modern.positions):
        if pos not in archaic_pos_to_idx:
            continue
        j = archaic_pos_to_idx[pos]
        if modern.ref_alleles[i] != archaic.ref_alleles[j]:
            continue
        if modern.alt_alleles[i] != archaic.alt_alleles[j]:
            continue
        keep_modern.append(i)
        keep_archaic.append(j)

    if len(keep_modern) == 0:
        raise ValueError("No overlapping SNPs between modern and archaic VCF")

    km = np.array(keep_modern)
    ka = np.array(keep_archaic)

    modern_sub = VariantData(
        chrom=modern.chrom,
        positions=modern.positions[km],
        genetic_pos=modern.genetic_pos[km],
        ref_alleles=[modern.ref_alleles[i] for i in km],
        alt_alleles=[modern.alt_alleles[i] for i in km],
        haplotypes=modern.haplotypes[:, km],
        sample_ids=modern.sample_ids,
    )
    archaic_sub = ArchaicVariantData(
        chrom=archaic.chrom,
        positions=archaic.positions[ka],
        ref_alleles=[archaic.ref_alleles[i] for i in ka],
        alt_alleles=[archaic.alt_alleles[i] for i in ka],
        genotypes=archaic.genotypes[ka, :],
        sample_id=archaic.sample_id,
    )

    logger.info(f"Intersection: {len(km)} shared SNPs")
    return modern_sub, archaic_sub
