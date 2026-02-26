"""
Evaluation metrics for ArchaicPainter (locked in evaluation-protocol.yaml).

Primary metrics:
  C1: Segment-level F1 @ IoU threshold 0.5
  C2: AUPRC (area under precision-recall curve)
  C3: FDR at 5% fixed recall threshold
  C4: Multi-source attribution accuracy (NEA vs DEN)

Secondary metrics:
  - Runtime (wall seconds)
  - Memory (peak RSS MB)
  - Archaic ancestry fraction

Segment matching uses 50% reciprocal overlap (IoU >= 0.5) as the true-positive
criterion, matching the evaluation-protocol.yaml specification.
"""
import numpy as np
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass
import time
import tracemalloc


@dataclass
class Segment:
    """A genomic segment with source label."""
    start: int
    end: int
    source: str  # "AMH", "NEA", "DEN", or "archaic" (source-agnostic)
    posterior: float = 1.0

    @property
    def length(self) -> int:
        return max(0, self.end - self.start)

    def overlaps(self, other: "Segment", iou_thresh: float = 0.5) -> bool:
        """True if reciprocal overlap (IoU) >= iou_thresh."""
        overlap_start = max(self.start, other.start)
        overlap_end = min(self.end, other.end)
        overlap_len = max(0, overlap_end - overlap_start)
        if overlap_len == 0:
            return False
        union_len = self.length + other.length - overlap_len
        if union_len <= 0:
            return False
        return (overlap_len / union_len) >= iou_thresh


def parse_bed_segments(
    bed_lines: List[str],
    source_col: Optional[int] = 4,
) -> List[Segment]:
    """
    Parse BED lines into Segment objects.

    Expects columns: chrom, start, end, [sample, source, ...]
    """
    segments = []
    for line in bed_lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        start = int(parts[1])
        end = int(parts[2])
        source = parts[source_col] if source_col is not None and len(parts) > source_col else "archaic"
        posterior = float(parts[6]) if len(parts) > 6 else 1.0
        segments.append(Segment(start=start, end=end, source=source, posterior=posterior))
    return segments


def segment_f1(
    true_segments: List[Segment],
    pred_segments: List[Segment],
    iou_thresh: float = 0.5,
    source_filter: Optional[str] = None,
) -> Dict[str, float]:
    """
    Compute segment-level Precision, Recall, and F1.

    Matching: each predicted segment is matched to at most one true segment
    using greedy IoU matching.

    Parameters
    ----------
    source_filter : str, optional
        If given, only evaluate segments with this source label (e.g., "NEA").
    """
    if source_filter:
        true_segments = [s for s in true_segments if s.source == source_filter]
        pred_segments = [s for s in pred_segments if s.source == source_filter]

    if not true_segments and not pred_segments:
        return {"precision": 1.0, "recall": 1.0, "f1": 1.0, "n_true": 0, "n_pred": 0}
    if not true_segments:
        return {"precision": 0.0, "recall": 1.0, "f1": 0.0, "n_true": 0, "n_pred": len(pred_segments)}
    if not pred_segments:
        return {"precision": 1.0, "recall": 0.0, "f1": 0.0, "n_true": len(true_segments), "n_pred": 0}

    matched_true = set()
    tp = 0
    for pred in pred_segments:
        for j, true in enumerate(true_segments):
            if j not in matched_true and pred.overlaps(true, iou_thresh):
                tp += 1
                matched_true.add(j)
                break

    precision = tp / len(pred_segments) if pred_segments else 0.0
    recall = tp / len(true_segments) if true_segments else 0.0
    f1 = (
        2 * precision * recall / (precision + recall)
        if (precision + recall) > 0
        else 0.0
    )
    return {
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "n_true": len(true_segments),
        "n_pred": len(pred_segments),
        "tp": tp,
    }


def compute_auprc(
    true_segments: List[Segment],
    pred_segments: List[Segment],
    iou_thresh: float = 0.5,
) -> float:
    """
    Area under the Precision-Recall curve.

    Sweeps posterior threshold from 0 to 1 to build the PR curve.
    """
    if not true_segments:
        return 0.0
    if not pred_segments:
        return 0.0

    # Sort predictions by posterior descending
    sorted_preds = sorted(pred_segments, key=lambda s: s.posterior, reverse=True)
    thresholds = sorted(set(s.posterior for s in sorted_preds), reverse=True)

    precisions = [1.0]
    recalls = [0.0]

    for thresh in thresholds:
        active_preds = [s for s in sorted_preds if s.posterior >= thresh]
        metrics = segment_f1(true_segments, active_preds, iou_thresh)
        precisions.append(metrics["precision"])
        recalls.append(metrics["recall"])

    precisions.append(0.0)
    recalls.append(1.0)

    # Trapezoidal integration (sklearn-style)
    auprc = float(np.trapz(precisions, recalls))
    return abs(auprc)


def fdr_at_fixed_recall(
    true_segments: List[Segment],
    pred_segments: List[Segment],
    target_recall: float = 0.05,
    iou_thresh: float = 0.5,
) -> float:
    """
    False Discovery Rate at the posterior threshold achieving >= target_recall.

    FDR = FP / (TP + FP) = 1 - Precision
    Returns FDR, or 1.0 if target recall cannot be reached.
    """
    if not true_segments or not pred_segments:
        return 1.0

    sorted_preds = sorted(pred_segments, key=lambda s: s.posterior, reverse=True)
    thresholds = sorted(set(s.posterior for s in sorted_preds), reverse=True)

    for thresh in thresholds:
        active_preds = [s for s in sorted_preds if s.posterior >= thresh]
        metrics = segment_f1(true_segments, active_preds, iou_thresh)
        if metrics["recall"] >= target_recall:
            fdr = 1.0 - metrics["precision"]
            return fdr

    return 1.0


def attribution_accuracy(
    true_segments: List[Segment],
    pred_segments: List[Segment],
    iou_thresh: float = 0.5,
) -> float:
    """
    Source attribution accuracy: among matched segments, fraction with correct
    source label (NEA vs DEN).

    Requires that true_segments and pred_segments both have non-trivial source labels.
    """
    correct = 0
    total = 0
    for pred in pred_segments:
        if pred.source not in ("NEA", "DEN"):
            continue
        for true in true_segments:
            if true.source not in ("NEA", "DEN"):
                continue
            if pred.overlaps(true, iou_thresh):
                total += 1
                if pred.source == true.source:
                    correct += 1
                break

    return correct / total if total > 0 else 0.0


def compute_ancestry_fraction(
    segments: List[Segment],
    chrom_length: int,
    source: Optional[str] = None,
) -> float:
    """
    Fraction of chromosome covered by archaic (or specific-source) segments.
    """
    if source is not None:
        segments = [s for s in segments if s.source == source]
    total_bp = sum(s.length for s in segments)
    return total_bp / chrom_length if chrom_length > 0 else 0.0


@dataclass
class EvaluationResult:
    """Complete evaluation result for one method on one simulation replicate."""
    method: str
    replicate: int
    seed: int
    f1_all: float
    f1_short: float         # segments < 50 kb — primary claim C1
    precision_all: float
    recall_all: float
    auprc: float            # C2
    fdr_5pct: float         # C3
    attribution_acc: float  # C4
    ancestry_fraction: float
    runtime_s: float
    memory_mb: float
    n_true: int
    n_pred: int


def evaluate_full(
    true_segments: List[Segment],
    pred_segments: List[Segment],
    chrom_length: int,
    method: str,
    replicate: int = 0,
    seed: int = 0,
    runtime_s: float = 0.0,
    memory_mb: float = 0.0,
) -> EvaluationResult:
    """
    Compute all locked metrics for one method/replicate combination.
    """
    all_metrics = segment_f1(true_segments, pred_segments)
    short_true = [s for s in true_segments if s.length < 50_000]
    short_pred = [s for s in pred_segments if s.length < 50_000]
    short_metrics = segment_f1(short_true, short_pred)

    auprc = compute_auprc(true_segments, pred_segments)
    fdr = fdr_at_fixed_recall(true_segments, pred_segments, target_recall=0.05)
    attr_acc = attribution_accuracy(true_segments, pred_segments)
    anc_frac = compute_ancestry_fraction(pred_segments, chrom_length)

    return EvaluationResult(
        method=method,
        replicate=replicate,
        seed=seed,
        f1_all=all_metrics["f1"],
        f1_short=short_metrics["f1"],
        precision_all=all_metrics["precision"],
        recall_all=all_metrics["recall"],
        auprc=auprc,
        fdr_5pct=fdr,
        attribution_acc=attr_acc,
        ancestry_fraction=anc_frac,
        runtime_s=runtime_s,
        memory_mb=memory_mb,
        n_true=all_metrics["n_true"],
        n_pred=all_metrics["n_pred"],
    )
