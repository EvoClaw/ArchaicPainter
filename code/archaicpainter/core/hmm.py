"""
Multi-source HMM for archaic introgression detection (ArchaicPainter).

States: AMH (0) | NEA (1) | DEN (2)

Transition model (per method-design.md):
  P(leave state s) = 1 - exp(-r_{i,i+1} * d_s)
  r = genetic distance in Morgans
  d_s = expected tract length = 1 / (r_mean * t_admix_s)

  Transition matrix at gap i:
    T[s, s'] = (1 - p_stay[s]) * pi[s'] / (1 - pi[s])  if s != s'
    T[s, s]  = p_stay[s]

Forward-backward and Viterbi algorithms implemented in log-space.
"""
import numpy as np
from typing import Optional, Tuple, Dict
import logging

logger = logging.getLogger(__name__)

# State indices
STATE_AMH = 0
STATE_NEA = 1
STATE_DEN = 2
N_STATES = 3
STATE_NAMES = ["AMH", "NEA", "DEN"]


class ArchaicHMM:
    """
    Three-state HMM for introgression detection.

    Parameters
    ----------
    pi_amh : float
        Prior probability of AMH state (typically 0.97-0.99).
    pi_nea : float
        Prior probability of NEA state.
    pi_den : float
        Prior probability of DEN state.
    t_admix_nea : float
        Generations since Neanderthal admixture (default: 1724, ~50 kya).
    t_admix_den : float
        Generations since Denisovan admixture (default: 1552, ~45 kya).
    r_mean_per_gen : float
        Mean recombination rate (Morgans/bp), used for tract length computation.
    """

    def __init__(
        self,
        pi_amh: float = 0.980,
        pi_nea: float = 0.018,
        pi_den: float = 0.002,
        t_admix_nea: float = 1724.0,
        t_admix_den: float = 1552.0,
        r_mean_per_gen: float = 1e-8,
    ):
        total = pi_amh + pi_nea + pi_den
        self.pi = np.array([pi_amh, pi_nea, pi_den], dtype=np.float64) / total
        # Clamp pi to avoid log(0); zero-probability states get log_pi = -1e30
        self.log_pi = np.where(self.pi > 0, np.log(np.maximum(self.pi, 1e-300)), -1e30)

        # Expected tract lengths in MORGANS = 1/t_gen
        # P(stay in state s across gap of r Morgans) = exp(-r / d_s) = exp(-r * t_s)
        #
        # d_nea = 1/t_admix_nea = 1/1724 M  (~58 kb, Neanderthal introgression)
        # d_den = 1/t_admix_den = 1/1552 M  (~64 kb, Denisovan introgression)
        #
        # d_amh: AMH is the stable background state covering ~98% of the genome.
        # Expected AMH block size between introgression events = 1/(pi_nea * t_admix_nea).
        # Using d_amh = 1/(pi_nea * t_admix_nea) gives p_stay_amh ≈ 0.999 per 4 kb gap,
        # ensuring AMH is stable while NEA/DEN transitions occur at the right rate.
        pi_nea_eff = max(self.pi[STATE_NEA], 1e-6)
        pi_den_eff = max(self.pi[STATE_DEN], 1e-6)
        pi_arch = pi_nea_eff + pi_den_eff
        self.d_amh = 1.0 / (pi_arch * max(t_admix_nea, t_admix_den))
        self.d_nea = 1.0 / t_admix_nea
        self.d_den = 1.0 / t_admix_den

        self.t_admix_nea = t_admix_nea
        self.t_admix_den = t_admix_den
        self.r_mean_per_gen = r_mean_per_gen

    def _build_transition_matrix(
        self,
        genetic_distances: np.ndarray,   # (L-1,) in Morgans
    ) -> np.ndarray:
        """
        Build transition matrix for each gap.

        Returns shape (L-1, N_STATES, N_STATES).
        T[gap, s, s'] = P(state at site i+1 = s' | state at site i = s)
        """
        L_minus1 = len(genetic_distances)

        # P(stay) per state per gap: each state has its own expected tract length
        p_stay_amh = np.exp(-genetic_distances / self.d_amh)
        p_stay_nea = np.exp(-genetic_distances / self.d_nea)
        p_stay_den = np.exp(-genetic_distances / self.d_den)

        p_stay = np.stack([p_stay_amh, p_stay_nea, p_stay_den], axis=1)  # (L-1, 3)

        # When leaving state s: jump to s' proportional to pi[s'] / (1 - pi[s])
        T = np.zeros((L_minus1, N_STATES, N_STATES), dtype=np.float64)
        for s in range(N_STATES):
            T[:, s, s] = p_stay[:, s]
            p_leave = 1.0 - p_stay[:, s]
            for sp in range(N_STATES):
                if sp != s:
                    # Normalise to sum-to-1 over off-diagonal entries
                    pi_sp_given_leave = self.pi[sp] / (1.0 - self.pi[s] + 1e-12)
                    T[:, s, sp] = p_leave * pi_sp_given_leave

        return T

    def forward_backward(
        self,
        log_emission: np.ndarray,        # (N_STATES, L)
        genetic_distances: np.ndarray,   # (L-1,) Morgans
    ) -> Tuple[np.ndarray, np.ndarray, float]:
        """
        Vectorized log-space forward-backward.

        Uses scipy.special.logsumexp with numpy broadcasting for speed:
        the inner loop is O(N^2) per step, vectorized over states.
        """
        from scipy.special import logsumexp

        L = log_emission.shape[1]
        log_T = np.log(self._build_transition_matrix(genetic_distances) + 1e-300)
        # log_T shape: (L-1, N, N)  — log_T[t, from, to]

        # Forward pass (vectorized)
        log_alpha = np.full((N_STATES, L), -np.inf, dtype=np.float64)
        log_alpha[:, 0] = self.log_pi + log_emission[:, 0]

        for t in range(1, L):
            # log_alpha[:, t-1] shape: (N,)
            # log_T[t-1] shape: (N, N) — axis 0=from, axis 1=to
            # We want logsumexp over 'from' for each 'to':
            # result[to] = logsumexp(log_alpha[from, t-1] + log_T[t-1, from, to])
            log_alpha[:, t] = (
                logsumexp(log_alpha[:, t-1][:, None] + log_T[t-1], axis=0)
                + log_emission[:, t]
            )

        log_likelihood = float(logsumexp(log_alpha[:, L-1]))

        # Backward pass (vectorized)
        log_beta = np.full((N_STATES, L), 0.0, dtype=np.float64)

        for t in range(L-2, -1, -1):
            # log_T[t] shape: (N, N) axis 0=from, axis 1=to
            # log_beta[s, t] = logsumexp over 'to' of log_T[t, s, to] + emit[to, t+1] + log_beta[to, t+1]
            log_beta[:, t] = logsumexp(
                log_T[t] + log_emission[:, t+1][None, :] + log_beta[:, t+1][None, :],
                axis=1,
            )

        # Posterior
        log_post = log_alpha + log_beta
        log_post -= logsumexp(log_post, axis=0, keepdims=True)
        posteriors = np.exp(np.clip(log_post, -500, 0))

        return posteriors, log_likelihood

    def viterbi(
        self,
        log_emission: np.ndarray,        # (N_STATES, L)
        genetic_distances: np.ndarray,   # (L-1,) Morgans
    ) -> np.ndarray:
        """
        Vectorized Viterbi decoding.
        """
        L = log_emission.shape[1]
        log_T = np.log(self._build_transition_matrix(genetic_distances) + 1e-300)

        dp = np.full((N_STATES, L), -np.inf, dtype=np.float64)
        ptr = np.zeros((N_STATES, L), dtype=np.int8)

        dp[:, 0] = self.log_pi + log_emission[:, 0]

        for t in range(1, L):
            # scores[from, to] = dp[from, t-1] + log_T[t-1, from, to]
            scores = dp[:, t-1][:, None] + log_T[t-1]   # (N, N)
            ptr[:, t] = np.argmax(scores, axis=0)         # best prev state for each current state
            dp[:, t] = scores[ptr[:, t], np.arange(N_STATES)] + log_emission[:, t]

        # Backtrack
        path = np.zeros(L, dtype=np.int8)
        path[L-1] = np.argmax(dp[:, L-1])
        for t in range(L-2, -1, -1):
            path[t] = ptr[path[t+1], t+1]

        return path


def decode_segments(
    path: np.ndarray,        # (L,) state indices
    positions: np.ndarray,   # (L,) bp positions
    posteriors: np.ndarray,  # (N_STATES, L)
    state_names: Optional[list] = None,
    min_length_bp: int = 10_000,
) -> list:
    """
    Convert Viterbi state path to BED-format segment list.

    Parameters
    ----------
    min_length_bp : int
        Minimum segment length filter (10 kb default per Phase 4b plan adjustments).

    Returns
    -------
    segments : list of dict
        Each dict: {chrom, start, end, state, length_bp, mean_posterior}
    """
    if state_names is None:
        state_names = ["AMH", "NEA", "DEN"]

    L = len(path)
    segments = []
    if L == 0:
        return segments

    seg_start = 0
    cur_state = path[0]

    for i in range(1, L + 1):
        if i == L or path[i] != cur_state:
            seg_end = positions[i - 1]
            seg_begin = positions[seg_start]
            length_bp = seg_end - seg_begin

            if length_bp >= min_length_bp:
                mean_post = float(posteriors[cur_state, seg_start:i].mean())
                segments.append({
                    "start": int(seg_begin),
                    "end": int(seg_end),
                    "state": state_names[cur_state],
                    "length_bp": int(length_bp),
                    "mean_posterior": round(mean_post, 4),
                })

            if i < L:
                seg_start = i
                cur_state = path[i]

    return segments
