"""
Random-number generation and path simulation for geometric Brownian motion.

Normals
-------
pseudo  : NumPy PCG64 pseudo-random standard normals
sobol   : scrambled Sobol' points mapped through the inverse normal CDF
          (randomised quasi-Monte Carlo)

Path construction
-----------------
exact   : log-Euler scheme, exact in distribution for GBM at every date
euler   : Euler-Maruyama on S itself (has O(dt) weak discretisation bias)
bridge  : Brownian-bridge ordering of the normals, so the first (best
          distributed) Sobol' coordinates fix the large-scale shape of the
          path - this is what makes QMC effective for path-dependent payoffs.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import norm, qmc


@dataclass(frozen=True)
class Market:
    S0: float = 100.0
    r: float = 0.05
    sigma: float = 0.2
    T: float = 1.0
    q: float = 0.0

    @property
    def disc(self) -> float:
        return float(np.exp(-self.r * self.T))


class NormalSource:
    """Streams standard-normal blocks of shape (n, dim)."""

    def __init__(self, kind: str, dim: int, seed: int):
        self.kind, self.dim = kind, dim
        if kind == "pseudo":
            self.rng = np.random.default_rng(seed)
        elif kind == "sobol":
            self.engine = qmc.Sobol(d=dim, scramble=True, seed=seed)
        else:
            raise ValueError(kind)

    def draw(self, n: int) -> np.ndarray:
        if self.kind == "pseudo":
            return self.rng.standard_normal((n, self.dim))
        u = self.engine.random(n)
        return norm.ppf(np.clip(u, 1e-16, 1 - 1e-16))


def brownian_bridge(Z: np.ndarray) -> np.ndarray:
    """Re-map i.i.d. normals so column 0 fixes W(T), column 1 W(T/2), ...
    Returns increments dW/sqrt(dt) (again i.i.d. N(0,1) in distribution)."""
    n, m = Z.shape
    t = np.arange(m + 1, dtype=float)          # time in units of dt, t[0] = 0
    W = np.zeros((n, m + 1))
    W[:, m] = np.sqrt(m) * Z[:, 0]
    k = 1
    queue = [(0, m)]
    while queue:
        nxt = []
        for left, right in queue:
            if right - left < 2:
                continue
            mid = (left + right) // 2
            a, b = t[mid] - t[left], t[right] - t[mid]
            W[:, mid] = (b * W[:, left] + a * W[:, right]) / (a + b) + np.sqrt(a * b / (a + b)) * Z[:, k]
            k += 1
            nxt += [(left, mid), (mid, right)]
        queue = nxt
    return np.diff(W, axis=1)


def gbm_paths(mkt: Market, Z: np.ndarray, scheme: str = "exact") -> np.ndarray:
    """Prices S(t_1..t_m) for t_i = i T/m driven by standard normals Z (n, m)."""
    n, m = Z.shape
    dt = mkt.T / m
    if scheme == "exact":
        inc = (mkt.r - mkt.q - 0.5 * mkt.sigma ** 2) * dt + mkt.sigma * np.sqrt(dt) * Z
        return mkt.S0 * np.exp(np.cumsum(inc, axis=1))
    if scheme == "euler":
        return mkt.S0 * np.cumprod(1.0 + (mkt.r - mkt.q) * dt + mkt.sigma * np.sqrt(dt) * Z, axis=1)
    raise ValueError(scheme)
