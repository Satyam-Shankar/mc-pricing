"""
Monte Carlo pricing engine with variance reduction.

Methods (all use the same number N of simulated paths, so comparisons are at
a fixed simulation budget):

plain               mean of N i.i.d. discounted payoffs
antithetic          N/2 pairs (Z, -Z); the unit of averaging is the pair mean
control             regression control variate  Y - b*(X - E[X]),
                    b* = Cov(Y, X) / Var(X) estimated from the same sample
antithetic+control  control variate applied to antithetic pair means
qmc                 randomised quasi-Monte Carlo: R independent scrambles of
                    a Sobol' sequence, N/R points each; SE from the spread of
                    the R estimates (optionally with Brownian-bridge paths)
qmc+control         as qmc, with the control variate inside each scramble

Simulation runs in batches so that memory stays bounded for any N; only
running sums of Y, X, Y^2, X^2, XY are kept, which is all the estimators need.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np
from scipy.stats import norm

from .simulate import Market, NormalSource, brownian_bridge, gbm_paths

METHODS = ["plain", "antithetic", "control", "antithetic+control", "qmc", "qmc+control"]


@dataclass
class MCResult:
    price: float
    se: float
    n_paths: int
    time_s: float
    method: str
    extra: dict = field(default_factory=dict)

    def ci(self, level: float = 0.95):
        z = norm.ppf(0.5 + level / 2)
        return self.price - z * self.se, self.price + z * self.se

    def __repr__(self):
        lo, hi = self.ci()
        return (f"MCResult({self.method}: {self.price:.5f} ± {self.se:.5f} "
                f"[{lo:.4f}, {hi:.4f}], N={self.n_paths}, {1000 * self.time_s:.1f} ms)")


class _Moments:
    """Running sums for mean / variance / covariance of (Y, X)."""

    def __init__(self):
        self.n = 0
        self.sy = self.sx = self.syy = self.sxx = self.sxy = 0.0

    def add(self, Y, X=None):
        self.n += len(Y)
        self.sy += Y.sum()
        self.syy += (Y * Y).sum()
        if X is not None:
            self.sx += X.sum()
            self.sxx += (X * X).sum()
            self.sxy += (X * Y).sum()

    def plain(self):
        my = self.sy / self.n
        vy = (self.syy - self.n * my * my) / (self.n - 1)
        return my, np.sqrt(max(vy, 0.0) / self.n)

    def control(self, EX):
        n = self.n
        my, mx = self.sy / n, self.sx / n
        vy = (self.syy - n * my * my) / (n - 1)
        vx = (self.sxx - n * mx * mx) / (n - 1)
        cxy = (self.sxy - n * mx * my) / (n - 1)
        b = cxy / vx if vx > 0 else 0.0
        est = my - b * (mx - EX)
        v_res = max(vy - 2 * b * cxy + b * b * vx, 0.0)
        rho = cxy / np.sqrt(vx * vy) if vx > 0 and vy > 0 else 0.0
        return est, np.sqrt(v_res / (n - 2)), b, rho


def price(option, mkt: Market, n_paths: int, method: str = "plain", m: int = 1,
          seed: int = 0, scheme: str = "exact", bridge: bool = False,
          n_scrambles: int = 16, batch: int = 1 << 16) -> MCResult:
    """Price `option` with n_paths simulated paths and m monitoring dates."""
    if method not in METHODS:
        raise ValueError(method)
    if not option.path_dependent and scheme == "exact":
        m = 1  # terminal value is all that is needed; simulate it exactly in one step
    t0 = time.perf_counter()
    use_cv = "control" in method
    disc = mkt.disc

    def unit_values(Z):
        if bridge and m > 1:
            Z = brownian_bridge(Z)
        S = gbm_paths(mkt, Z, scheme)
        Y = disc * option.payoff(S)
        if use_cv:
            Xr, EX = option.control(S, mkt, m)
            return Y, disc * Xr, EX
        return Y, None, None

    if method.startswith("qmc"):
        R = n_scrambles
        per = max(2, n_paths // R)
        ests = []
        info = {}
        for r in range(R):
            src = NormalSource("sobol", m, seed=seed * 1000 + r)
            mom, EX = _Moments(), None
            done = 0
            while done < per:
                k = min(batch, per - done)
                Y, X, EX = unit_values(src.draw(k))
                mom.add(Y, X)
                done += k
            if use_cv:
                est, _, b, rho = mom.control(EX)
                info = {"b": b, "rho": rho}
            else:
                est, _ = mom.plain()
            ests.append(est)
        ests = np.array(ests)
        return MCResult(float(ests.mean()), float(ests.std(ddof=1) / np.sqrt(R)), per * R,
                        time.perf_counter() - t0, method, info)

    src = NormalSource("pseudo", m, seed)
    anti = method.startswith("antithetic")
    units = n_paths // 2 if anti else n_paths
    mom, EX = _Moments(), None
    done = 0
    while done < units:
        k = min(batch, units - done)
        Z = src.draw(k)
        Y, X, EX = unit_values(Z)
        if anti:
            Y2, X2, _ = unit_values(-Z)
            Y = 0.5 * (Y + Y2)
            X = 0.5 * (X + X2) if X is not None else None
        mom.add(Y, X)
        done += k
    if use_cv:
        est, se, b, rho = mom.control(EX)
        extra = {"b": b, "rho": rho}
    else:
        (est, se), extra = mom.plain(), {}
    return MCResult(float(est), float(se), units * (2 if anti else 1),
                    time.perf_counter() - t0, method, extra)
