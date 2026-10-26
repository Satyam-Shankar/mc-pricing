"""
Option contracts. Each contract knows its payoff on simulated monitoring
dates and, for the control-variate estimator, a correlated payoff with a
known expectation.

    payoff(S)          -> undiscounted payoff per path, S has shape (n, m)
    control(S, mkt, m) -> (undiscounted control payoff per path, its exact
                           *discounted* expectation)
    exact(mkt, m)      -> analytical price if one exists, else None
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .analytic import (black_scholes, down_and_out_call, down_and_out_call_discrete,
                       geometric_asian)
from .simulate import Market


@dataclass(frozen=True)
class EuropeanCall:
    K: float
    name: str = "European call"
    path_dependent: bool = False

    def payoff(self, S):
        return np.maximum(S[:, -1] - self.K, 0.0)

    def control(self, S, mkt: Market, m: int):
        # discounted stock price is a martingale: E[e^{-rT} S_T] = S0 e^{-qT}
        return S[:, -1], mkt.S0 * np.exp(-mkt.q * mkt.T)

    def exact(self, mkt: Market, m: int = 1):
        return float(black_scholes(mkt.S0, self.K, mkt.r, mkt.sigma, mkt.T, mkt.q, call=True))


@dataclass(frozen=True)
class EuropeanPut:
    K: float
    name: str = "European put"
    path_dependent: bool = False

    def payoff(self, S):
        return np.maximum(self.K - S[:, -1], 0.0)

    def control(self, S, mkt: Market, m: int):
        return S[:, -1], mkt.S0 * np.exp(-mkt.q * mkt.T)

    def exact(self, mkt: Market, m: int = 1):
        return float(black_scholes(mkt.S0, self.K, mkt.r, mkt.sigma, mkt.T, mkt.q, call=False))


@dataclass(frozen=True)
class ArithmeticAsianCall:
    """Payoff max(mean_i S(t_i) - K, 0). No closed form; the geometric-average
    Asian call (closed form) is an almost perfectly correlated control."""
    K: float
    name: str = "Arithmetic Asian call"
    path_dependent: bool = True

    def payoff(self, S):
        return np.maximum(S.mean(axis=1) - self.K, 0.0)

    def control(self, S, mkt: Market, m: int):
        G = np.exp(np.log(S).mean(axis=1))
        return np.maximum(G - self.K, 0.0), float(
            geometric_asian(mkt.S0, self.K, mkt.r, mkt.sigma, mkt.T, m, mkt.q))

    def exact(self, mkt: Market, m: int = 1):
        return None


@dataclass(frozen=True)
class GeometricAsianCall:
    K: float
    name: str = "Geometric Asian call"
    path_dependent: bool = True

    def payoff(self, S):
        return np.maximum(np.exp(np.log(S).mean(axis=1)) - self.K, 0.0)

    def control(self, S, mkt: Market, m: int):
        return S[:, -1], mkt.S0 * np.exp(-mkt.q * mkt.T)

    def exact(self, mkt: Market, m: int = 1):
        return float(geometric_asian(mkt.S0, self.K, mkt.r, mkt.sigma, mkt.T, m, mkt.q))


@dataclass(frozen=True)
class DownAndOutCall:
    """Knocked out if S(t_i) <= H at any monitoring date. Control: the
    vanilla European call with the same strike."""
    K: float
    H: float
    name: str = "Down-and-out call"
    path_dependent: bool = True

    def payoff(self, S):
        alive = S.min(axis=1) > self.H
        return np.where(alive, np.maximum(S[:, -1] - self.K, 0.0), 0.0)

    def control(self, S, mkt: Market, m: int):
        return np.maximum(S[:, -1] - self.K, 0.0), float(
            black_scholes(mkt.S0, self.K, mkt.r, mkt.sigma, mkt.T, mkt.q))

    def exact(self, mkt: Market, m: int = 1):
        """Discrete-monitoring price via the BGK continuity correction
        (an accurate approximation, not an identity)."""
        return down_and_out_call_discrete(mkt.S0, self.K, self.H, mkt.r, mkt.sigma, mkt.T, m, mkt.q)

    def continuous(self, mkt: Market):
        return down_and_out_call(mkt.S0, self.K, self.H, mkt.r, mkt.sigma, mkt.T, mkt.q)
