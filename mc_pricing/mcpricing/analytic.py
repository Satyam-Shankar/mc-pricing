"""
Closed-form prices used as ground truth.

All formulas assume geometric Brownian motion under the risk-neutral measure,

    dS_t = (r - q) S_t dt + sigma S_t dW_t,

with constant interest rate r, dividend yield q and volatility sigma.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import norm

BGK_BETA = 0.5825971579390106  # -zeta(1/2)/sqrt(2*pi), Broadie-Glasserman-Kou


def black_scholes(S0, K, r, sigma, T, q=0.0, call=True):
    """European call/put price (Black-Scholes-Merton)."""
    S0, K, sigma, T = map(np.asarray, (S0, K, sigma, T))
    sq = sigma * np.sqrt(T)
    d1 = (np.log(S0 / K) + (r - q + 0.5 * sigma ** 2) * T) / sq
    d2 = d1 - sq
    if call:
        return S0 * np.exp(-q * T) * norm.cdf(d1) - K * np.exp(-r * T) * norm.cdf(d2)
    return K * np.exp(-r * T) * norm.cdf(-d2) - S0 * np.exp(-q * T) * norm.cdf(-d1)


def bs_delta(S0, K, r, sigma, T, q=0.0, call=True):
    d1 = (np.log(S0 / K) + (r - q + 0.5 * sigma ** 2) * T) / (sigma * np.sqrt(T))
    return np.exp(-q * T) * (norm.cdf(d1) if call else norm.cdf(d1) - 1)


def geometric_asian(S0, K, r, sigma, T, m, q=0.0, call=True):
    """Discretely monitored geometric-average Asian option, averaging over
    S(t_i), t_i = i T / m, i = 1..m. log G is normal, so the price is exact."""
    mu = np.log(S0) + (r - q - 0.5 * sigma ** 2) * T * (m + 1) / (2 * m)
    var = sigma ** 2 * T * (m + 1) * (2 * m + 1) / (6 * m ** 2)
    sd = np.sqrt(var)
    d1 = (mu - np.log(K) + var) / sd
    d2 = d1 - sd
    disc = np.exp(-r * T)
    if call:
        return disc * (np.exp(mu + 0.5 * var) * norm.cdf(d1) - K * norm.cdf(d2))
    return disc * (K * norm.cdf(-d2) - np.exp(mu + 0.5 * var) * norm.cdf(-d1))


def down_and_out_call(S0, K, H, r, sigma, T, q=0.0):
    """Continuously monitored down-and-out call, barrier H <= K, no rebate
    (Merton 1973; Reiner & Rubinstein 1991)."""
    if H >= S0:
        return 0.0
    if H > K:
        raise NotImplementedError("formula implemented for H <= K")
    lam = (r - q + 0.5 * sigma ** 2) / sigma ** 2
    sq = sigma * np.sqrt(T)
    y = np.log(H ** 2 / (S0 * K)) / sq + lam * sq
    c_di = (S0 * np.exp(-q * T) * (H / S0) ** (2 * lam) * norm.cdf(y)
            - K * np.exp(-r * T) * (H / S0) ** (2 * lam - 2) * norm.cdf(y - sq))
    return float(black_scholes(S0, K, r, sigma, T, q) - c_di)


def down_and_out_call_discrete(S0, K, H, r, sigma, T, m, q=0.0):
    """Approximation for *discrete* monitoring at m dates: the continuous
    formula with the barrier shifted to H * exp(-beta sigma sqrt(T/m))
    (Broadie, Glasserman & Kou 1997)."""
    return down_and_out_call(S0, K, H * np.exp(-BGK_BETA * sigma * np.sqrt(T / m)), r, sigma, T, q)
