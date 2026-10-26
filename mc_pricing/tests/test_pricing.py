"""Correctness tests. Run with `python tests/test_pricing.py` (or pytest)."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from mcpricing import (METHODS, ArithmeticAsianCall, DownAndOutCall, EuropeanCall,  # noqa: E402
                       EuropeanPut, GeometricAsianCall, Market, black_scholes,
                       down_and_out_call, geometric_asian, price)
from mcpricing.simulate import brownian_bridge  # noqa: E402


def test_black_scholes_reference_value_and_parity():
    # Hull, textbook value: S=K=100, r=5%, sigma=20%, T=1 -> 10.4506
    c = black_scholes(100, 100, 0.05, 0.2, 1.0)
    assert abs(c - 10.4506) < 1e-4
    for K in (80, 100, 130):
        c = black_scholes(100, K, 0.03, 0.3, 2.0, q=0.01)
        p = black_scholes(100, K, 0.03, 0.3, 2.0, q=0.01, call=False)
        assert np.isclose(c - p, 100 * np.exp(-0.02) - K * np.exp(-0.06))


def test_geometric_asian_with_one_date_is_black_scholes():
    assert np.isclose(geometric_asian(100, 95, 0.05, 0.25, 1.5, m=1),
                      black_scholes(100, 95, 0.05, 0.25, 1.5))


def test_barrier_limits():
    # a barrier far below spot is (almost) never hit -> vanilla price
    assert np.isclose(down_and_out_call(100, 100, 1e-6, 0.05, 0.2, 1.0),
                      black_scholes(100, 100, 0.05, 0.2, 1.0))
    # price decreases as the barrier rises towards spot
    vals = [down_and_out_call(100, 100, h, 0.05, 0.2, 1.0) for h in (60, 80, 90, 95, 99)]
    assert all(a > b for a, b in zip(vals, vals[1:]))


def test_all_methods_unbiased_against_closed_forms():
    mkt = Market(sigma=0.3, T=0.75)
    for opt in (EuropeanCall(105), EuropeanPut(95)):
        ex = opt.exact(mkt)
        for meth in METHODS:
            r = price(opt, mkt, 2 ** 15, meth, seed=11)
            assert abs(r.price - ex) < 4 * r.se + 1e-4, (opt.name, meth, r, ex)
    opt = GeometricAsianCall(100)
    ex = opt.exact(mkt, 12)
    for meth in METHODS:
        r = price(opt, mkt, 2 ** 14, meth, m=12, seed=5, bridge=meth.startswith("qmc"))
        assert abs(r.price - ex) < 4 * r.se + 1e-4, (meth, r, ex)


def test_discrete_barrier_matches_bgk_and_asian_cv_is_strong():
    mkt = Market()
    opt = DownAndOutCall(100, 90)
    r = price(opt, mkt, 2 ** 16, "antithetic+control", m=50, seed=3)
    assert abs(r.price - opt.exact(mkt, 50)) < 4 * r.se + 0.01
    a = ArithmeticAsianCall(100)
    p = price(a, mkt, 2 ** 13, "plain", m=12, seed=1)
    c = price(a, mkt, 2 ** 13, "control", m=12, seed=1)
    assert c.se < p.se / 10 and c.extra["rho"] > 0.99


def test_brownian_bridge_gives_iid_normals():
    Z = np.random.default_rng(0).standard_normal((200_000, 12))
    B = brownian_bridge(Z)
    C = np.cov(B.T)
    assert np.allclose(C, np.eye(12), atol=0.02)


def test_standard_error_definition():
    r = price(EuropeanCall(100), Market(), 1000, "plain", seed=4)
    Z = np.random.default_rng(4).standard_normal((1000, 1))
    mkt = Market()
    ST = mkt.S0 * np.exp((mkt.r - 0.5 * mkt.sigma ** 2) * mkt.T + mkt.sigma * np.sqrt(mkt.T) * Z[:, 0])
    Y = np.exp(-mkt.r * mkt.T) * np.maximum(ST - 100, 0)
    assert np.isclose(r.price, Y.mean()) and np.isclose(r.se, Y.std(ddof=1) / np.sqrt(1000))


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print(f"PASS {name}")
