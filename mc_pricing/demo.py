"""
Quick tour: price an at-the-money European call, an arithmetic Asian call and
a down-and-out call with every estimator at the same budget.

    python demo.py
"""
from mcpricing import (METHODS, ArithmeticAsianCall, DownAndOutCall, EuropeanCall,
                       Market, price)

mkt = Market(S0=100, r=0.05, sigma=0.2, T=1.0)
N = 2 ** 16

cases = [
    (EuropeanCall(K=100), 1, "Black-Scholes"),
    (ArithmeticAsianCall(K=100), 52, "no closed form"),
    (DownAndOutCall(K=100, H=90), 52, "BGK discrete-monitoring formula"),
]
for opt, m, ref_name in cases:
    ref = opt.exact(mkt, m)
    ref_txt = f"{ref:.5f} ({ref_name})" if ref is not None else ref_name
    print(f"\n{opt.name}, m = {m} monitoring dates, N = {N:,} paths   reference: {ref_txt}")
    print(f"  {'method':<20}{'price':>10}{'std err':>10}{'95% CI':>22}{'ms':>8}{'VRF':>9}")
    base = None
    for meth in METHODS:
        r = price(opt, mkt, N, meth, m=m, seed=1, bridge=meth.startswith("qmc"))
        base = base or r.se
        lo, hi = r.ci()
        print(f"  {meth:<20}{r.price:>10.5f}{r.se:>10.5f}   [{lo:8.4f}, {hi:8.4f}]"
              f"{1000 * r.time_s:>8.1f}{(base / r.se) ** 2:>9.1f}")
