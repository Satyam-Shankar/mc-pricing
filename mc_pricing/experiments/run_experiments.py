"""
Numerical experiments.

E1 validation   : MC vs Black-Scholes over 90 contracts (z-scores) and the
                  empirical coverage of 95% confidence intervals vs N
E2 convergence  : RMSE vs number of paths for every estimator (rate fitting)
E3 efficiency   : variance-reduction factor and time-adjusted efficiency over a
                  grid of moneyness x volatility x maturity at a fixed budget
E4 asian        : arithmetic Asian call - estimator comparison vs monitoring dates
E5 barrier      : down-and-out call - discretisation bias vs monitoring
                  frequency, continuity correction, and the optimal split of a
                  fixed budget N*m between paths and time steps
E6 euler        : Euler vs exact scheme for a European call - weak-error bias
                  vs time steps and cost-accuracy at fixed budget

    python experiments/run_experiments.py [--quick] [--only e2 e5]
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT)

from mcpricing import (METHODS, ArithmeticAsianCall, DownAndOutCall, EuropeanCall,  # noqa: E402
                       EuropeanPut, Market, price)
from mcpricing.analytic import BGK_BETA  # noqa: E402

OUT = os.path.join(ROOT, "results")
Q = False  # quick mode


def save(df, name):
    os.makedirs(OUT, exist_ok=True)
    df.to_csv(os.path.join(OUT, name + ".csv"), index=False)
    print(f"  -> results/{name}.csv ({len(df)} rows)", flush=True)


# --------------------------------------------------------------------------- E1
def e1_validation():
    rows = []
    Ks = [80, 90, 100, 110, 120]
    sigmas = [0.1, 0.2, 0.4]
    Ts = [0.25, 1.0, 2.0]
    N = 20_000 if Q else 200_000
    seed = 0
    for K in Ks:
        for s in sigmas:
            for T in Ts:
                for opt in (EuropeanCall(K), EuropeanPut(K)):
                    mkt = Market(sigma=s, T=T)
                    exact = opt.exact(mkt)
                    for meth in ("plain", "antithetic+control"):
                        seed += 1
                        r = price(opt, mkt, N, meth, seed=seed)
                        rows.append(dict(option=opt.name, K=K, sigma=s, T=T, method=meth, N=N,
                                         exact=exact, mc=r.price, se=r.se,
                                         z=(r.price - exact) / r.se if r.se > 0 else 0.0))
    save(pd.DataFrame(rows), "e1_validation")

    # CI coverage: repeat the experiment many times and count how often the
    # 95% interval contains the true price
    rows = []
    R = 200 if Q else 2000
    cases = [("ATM call", EuropeanCall(100), Market()),
             ("deep OTM call", EuropeanCall(160), Market(sigma=0.2, T=0.5))]
    for label, opt, mkt in cases:
        exact = opt.exact(mkt)
        for N in [50, 100, 500, 1000, 10_000]:
            for meth in ("plain", "antithetic", "control"):
                hits = []
                for rep in range(R):
                    r = price(opt, mkt, N, meth, seed=10_000 + rep)
                    lo, hi = r.ci()
                    hits.append(lo <= exact <= hi)
                rows.append(dict(case=label, N=N, method=meth, coverage=np.mean(hits), reps=R,
                                 exact=exact))
    save(pd.DataFrame(rows), "e1_coverage")


# --------------------------------------------------------------------------- E2
def e2_convergence():
    mkt = Market()
    opt = EuropeanCall(100)
    exact = opt.exact(mkt)
    rows = []
    exps = range(7, 17) if Q else range(7, 21)
    for e in exps:
        N = 2 ** e
        R = 20 if Q else (100 if N <= 2 ** 16 else 30)
        for meth in METHODS:
            for rep in range(R):
                r = price(opt, mkt, N, meth, seed=rep + 1)
                rows.append(dict(N=N, method=meth, rep=rep, price=r.price, se=r.se,
                                 err=r.price - exact, time_s=r.time_s))
    save(pd.DataFrame(rows), "e2_convergence")


# --------------------------------------------------------------------------- E3
def e3_efficiency():
    rows = []
    N = 2 ** 14 if Q else 2 ** 17
    seed = 0
    for K in [80, 90, 100, 110, 120, 140]:
        for s in [0.1, 0.2, 0.3, 0.5, 0.8]:
            for T in [0.1, 0.5, 1.0, 2.0, 5.0]:
                mkt = Market(sigma=s, T=T)
                opt = EuropeanCall(K)
                exact = opt.exact(mkt)
                seed += 1
                for meth in METHODS:
                    r = price(opt, mkt, N, meth, seed=seed)
                    rows.append(dict(K=K, moneyness=K / mkt.S0, sigma=s, T=T, method=meth, N=N,
                                     exact=exact, price=r.price, se=r.se, time_s=r.time_s,
                                     rho=r.extra.get("rho", np.nan)))
    save(pd.DataFrame(rows), "e3_efficiency")


# --------------------------------------------------------------------------- E4
def e4_asian():
    mkt = Market()
    rows = []
    N = 2 ** 13 if Q else 2 ** 16
    for m in [12, 52, 252]:
        for K in [90, 100, 110]:
            opt = ArithmeticAsianCall(K)
            ref = price(opt, mkt, 2 ** 16 if Q else 2 ** 20, "qmc+control", m=m, seed=999, bridge=True)
            for meth in METHODS:
                for bridge in ([False, True] if meth.startswith("qmc") else [False]):
                    r = price(opt, mkt, N, meth, m=m, seed=7, bridge=bridge)
                    rows.append(dict(m=m, K=K, method=meth + (" (bridge)" if bridge else ""),
                                     N=N, price=r.price, se=r.se, time_s=r.time_s,
                                     ref=ref.price, ref_se=ref.se,
                                     rho=r.extra.get("rho", np.nan)))
            print(f"  asian m={m} K={K} ref={ref.price:.5f} ± {ref.se:.5f}", flush=True)
    save(pd.DataFrame(rows), "e4_asian")


# --------------------------------------------------------------------------- E5
def e5_barrier():
    mkt = Market()
    K, H = 100, 90
    cont = DownAndOutCall(K, H).continuous(mkt)
    rows = []
    # (a) discretisation bias: discrete-monitoring MC price vs continuous price
    N = 2 ** 15 if Q else 2 ** 18
    for m in [4, 8, 16, 32, 64, 128, 256, 512]:
        for corrected in (False, True):
            Hm = H * np.exp(BGK_BETA * mkt.sigma * np.sqrt(mkt.T / m)) if corrected else H
            r = price(DownAndOutCall(K, Hm), mkt, N, "antithetic+control", m=m, seed=m)
            rows.append(dict(part="bias", m=m, corrected=corrected, N=N, price=r.price, se=r.se,
                             continuous=cont, bgk_discrete=DownAndOutCall(K, H).exact(mkt, m),
                             time_s=r.time_s))
        print(f"  barrier bias m={m}", flush=True)
    # (b) fixed budget N*m: how to split between paths and steps?
    B = 2 ** 20 if Q else 2 ** 22
    R = 10 if Q else 25
    for m in [4, 8, 16, 32, 64, 128, 256]:
        Nm = B // m
        for corrected in (False, True):
            Hm = H * np.exp(BGK_BETA * mkt.sigma * np.sqrt(mkt.T / m)) if corrected else H
            for rep in range(R):
                r = price(DownAndOutCall(K, Hm), mkt, Nm, "plain", m=m, seed=100 + rep)
                rows.append(dict(part="budget", m=m, corrected=corrected, N=Nm, rep=rep,
                                 price=r.price, se=r.se, continuous=cont, time_s=r.time_s))
        print(f"  barrier budget m={m}", flush=True)
    save(pd.DataFrame(rows), "e5_barrier")


# --------------------------------------------------------------------------- E6
def e6_euler():
    mkt = Market(sigma=0.4)
    opt = EuropeanCall(100)
    exact = opt.exact(mkt)
    rows = []
    N = 2 ** 16 if Q else 2 ** 20
    for m in [1, 2, 4, 8, 16, 32, 64]:
        for scheme in ("euler", "exact"):
            # common random numbers would hide nothing here; use antithetic+control
            # so the statistical error is small relative to the bias
            r = price(opt, mkt, N, "antithetic+control", m=m, seed=3, scheme=scheme)
            rows.append(dict(m=m, scheme=scheme, N=N, price=r.price, se=r.se, exact=exact,
                             time_s=r.time_s))
        print(f"  euler m={m}", flush=True)
    save(pd.DataFrame(rows), "e6_euler")


EXPS = {"e1": e1_validation, "e2": e2_convergence, "e3": e3_efficiency,
        "e4": e4_asian, "e5": e5_barrier, "e6": e6_euler}


def main():
    global Q
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--only", nargs="*")
    args = ap.parse_args()
    Q = args.quick
    for k, f in EXPS.items():
        if args.only and k not in args.only:
            continue
        t0 = time.time()
        print(f"[{k}] {f.__doc__ or f.__name__}", flush=True)
        f()
        print(f"[{k}] done in {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
