"""
Figures and summary tables from results/*.csv.

    python experiments/analyze.py
"""
from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import norm  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RES, FIG = os.path.join(ROOT, "results"), os.path.join(ROOT, "figures")
METHODS = ["plain", "antithetic", "control", "antithetic+control", "qmc", "qmc+control"]
COL = {"plain": "#6b7785", "antithetic": "#4f8fc0", "control": "#e08a2b",
       "antithetic+control": "#b5541c", "qmc": "#3a9a6a", "qmc+control": "#1e5e3c",
       "qmc (bridge)": "#7cc39b", "qmc+control (bridge)": "#2f7d55"}
plt.rcParams.update({"figure.dpi": 130, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.alpha": 0.25, "font.size": 9.5,
                     "legend.frameon": False})


def rd(name):
    p = os.path.join(RES, name + ".csv")
    return pd.read_csv(p) if os.path.exists(p) else None


def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, name), bbox_inches="tight")
    plt.close(fig)


def slope(x, y):
    b, a = np.polyfit(np.log(x), np.log(y), 1)
    return b, np.exp(a)


md = ["# Results summary", ""]


# ------------------------------------------------------------------- E1
def e1():
    global md
    v, c = rd("e1_validation"), rd("e1_coverage")
    if v is None:
        return
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    p = v[v.method == "plain"]
    ax[0].errorbar(p.exact, p.mc, yerr=1.96 * p.se, fmt="o", ms=3, color="#4f8fc0", alpha=0.8)
    lim = [0, p.exact.max() * 1.05]
    ax[0].plot(lim, lim, "k-", lw=0.8)
    ax[0].set_xlabel("Black-Scholes price")
    ax[0].set_ylabel("Monte Carlo price (±95% CI)")
    ax[0].set_title(f"90 calls & puts, N = {int(p.N.iloc[0]):,}")
    xs = np.linspace(-4, 4, 200)
    for meth, col in (("plain", "#4f8fc0"), ("antithetic+control", "#b5541c")):
        ax[1].hist(v[v.method == meth].z, bins=np.linspace(-4, 4, 17), density=True, alpha=0.5,
                   color=col, label=meth)
    ax[1].plot(xs, norm.pdf(xs), "k-", lw=1, label="N(0,1)")
    ax[1].set_xlabel("z = (MC - BS) / SE")
    ax[1].set_title("Standardised errors")
    ax[1].legend(fontsize=8)
    if c is not None:
        for (case, meth), g in c.groupby(["case", "method"]):
            ax[2].plot(g.N, g.coverage, marker="o", ms=3, color=COL[meth],
                       ls="-" if case.startswith("ATM") else "--", label=f"{case}, {meth}")
        R = int(c.reps.iloc[0])
        band = 1.96 * np.sqrt(0.95 * 0.05 / R)
        ax[2].axhspan(0.95 - band, 0.95 + band, color="k", alpha=0.08)
        ax[2].axhline(0.95, color="k", lw=0.8)
        ax[2].set_xscale("log")
        ax[2].set_xlabel("paths N")
        ax[2].set_ylabel("empirical coverage")
        ax[2].set_title(f"95% CI coverage ({R} repetitions)")
        ax[2].legend(fontsize=6.5, ncol=1)
    save(fig, "fig1_validation.png")

    md.extend(["## E1 - validation against Black-Scholes", ""])
    s = v.groupby("method").z.agg(mean="mean", sd="std",
                                  frac_outside_95=lambda z: np.mean(np.abs(z) > 1.96),
                                  max_abs="max")
    s["max_abs"] = v.groupby("method").z.apply(lambda z: np.abs(z).max())
    md += [f"{len(v) // 2} contracts (K = 80-120, sigma = 0.1-0.4, T = 0.25-2, calls and puts). "
           "If the estimator is unbiased with correct standard errors, z ~ N(0,1): mean 0, sd 1, "
           "about 5% outside ±1.96.", "", s.to_markdown(floatfmt=".3f"), ""]
    vv = v[v.exact >= 1.0]
    rel = vv.assign(rel=np.abs(vv.mc - vv.exact) / vv.exact).groupby("method").rel.agg(["mean", "max"])
    md += [f"Absolute relative pricing error ({len(vv) // 2} contracts with price >= 1; "
           "for near-worthless options relative error is not meaningful):", "", (100 * rel).rename(
        columns={"mean": "mean %", "max": "max %"}).to_markdown(floatfmt=".3f"), ""]
    if c is not None:
        t = c.pivot_table(index=["case", "method"], columns="N", values="coverage")
        md += ["### Empirical coverage of the nominal 95% confidence interval", "",
               t.to_markdown(floatfmt=".3f"), ""]


# ------------------------------------------------------------------- E2
def e2():
    global md
    d = rd("e2_convergence")
    if d is None:
        return
    g = d.groupby(["method", "N"]).agg(rmse=("err", lambda e: np.sqrt(np.mean(e ** 2))),
                                       se=("se", "mean"), t=("time_s", "mean")).reset_index()
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.2))
    rows = []
    for meth in METHODS:
        h = g[g.method == meth]
        k, a = slope(h.N, h.rmse)
        rows.append(dict(method=meth, rate=k, rmse_at_2e16=float(
            h.loc[h.N == 2 ** 16, "rmse"].iloc[0]) if (h.N == 2 ** 16).any() else np.nan))
        ax[0].plot(h.N, h.rmse, "o-", ms=3, color=COL[meth], label=f"{meth}  (slope {k:.2f})")
        ax[1].plot(h.t * 1000, h.rmse, "o-", ms=3, color=COL[meth], label=meth)
    Ns = np.array([g.N.min(), g.N.max()], float)
    base = g[(g.method == "plain") & (g.N == g.N.min())].rmse.iloc[0]
    ax[0].plot(Ns, base * (Ns / Ns[0]) ** -0.5, "k:", lw=1, label="reference $N^{-1/2}$")
    ax[0].plot(Ns, base * (Ns / Ns[0]) ** -1.0, "k--", lw=0.8, label="reference $N^{-1}$")
    for a in ax:
        a.set_xscale("log")
        a.set_yscale("log")
        a.set_ylabel("RMSE vs Black-Scholes")
    ax[0].set_xlabel("paths N")
    ax[0].set_title("Convergence rate (ATM call, S=K=100, σ=0.2, T=1)")
    ax[1].set_xlabel("wall-clock time per estimate (ms)")
    ax[1].set_title("Error vs computational cost")
    ax[0].legend(fontsize=7.5)
    save(fig, "fig2_convergence.png")
    md.extend(["## E2 - convergence", "",
               "Fitted log-log slope of RMSE vs N (theory: -0.5 for Monte Carlo, "
               "close to -1 for randomised QMC on smooth problems).", "",
               pd.DataFrame(rows).to_markdown(index=False, floatfmt=".4f"), ""])
    # CI calibration: mean reported SE vs realised RMSE
    cal = g.assign(ratio=g.se / g.rmse).groupby("method").ratio.median()
    md += ["Median ratio (reported SE) / (realised RMSE) across N - close to 1 means the "
           "standard errors are honest:", "", cal.to_frame("SE / RMSE").to_markdown(floatfmt=".2f"), ""]


# ------------------------------------------------------------------- E3
def e3():
    global md
    d = rd("e3_efficiency")
    if d is None:
        return
    d["var"] = d.se ** 2 * d.N
    key = ["K", "sigma", "T"]
    base = d[d.method == "plain"].set_index(key)[["var", "se", "time_s"]]
    d = d.join(base, on=key, rsuffix="_plain")
    d["vrf"] = d.var_plain / d["var"]
    d["eff"] = (d.se_plain ** 2 * d.time_s_plain) / (d.se ** 2 * d.time_s)
    d["rel_se_plain"] = d.se_plain / d.exact
    meths = ["antithetic", "control", "antithetic+control", "qmc", "qmc+control"]

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    sub = d[d["T"] == 1.0]
    for ax, meth in zip(axes, ["antithetic", "control", "qmc+control"]):
        h = sub[sub.method == meth].pivot(index="sigma", columns="moneyness", values="vrf")
        im = ax.imshow(np.log10(h.values), cmap="viridis", origin="lower", aspect="auto")
        ax.set_xticks(range(h.shape[1]), [f"{x:.1f}" for x in h.columns])
        ax.set_yticks(range(h.shape[0]), h.index)
        for i in range(h.shape[0]):
            for j in range(h.shape[1]):
                val = h.values[i, j]
                ax.text(j, i, f"{val:.0f}" if val >= 10 else f"{val:.1f}", ha="center",
                        va="center", fontsize=7.5, color="k" if (np.log10(val) - np.log10(h.values).min()) > 0.6 * np.ptp(np.log10(h.values)) else "w")
        ax.set_xlabel("moneyness K / S0")
        ax.set_ylabel("volatility σ")
        ax.set_title(f"{meth}: variance-reduction factor")
        ax.grid(False)
        fig.colorbar(im, ax=ax, label="log10 VRF", shrink=0.8)
    fig.suptitle("Variance reduction vs plain MC (European call, T = 1)", y=1.02)
    save(fig, "fig3_vrf_heatmaps.png")

    fig, ax = plt.subplots(1, 3, figsize=(14, 3.9))
    for meth in ["plain"] + meths:
        h = d[d.method == meth]
        for a, x, flt in ((ax[0], "moneyness", (h.sigma == 0.2) & (h["T"] == 1.0)),
                          (ax[1], "sigma", (h.K == 100) & (h["T"] == 1.0)),
                          (ax[2], "T", (h.K == 100) & (h.sigma == 0.2))):
            q = h[flt].sort_values(x)
            a.plot(q[x], q.se / q.exact * 100, "o-", ms=3, color=COL[meth], label=meth)
    for a, x in zip(ax, ["moneyness K/S0  (σ=0.2, T=1)", "volatility σ  (ATM, T=1)",
                          "maturity T  (ATM, σ=0.2)"]):
        a.set_xlabel(x)
        a.set_yscale("log")
        a.set_ylabel("relative standard error (%)")
    ax[1].set_xscale("linear")
    ax[2].set_xscale("log")
    ax[0].legend(fontsize=7.5)
    fig.suptitle(f"Relative estimation error at fixed budget N = {int(d.N.iloc[0]):,}", y=1.02)
    save(fig, "fig4_error_sensitivity.png")

    # exclude contracts that are essentially worthless (plain MC sees almost
    # no positive payoffs, so its variance - the denominator - is meaningless)
    ok = (d.exact >= 0.01 * d.K / d.moneyness) & (d.se_plain > 0)
    s = d[ok & (d.method != "plain")].groupby("method")[["vrf", "eff"]].quantile(
        [0.1, 0.5, 0.9]).unstack()
    s.columns = [f"{a} p{int(100 * b)}" for a, b in s.columns]
    md += ["## E3 - efficiency over 150 contracts (K/S0 0.8-1.4, σ 0.1-0.8, T 0.1-5)", "",
           f"{int(ok.sum() / 6)} contracts with price >= 1% of spot; 10th / 50th / 90th percentiles.", "",
           "VRF = Var(plain)/Var(method) per path; efficiency also charges for run time: "
           "(SE² × time)_plain / (SE² × time)_method.", "",
           s.reindex(meths).to_markdown(floatfmt=".1f"), ""]
    by_m = d[d.method.isin(["antithetic", "control"])].pivot_table(
        index="moneyness", columns="method", values="vrf", aggfunc="median")
    md += ["Median VRF by moneyness (all σ, T):", "", by_m.to_markdown(floatfmt=".2f"), ""]
    p = d[d.method == "plain"]
    md += ["Relative SE of plain MC (%, median over the other two dimensions):", "",
           (100 * p.groupby("sigma").rel_se_plain.median()).to_frame("by σ").T.to_markdown(floatfmt=".2f"), "",
           (100 * p.groupby("T").rel_se_plain.median()).to_frame("by T").T.to_markdown(floatfmt=".2f"), "",
           (100 * p.groupby("moneyness").rel_se_plain.median()).to_frame("by K/S0").T.to_markdown(floatfmt=".2f"), ""]


# ------------------------------------------------------------------- E4
def e4():
    global md
    d = rd("e4_asian")
    if d is None:
        return
    d["var"] = d.se ** 2 * d.N
    base = d[d.method == "plain"].set_index(["m", "K"])[["var", "time_s", "se"]]
    d = d.join(base, on=["m", "K"], rsuffix="_p")
    d["vrf"] = d.var_p / d["var"]
    d["eff"] = (d.se_p ** 2 * d.time_s_p) / (d.se ** 2 * d.time_s)
    d["z"] = (d.price - d.ref) / np.sqrt(d.se ** 2 + d.ref_se ** 2)
    order = ["plain", "antithetic", "control", "antithetic+control", "qmc", "qmc (bridge)",
             "qmc+control", "qmc+control (bridge)"]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    h = d[d.K == 100]
    w = 0.25
    for i, m in enumerate(sorted(h.m.unique())):
        q = h[h.m == m].set_index("method").reindex(order)
        x = np.arange(len(order))
        ax[0].bar(x + (i - 1) * w, q.vrf, w, label=f"m = {m} dates")
        ax[1].bar(x + (i - 1) * w, q.eff, w, label=f"m = {m} dates")
    for a, t in zip(ax, ["variance-reduction factor", "time-adjusted efficiency gain"]):
        a.set_xticks(np.arange(len(order)), order, rotation=35, ha="right")
        a.set_yscale("log")
        a.set_title(f"Arithmetic Asian call (ATM): {t}")
        a.axhline(1, color="k", lw=0.8)
    ax[0].legend(fontsize=8)
    save(fig, "fig5_asian.png")
    t = d.pivot_table(index="method", columns="m", values="vrf", aggfunc="median").reindex(order)
    refs = d.groupby(["m", "K"])[["ref", "ref_se"]].first().reset_index()
    md += ["## E4 - arithmetic Asian call (no closed form)", "",
           "Reference prices (QMC + control variate, Brownian bridge, 2^20 paths):", "",
           refs.to_markdown(index=False, floatfmt=".5f"), "",
           "Median variance-reduction factor over K in {90, 100, 110}:", "",
           t.to_markdown(floatfmt=".1f"), "",
           f"Consistency with the reference: mean |z| = {np.abs(d.z).mean():.2f}, "
           f"max |z| = {np.abs(d.z).max():.2f} over {len(d)} estimates.", ""]


# ------------------------------------------------------------------- E5
def e5():
    global md
    d = rd("e5_barrier")
    if d is None:
        return
    b = d[d.part == "bias"]
    cont = b.continuous.iloc[0]
    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    for corr, col, lab in ((False, "#4f8fc0", "discrete monitoring, barrier H"),
                           (True, "#b5541c", "barrier shifted by BGK correction")):
        q = b[b.corrected == corr]
        ax[0].errorbar(q.m, q.price, yerr=1.96 * q.se, marker="o", ms=3, color=col, label=lab, capsize=2)
        ax[1].plot(q.m, np.abs(q.price - cont), "o-", ms=3, color=col, label=lab)
    q = b[~b.corrected]
    ax[0].plot(q.m, q.bgk_discrete, "k:", lw=1, label="BGK analytic, discrete")
    ax[0].axhline(cont, color="k", lw=0.9, label="continuous (analytic)")
    k, a = slope(q.m, np.abs(q.price - cont))
    ax[1].plot(q.m, a * q.m ** k, "k:", lw=0.8, label=f"fit slope {k:.2f}")
    for a_ in ax[:2]:
        a_.set_xscale("log", base=2)
        a_.set_xlabel("monitoring dates m")
    ax[0].set_ylabel("down-and-out call price")
    ax[0].set_title("Price vs monitoring frequency")
    ax[1].set_yscale("log")
    ax[1].set_ylabel("|MC price - continuous price|")
    ax[1].set_title("Discretisation bias")
    ax[0].legend(fontsize=7)
    ax[1].legend(fontsize=7)
    u = d[d.part == "budget"]
    g = u.groupby(["m", "corrected"]).apply(
        lambda x: pd.Series(dict(rmse=np.sqrt(np.mean((x.price - x.continuous) ** 2)),
                                 bias=abs(x.price.mean() - x.continuous.iloc[0]),
                                 sd=x.price.std(), N=x.N.iloc[0])), include_groups=False).reset_index()
    for corr, col in ((False, "#4f8fc0"), (True, "#b5541c")):
        q = g[g.corrected == corr]
        ax[2].plot(q.m, q.rmse, "o-", color=col, ms=3, label=("corrected" if corr else "uncorrected") + " RMSE")
        ax[2].plot(q.m, q.bias, ":", color=col, lw=1, label=("corrected" if corr else "uncorrected") + " |bias|")
    ax[2].plot(g[~g.corrected].m, g[~g.corrected].sd, "--", color="#6b7785", lw=1, label="statistical sd")
    ax[2].set_xscale("log", base=2)
    ax[2].set_yscale("log")
    ax[2].set_xlabel("monitoring dates m   (paths N = budget / m)")
    ax[2].set_ylabel("error vs continuous price")
    ax[2].set_title(f"Fixed budget N × m = {int(u.N.iloc[0] * u.m.iloc[0]):,}")
    ax[2].legend(fontsize=7)
    save(fig, "fig6_barrier.png")
    md += ["## E5 - down-and-out call (K=100, H=90): discretisation and budget allocation", "",
           f"Continuous-monitoring analytic price: {cont:.5f}. Bias of discrete monitoring "
           f"decays like m^{k:.2f} (theory m^-0.5).", "",
           b[["m", "corrected", "price", "se", "bgk_discrete"]].assign(
               bias=b.price - cont).to_markdown(index=False, floatfmt=".5f"), "",
           "Fixed budget N × m - RMSE against the continuous price:", "",
           g.pivot(index="m", columns="corrected", values="rmse").rename(
               columns={False: "uncorrected", True: "BGK-corrected"}).to_markdown(floatfmt=".5f"), ""]


# ------------------------------------------------------------------- E6
def e6():
    global md
    d = rd("e6_euler")
    if d is None:
        return
    exact = d.exact.iloc[0]
    fig, ax = plt.subplots(figsize=(6, 4))
    for scheme, col in (("euler", "#4f8fc0"), ("exact", "#b5541c")):
        q = d[d.scheme == scheme]
        ax.errorbar(q.m, np.abs(q.price - exact), yerr=1.96 * q.se, marker="o", ms=3,
                    color=col, label=f"{scheme} scheme", capsize=2)
    q = d[d.scheme == "euler"]
    qf = q[q.m >= 4]  # m = 1, 2 are pre-asymptotic (the error changes sign)
    k, a = slope(qf.m, np.abs(qf.price - exact))
    ax.plot(qf.m, a * qf.m.astype(float) ** k, "k:", lw=0.9, label=f"Euler fit slope {k:.2f} (m ≥ 4)")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("time steps m")
    ax.set_ylabel("|MC price - Black-Scholes|")
    ax.set_title("Weak discretisation error (European call, σ = 0.4)")
    ax.legend(fontsize=8)
    save(fig, "fig7_euler.png")
    md += ["## E6 - time-stepping scheme", "",
           f"Euler weak error decays like m^{k:.2f} for m >= 4 (theory m^-1); the exact log scheme has no "
           "discretisation bias, so one step suffices for European payoffs.", "",
           d.assign(error=d.price - exact).pivot(index="m", columns="scheme", values="error")
           .to_markdown(floatfmt=".5f"), ""]


def main():
    os.makedirs(FIG, exist_ok=True)
    for f in (e1, e2, e3, e4, e5, e6):
        f()
    open(os.path.join(RES, "summary.md"), "w").write("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
