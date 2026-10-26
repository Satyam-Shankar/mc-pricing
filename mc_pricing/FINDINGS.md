# Findings

Base case: S₀ = 100, r = 5%, q = 0. Full tables are in `results/summary.md`;
figures are in `figures/`. Total run time of the experiment suite was about
5 minutes on 2 CPU cores.

**1. The estimators are unbiased and their confidence intervals are honest — except for rare-event payoffs.**

Across 90 European contracts priced with 200,000 paths, the standardised
errors z = (MC − Black–Scholes)/SE behave like N(0, 1), as they should:

| Method | Mean z | SD of z | Share outside ±1.96 | Mean relative error* | Max relative error* |
|---|---|---|---|---|---|
| Plain | 0.02 | 1.00 | 5.6% | 0.27% | 1.3% |
| Antithetic + control | 0.11 | 0.96 | 4.4% | 0.04% | 0.18% |

\*Over the 77 contracts priced at 1 or more.

Over 2000 repetitions, the nominal 95% interval covered the true price as follows:

| Case | N = 50 | N = 100 | N = 1,000 | N = 10,000 |
|---|---|---|---|---|
| At-the-money call | 92–94% | 94–95% | 95% | 95% |
| Deep out-of-the-money call (K = 160) | 2–3% | 5% | 40% | 83% |

The deep out-of-the-money failure happens because most samples contain few or
no paths that finish in the money. The sample standard deviation then badly
underestimates the true one, so the interval is far too narrow.
Antithetic and control variates do not fix this. It is the case for
importance sampling.

**2. Convergence follows theory: every Monte Carlo method at N^(−1/2), randomised QMC at about N^(−1).**

Fitted log-log slopes of RMSE against N, from 128 to 1,048,576 paths:

| Method | Slope | RMSE at N = 65,536 |
|---|---|---|
| Plain | −0.503 | 0.059 |
| Antithetic | −0.503 | 0.041 |
| Control | −0.508 | 0.023 |
| Antithetic + control | −0.522 | 0.012 |
| QMC | −1.025 | 0.0007 |
| QMC + control | −1.059 | 0.0003 |

Variance reduction changes the *constant*, not the *rate*:
antithetic + control is worth about 22× more paths of plain Monte Carlo.
Only QMC changes the rate. The reported standard errors matched the realised
RMSE to within 7% for every method.

**3. How much variance reduction helps depends strongly on the contract.**

Over the 131 contracts priced at 1% of spot or more (median, with the
10th–90th percentile range):

| Method | Variance-reduction factor | Time-adjusted efficiency |
|---|---|---|
| Antithetic | 1.4 (1.1–12) | 1.8 |
| Control (stock price) | 13 (3–155) | 12 |
| Antithetic + control | 96 (27–1,252) | 130 |
| QMC | 6,240 | 1,590 |
| QMC + control | 62,000 | 16,900 |

Both antithetic and control variates work best for in-the-money,
low-volatility calls, whose payoff is nearly linear in the stock price. At
K/S₀ = 0.8 and σ = 0.1, antithetic reaches 85× and control 3,032×. Both fall
to about 1× for out-of-the-money options, whose payoff is mostly zero.

At a fixed budget, the relative standard error of plain Monte Carlo grows:
- with volatility: 0.30% → 0.78% as σ goes from 0.1 to 0.8;
- with out-of-the-moneyness: 0.28% → 0.97% as K/S₀ goes from 0.8 to 1.4;
- only mildly with maturity: 0.40% → 0.60% as T goes from 0.1 to 5.

QMC is about 4× slower per path, because of the inverse-normal transform and
scrambling. Its time-adjusted gain is therefore about 4× smaller than its
per-path gain, but still the largest of any method.

**4. Arithmetic Asian options need the right control, and QMC needs a Brownian bridge.**

Reference prices were computed with QMC + control, a Brownian bridge and 2²⁰
paths, giving standard errors of about 3×10⁻⁵. For example, at the money with
m = 52 the price is 5.85384.

Variance-reduction factors, at 12 / 52 / 252 monitoring dates:

| Method | m = 12 | m = 52 | m = 252 |
|---|---|---|---|
| Antithetic | 2.1 | 2.1 | 2.1 |
| Control (geometric Asian) | ~1,300 | ~1,300 | ~1,300 |
| QMC, standard path construction | 334 | 63 | 22 |
| QMC, Brownian bridge | 1,100–1,500 | 1,100–1,500 | 1,100–1,500 |
| QMC + control + bridge | 27,000–59,000 | 27,000–59,000 | 27,000–59,000 |

- The geometric-average Asian option, which has a closed form, is an almost perfectly correlated control, and its benefit does not depend on m.
- QMC with standard path construction degrades as the dimension m grows.
- With a Brownian bridge it stays stable across m.

All 72 estimates were consistent with the reference (mean |z| = 0.90).

**5. For barrier options, discretisation bias dominates, and a continuity correction beats brute force.**

Monitoring a down-and-out call (K = 100, H = 90) at m discrete dates
overprices it relative to continuous monitoring. The bias falls from 1.30 at
m = 4 to 0.17 at m = 512. The fitted rate is m^(−0.42), approaching the
theoretical m^(−1/2).

Shifting the barrier by the Broadie–Glasserman–Kou factor removes almost all
of this bias: below 0.02 for every m ≥ 16.

Under a fixed budget of N × m = 4.2 million path-steps:

| Approach | Best m | RMSE vs continuous price |
|---|---|---|
| Uncorrected | 256 | 0.287 (still bias-dominated) |
| BGK-corrected | 16–32 | 0.026 |

That is 11× lower error at the same cost. With the correction, the best
allocation uses few time steps and many paths. Without it, no affordable
number of time steps gets close.

**6. Choose the time-stepping scheme with care.**

For a European call (σ = 0.4), the Euler–Maruyama scheme has a weak error that
decays like m^(−0.95), against a theoretical m^(−1):

| Time steps m | Euler error |
|---|---|
| 4 | 0.164 |
| 64 | 0.013 |

The exact log-normal scheme has no discretisation bias: one step gave an error
of 0.002, which is just statistical noise. Wherever exact simulation exists,
each time step spent on Euler is a waste of budget.

**Overall cost/accuracy guidance**
1. Use exact simulation where available.
2. Fix the discretisation bias analytically (for example with a continuity correction) rather than by adding time steps.
3. Pick a control variate with a closed form that is highly correlated with the payoff (correlation above 0.99 gives a VRF over 1,000).
4. Add randomised QMC with a Brownian bridge for smooth payoffs.
5. Treat confidence intervals for rare-event payoffs with suspicion.
