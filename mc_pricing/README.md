# Monte Carlo Methods for Derivative Pricing and Numerical Error Analysis

*Independent research project — Probability, Numerical Methods & Quantitative Finance*

A Monte Carlo pricing library for European and path-dependent options under
Black–Scholes dynamics, together with a set of numerical experiments on the
**error** of Monte Carlo estimates:

* how fast the error shrinks with simulation size;
* whether the reported confidence intervals are honest;
* how much variance reduction and quasi-Monte Carlo buy at a fixed computational budget;
* how time discretisation trades off against the number of paths.

Every estimate is validated against an analytical price wherever one exists.

---

## 1. Model and ground truth — `mcpricing/analytic.py`

Under the risk-neutral measure the stock follows geometric Brownian motion,
`dS = (r − q) S dt + σ S dW`, which can be sampled exactly:
`S(t+Δt) = S(t) exp((r − q − σ²/2)Δt + σ√Δt Z)`.

| Contract | Benchmark price |
|---|---|
| European call / put | Black–Scholes–Merton (exact) |
| Geometric-average Asian call (discrete) | Kemna–Vorst-type closed form (exact; log of the average is normal) |
| Arithmetic-average Asian call | no closed form → high-precision QMC reference with its own standard error |
| Down-and-out call, continuous monitoring | Merton / Reiner–Rubinstein (exact) |
| Down-and-out call, discrete monitoring | Broadie–Glasserman–Kou continuity correction (accurate approximation) |

## 2. Estimators — `mcpricing/engine.py`

Every method uses the **same number of simulated paths N**, so comparisons are
at equal simulation budget. Run time is measured too, so efficiency can also
be compared at equal *time*.

| Method | Idea |
|---|---|
| `plain` | mean of i.i.d. discounted payoffs; SE = s/√N; 95% CI = mean ± 1.96 SE |
| `antithetic` | pairs (Z, −Z); averages each pair; works when the payoff is monotone in Z |
| `control` | `Y − b(X − E[X])` with the optimal b = Cov(Y,X)/Var(X) |
| `antithetic+control` | control variate applied to antithetic pair means |
| `qmc` | randomised quasi-Monte Carlo: 16 independent scrambled Sobol' sequences; SE from their spread |
| `qmc+control` | both |

The control variates are:
- **European options:** the discounted terminal stock price, whose expectation is S₀e^(−qT).
- **Arithmetic Asian:** the geometric Asian payoff, which has a closed form.
- **Barrier:** the vanilla European call.

Path-dependent QMC can optionally use a **Brownian-bridge** construction, so
that the best-distributed Sobol' coordinates fix the overall shape of each path.

Simulation runs in memory-bounded batches. Only running sums are stored, so N
can be in the millions with hundreds of monitoring dates.

## 3. Experiments — `experiments/run_experiments.py`

| ID | Question | Design |
|---|---|---|
| **E1** | Are the estimates and their confidence intervals correct? | 90 European contracts vs Black–Scholes (z-scores); empirical 95% CI coverage from 2000 repetitions, N = 50 … 10⁴, at-the-money vs deep out-of-the-money |
| **E2** | Convergence rate | RMSE vs N = 2⁷ … 2²⁰ (30–100 repetitions each), all six estimators, with fitted log-log slopes |
| **E3** | When does variance reduction help, and does it pay for its run time? | 150 contracts: K/S₀ ∈ [0.8, 1.4], σ ∈ [0.1, 0.8], T ∈ [0.1, 5]; variance-reduction factor and time-adjusted efficiency |
| **E4** | Path-dependent: arithmetic Asian | all estimators (QMC with and without Brownian bridge) for m = 12, 52, 252 monitoring dates |
| **E5** | Path-dependent: barrier | (a) discretisation bias of discrete monitoring vs m, with and without the BGK barrier shift; (b) optimal split of a fixed budget N·m between paths and time steps |
| **E6** | Time-stepping scheme | Euler–Maruyama vs exact scheme: weak error vs number of steps |

Definitions used throughout:
- **Variance-reduction factor (VRF)** = Var_plain / Var_method per path.
- **Time-adjusted efficiency** = (SE² × time)_plain / (SE² × time)_method. A value of 10 means the method reaches a given accuracy 10× faster.

## 4. Running

Run these one at a time (copy-paste safe in bash and zsh):

```bash
pip install -r requirements.txt
python tests/test_pricing.py
python demo.py
python experiments/run_experiments.py
python experiments/analyze.py
```

| Command | What it does | Time |
|---|---|---|
| `tests/test_pricing.py` | Black–Scholes reference value and put–call parity, formula limits, all estimators unbiased against closed forms, Brownian-bridge distribution, SE definition | ~2 s |
| `demo.py` | prices a European, an Asian and a barrier option with every method | ~3 s |
| `experiments/run_experiments.py` | E1–E6, writes `results/*.csv`; add `--quick` for a 30-second smoke test or `--only e2 e5` for a subset | ~5–10 min |
| `experiments/analyze.py` | writes `figures/` and `results/summary.md` | ~5 s |

## 5. Layout

```
mcpricing/
  analytic.py   Black-Scholes, geometric Asian, continuous & BGK barrier formulas
  simulate.py   Market, pseudo-random / Sobol' normals, Brownian bridge, exact & Euler paths
  options.py    European, Asian (arithmetic, geometric), down-and-out contracts + control variates
  engine.py     batched estimators: plain, antithetic, control, QMC, combinations
experiments/
  run_experiments.py  analyze.py
tests/test_pricing.py
demo.py
```

## 6. Findings

See `FINDINGS.md`, `results/summary.md` and `figures/`.

## References
* Glasserman, P. (2003). *Monte Carlo Methods in Financial Engineering*. Springer.
* Boyle, P., Broadie, M. & Glasserman, P. (1997). Monte Carlo methods for security pricing. *J. Econ. Dynamics & Control*.
* Kemna, A. & Vorst, A. (1990). A pricing method for options based on average asset values. *J. Banking & Finance*.
* Broadie, M., Glasserman, P. & Kou, S. (1997). A continuity correction for discrete barrier options. *Mathematical Finance*.
* Owen, A. B. (1997). Scrambled net variance for integrals of smooth functions. *Annals of Statistics*.
* Kloeden, P. & Platen, E. (1992). *Numerical Solution of Stochastic Differential Equations*. Springer.
