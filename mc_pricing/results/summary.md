# Results summary

## E1 - validation against Black-Scholes

90 contracts (K = 80-120, sigma = 0.1-0.4, T = 0.25-2, calls and puts). If the estimator is unbiased with correct standard errors, z ~ N(0,1): mean 0, sd 1, about 5% outside ±1.96.

| method             |   mean |    sd |   frac_outside_95 |   max_abs |
|:-------------------|-------:|------:|------------------:|----------:|
| antithetic+control |  0.111 | 0.964 |             0.044 |     2.517 |
| plain              |  0.021 | 1.003 |             0.056 |     2.441 |

Absolute relative pricing error (77 contracts with price >= 1; for near-worthless options relative error is not meaningful):

| method             |   mean % |   max % |
|:-------------------|---------:|--------:|
| antithetic+control |    0.042 |   0.182 |
| plain              |    0.266 |   1.324 |

### Empirical coverage of the nominal 95% confidence interval

|                                 |    50 |   100 |   500 |   1000 |   10000 |
|:--------------------------------|------:|------:|------:|-------:|--------:|
| ('ATM call', 'antithetic')      | 0.923 | 0.940 | 0.949 |  0.952 |   0.950 |
| ('ATM call', 'control')         | 0.915 | 0.948 | 0.953 |  0.956 |   0.950 |
| ('ATM call', 'plain')           | 0.935 | 0.947 | 0.947 |  0.952 |   0.945 |
| ('deep OTM call', 'antithetic') | 0.030 | 0.055 | 0.252 |  0.399 |   0.841 |
| ('deep OTM call', 'control')    | 0.024 | 0.053 | 0.241 |  0.400 |   0.831 |
| ('deep OTM call', 'plain')      | 0.024 | 0.054 | 0.241 |  0.401 |   0.833 |

## E2 - convergence

Fitted log-log slope of RMSE vs N (theory: -0.5 for Monte Carlo, close to -1 for randomised QMC on smooth problems).

| method             |    rate |   rmse_at_2e16 |
|:-------------------|--------:|---------------:|
| plain              | -0.5033 |         0.0589 |
| antithetic         | -0.5031 |         0.0409 |
| control            | -0.5081 |         0.0228 |
| antithetic+control | -0.5224 |         0.0124 |
| qmc                | -1.0248 |         0.0007 |
| qmc+control        | -1.0593 |         0.0003 |

Median ratio (reported SE) / (realised RMSE) across N - close to 1 means the standard errors are honest:

| method             |   SE / RMSE |
|:-------------------|------------:|
| antithetic         |        1.04 |
| antithetic+control |        0.97 |
| control            |        1.04 |
| plain              |        0.98 |
| qmc                |        0.93 |
| qmc+control        |        0.94 |

## E3 - efficiency over 150 contracts (K/S0 0.8-1.4, σ 0.1-0.8, T 0.1-5)

131 contracts with price >= 1% of spot; 10th / 50th / 90th percentiles.

VRF = Var(plain)/Var(method) per path; efficiency also charges for run time: (SE² × time)_plain / (SE² × time)_method.

| method             |   vrf p10 |   vrf p50 |   vrf p90 |   eff p10 |   eff p50 |   eff p90 |
|:-------------------|----------:|----------:|----------:|----------:|----------:|----------:|
| antithetic         |       1.1 |       1.4 |      12.1 |       1.3 |       1.6 |      14.9 |
| control            |       3.1 |      12.8 |     155.0 |       2.9 |      12.4 |     147.0 |
| antithetic+control |      26.6 |      95.5 |    1252.0 |      29.9 |     104.3 |    1385.2 |
| qmc                |      92.0 |    6239.9 |   59024.3 |      34.4 |    2472.7 |   20675.3 |
| qmc+control        |    6776.8 |   62052.1 |  536930.0 |    2493.2 |   22005.0 |  189904.1 |

Median VRF by moneyness (all σ, T):

|   moneyness |   antithetic |   control |
|------------:|-------------:|----------:|
|        0.80 |         4.05 |     57.85 |
|        0.90 |         2.35 |     20.03 |
|        1.00 |         1.68 |      8.82 |
|        1.10 |         1.27 |      4.68 |
|        1.20 |         1.17 |      3.50 |
|        1.40 |         1.09 |      2.60 |

Relative SE of plain MC (%, median over the other two dimensions):

|      |   0.1 |   0.2 |   0.3 |   0.5 |   0.8 |
|:-----|------:|------:|------:|------:|------:|
| by σ |  0.30 |  0.40 |  0.51 |  0.65 |  0.78 |

|      |   0.1 |   0.5 |   1.0 |   2.0 |   5.0 |
|:-----|------:|------:|------:|------:|------:|
| by T |  0.40 |  0.54 |  0.58 |  0.59 |  0.60 |

|         |   0.8 |   0.9 |   1.0 |   1.1 |   1.2 |   1.4 |
|:--------|------:|------:|------:|------:|------:|------:|
| by K/S0 |  0.28 |  0.34 |  0.44 |  0.60 |  0.73 |  0.97 |

## E4 - arithmetic Asian call (no closed form)

Reference prices (QMC + control variate, Brownian bridge, 2^20 paths):

|         m |         K |      ref |   ref_se |
|----------:|----------:|---------:|---------:|
|  12.00000 |  90.00000 | 12.91990 |  0.00002 |
|  12.00000 | 100.00000 |  6.15606 |  0.00003 |
|  12.00000 | 110.00000 |  2.29023 |  0.00004 |
|  52.00000 |  90.00000 | 12.67032 |  0.00004 |
|  52.00000 | 100.00000 |  5.85384 |  0.00003 |
|  52.00000 | 110.00000 |  2.05826 |  0.00005 |
| 252.00000 |  90.00000 | 12.61124 |  0.00003 |
| 252.00000 | 100.00000 |  5.78184 |  0.00003 |
| 252.00000 | 110.00000 |  2.00395 |  0.00004 |

Median variance-reduction factor over K in {90, 100, 110}:

| method               |      12 |      52 |     252 |
|:---------------------|--------:|--------:|--------:|
| plain                |     1.0 |     1.0 |     1.0 |
| antithetic           |     2.1 |     2.1 |     2.1 |
| control              |  1307.4 |  1272.7 |  1342.1 |
| antithetic+control   |  1223.7 |  1181.9 |  1191.6 |
| qmc                  |   334.0 |    63.3 |    21.9 |
| qmc (bridge)         |  1526.0 |  1117.3 |  1399.8 |
| qmc+control          | 13172.4 |  3282.8 |  2837.4 |
| qmc+control (bridge) | 55334.5 | 59325.3 | 26653.7 |

Consistency with the reference: mean |z| = 0.90, max |z| = 2.79 over 72 estimates.

## E5 - down-and-out call (K=100, H=90): discretisation and budget allocation

Continuous-monitoring analytic price: 8.66547. Bias of discrete monitoring decays like m^-0.42 (theory m^-0.5).

|   m | corrected   |   price |      se |   bgk_discrete |     bias |
|----:|:------------|--------:|--------:|---------------:|---------:|
|   4 | False       | 9.96123 | 0.00546 |        9.96255 |  1.29575 |
|   4 | True        | 8.92795 | 0.00979 |        9.96255 |  0.26248 |
|   8 | False       | 9.71385 | 0.00684 |        9.71833 |  1.04838 |
|   8 | True        | 8.71446 | 0.01067 |        9.71833 |  0.04899 |
|  16 | False       | 9.48552 | 0.00792 |        9.48775 |  0.82005 |
|  16 | True        | 8.65960 | 0.01093 |        9.48775 | -0.00587 |
|  32 | False       | 9.29183 | 0.00874 |        9.28962 |  0.62636 |
|  32 | True        | 8.67711 | 0.01089 |        9.28962 |  0.01164 |
|  64 | False       | 9.11193 | 0.00936 |        9.12959 |  0.44646 |
|  64 | True        | 8.64598 | 0.01089 |        9.12959 | -0.01949 |
| 128 | False       | 8.99877 | 0.00982 |        9.00556 |  0.33330 |
| 128 | True        | 8.66271 | 0.01090 |        9.00556 | -0.00276 |
| 256 | False       | 8.91831 | 0.01005 |        8.91209 |  0.25284 |
| 256 | True        | 8.67486 | 0.01083 |        8.91209 |  0.00939 |
| 512 | False       | 8.83537 | 0.01040 |        8.84299 |  0.16990 |
| 512 | True        | 8.66240 | 0.01095 |        8.84299 | -0.00307 |

Fixed budget N × m - RMSE against the continuous price:

|   m |   uncorrected |   BGK-corrected |
|----:|--------------:|----------------:|
|   4 |       1.29924 |         0.26540 |
|   8 |       1.04119 |         0.03276 |
|  16 |       0.82175 |         0.02576 |
|  32 |       0.62970 |         0.02581 |
|  64 |       0.47274 |         0.04425 |
| 128 |       0.34641 |         0.06847 |
| 256 |       0.28654 |         0.11848 |

## E6 - time-stepping scheme

Euler weak error decays like m^-0.95 for m >= 4 (theory m^-1); the exact log scheme has no discretisation bias, so one step suffices for European payoffs.

|   m |    euler |   exact |
|----:|---------:|--------:|
|   1 | -0.33119 | 0.00190 |
|   2 |  0.17723 | 0.00190 |
|   4 |  0.16424 | 0.00190 |
|   8 |  0.09778 | 0.00190 |
|  16 |  0.05211 | 0.00190 |
|  32 |  0.02156 | 0.00190 |
|  64 |  0.01297 | 0.00190 |
