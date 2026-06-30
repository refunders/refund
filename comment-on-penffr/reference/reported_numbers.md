# Numbers as reported in Tamo Tchomgui et al. (2026), ADAC

Read directly from the published PDF (tables are images; values transcribed from rendered pages).

## Table 2 — Simulation (concurrent model), average (sd)
Columns: Test-set MRPE (×10³) · R̃² · CovP · MSE (of estimated parameters)

| Scenario | Method | MRPE×10³ | R̃² | CovP | MSE |
|---|---|---|---|---|---|
| S1 (n=200, σ²=4) | FFR    | 55.86 (0.25) | 0.922 (0.0016) | 11.8% (0.036) | 0.109 (0.006) |
|                  | PenFFR | 55.48 (0.33) | 0.906 (0.0007) | 10.6% (0.058) | 0.088 (0.003) |
|                  | pffr   | 91.70 (0.07) | 0.882 (0.0003) | –             | 0.085 (0.002) |
| S2 (n=200, σ²=1) | FFR    | 54.59 (0.20) | 0.925 (0.0004) | 13.4% (0.058) | 0.079 (0.007) |
|                  | PenFFR | 53.78 (0.17) | 0.912 (0.0006) | 11.7% (0.101) | 0.073 (0.009) |
|                  | pffr   | 91.65 (0.06) | 0.883 (0.0001) | –             | 0.074 (0.001) |
| S3 (n=500, σ²=4) | FFR    | 54.52 (0.13) | 0.924 (0.0008) | 2.6% (0.005)  | 0.077 (0.004) |
|                  | PenFFR | 54.45 (0.07) | 0.909 (0.0004) | 4.9% (0.029)  | 0.070 (0.006) |
|                  | pffr   | 91.48 (0.05) | 0.883 (0.0001) | –             | 0.074 (0.001) |
| S4 (n=500, σ²=1) | FFR    | 54.00 (0.07) | 0.926 (0.0003) | 9.0% (0.013)  | 0.058 (0.003) |
|                  | PenFFR | 53.88 (0.09) | 0.921 (0.0003) | 8.6% (0.015)  | 0.054 (0.005) |
|                  | pffr   | 91.35 (0.02) | 0.883 (0.0001) | –             | 0.074 (0.001) |

Observations:
- pffr MRPE is essentially **constant (~91.5)** across all 4 scenarios — independent of n and noise level, which is implausible for a real predictor and points to a systematic prediction-setup artifact (e.g. intercept/scale handling), not a genuine accuracy gap.
- pffr's coefficient-estimation **MSE (0.074–0.085) is competitive** with PenFFR's (0.054–0.088) and better than FFR's in S1 — i.e. pffr estimates β about as well, but its *predictions* are reported ~70% worse. Estimation OK + prediction bad ⇒ the prediction step is the problem.
- CovP is "–" for pffr (their conformal CI is not applied to pffr) — fine, not a pffr deficiency.

## Table 3 — Basis settings AS REPORTED (cf. actual code in authors_RD_Canada.Rmd)
Columns: Type of basis · X_i^ℓ(t) (covariate) · β_ℓ(t) (coefficient)

| Method | CW basis | CW X | CW β | HO basis | HO X | HO β |
|---|---|---|---|---|---|---|
| Integral PenFFR/FFR | cubic B-splines | 100 | 10 | cubic B-splines | 40 | 6 |
| Concurrent PenFFR/FFR | cubic B-splines | 100 | 40 | cubic B-splines | 40 | 20 |
| Integral pffr | cubic B-splines | 100 | 10 | cubic B-splines | 40 | 6 |
| Concurrent pffr | cubic B-splines | 100 | 40 | cubic B-splines | 40 | 20 |
| wSigcomp | wavelets+SVD | 40 | 80 | wavelets+SVD | 20 | 40 |

**Mismatch with the released code (authors_RD_Canada.Rmd):** the integral pffr call uses
`ff(tmp.temp, basistype="te", splinepars=list(bs="cr", k=10))` + `bs.yindex=list(bs="cr",k=50)` + `bs.int=list(bs="cr",k=50)`, and the concurrent pffr uses `bs.yindex=list(bs="cr",k=50)`.
That is cubic **regression** splines (cr), not B-splines, with k=10 (predictor margin) / k=50 (response, intercept) — not the "100/40" cubic B-splines claimed in Table 3. The pffr fits also include `lat` and `lon` as **functional (time-varying, k=50) terms** by passing them as n×365 matrices, whereas PenFFR receives them as scalars (one coefficient each).

## Table 4 — Canadian Weather, ISE avg (sd); best in bold
| Method | ISE |
|---|---|
| **Integral PenFFR** | **33.66 (22.99)** |
| Concurrent PenFFR | 36.40 (40.42) |
| Integral FFR | 34.63 (26.03) |
| Concurrent FFR | 36.50 (40.51) |
| Integral pffr | 41.37 (48.91) |
| **Concurrent pffr** | **89.31 (52.03)**  ← >2× worse than all others |
| wSigcomp | 45.37 (52.45) |
| OPFFR | 40.28 (45.76) |
| FDA | 44.16 (56.95) |
| FPCA | 45.51 (45.78) |

## Table 5 — Hawaii Ocean, ISE (×10²) avg (sd); best in bold
| Method | ISE×10² |
|---|---|
| Integral PenFFR | 0.57 (0.74) |
| Concurrent PenFFR | 1.83 (0.88) |
| Integral pffr | 2.37 (1.55) |
| **Concurrent pffr** | **0.52 (0.26)**  ← BEST in their own table |
| wSigcomp | 4.79 (4.46) |

**Contradiction:** the text (Sect. 6.2) states "our method once again proves it outperforms all the other methods", but their own Table 5 bolds **Concurrent pffr (0.52)** as the best, ahead of Integral PenFFR (0.57) and far ahead of Concurrent PenFFR (1.83).
