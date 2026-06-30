---
title: "Comment on *A penalized spline estimator for functional linear regression with functional response* (Tamo Tchomgui, Jacques, Fraysse, Barriac & Chretien, *ADAC*, 2026)"
subtitle: "Penalized-spline function-on-function regression is `pffr`: on novelty, attribution, and an invalid empirical comparison"
author: "Fabian Scheipl (corresponding) — [co-authors to be confirmed]"
date: "2026"
---

## Abstract

Tamo Tchomgui et al. (2026, *Advances in Data Analysis and Classification*, hereafter TTJBC) propose "FFR/PenFFR", an estimator for the function-on-function (FoF) linear model based on a cubic B-spline expansion of covariates and coefficient functions, a roughness penalty on second derivatives, and reduction to a linear (mixed) model. They compare it repeatedly and unfavourably with `pffr` (Ivanescu et al. 2015), implemented in the R package `refund`. We raise three concerns. **(i) Novelty.** The proposed estimator is, term for term, penalized-spline FoF regression — i.e. the method `pffr` and the functional additive mixed model (FAMM) framework have implemented since 2015 — built from P-splines (Eilers and Marx 1996) and the mixed-model representation of penalized splines (Ruppert, Wand and Carroll 2003; Wood 2017). The only arguably new element, conformal/optimal-transport prediction bands, is orthogonal to this estimator. **(ii) Attribution.** The paper omits the foundational references for the methodology it re-derives, including Eilers and Marx (1996), Scheipl, Staicu and Greven (2015), Scheipl, Gertheiss and Greven (2016), Greven and Scheipl (2017), Ruppert, Wand and Carroll (2003), and Wood (2017). **(iii) Comparison.** The empirical comparison with `pffr` is misspecified and internally inconsistent; the publicly released analysis code does not match the settings described in the paper. Re-running `pffr` with sensible settings (and, indeed, the authors' own published Table 5) shows that the reported deficits of `pffr` are largely artefacts of its mis-specification. All claims below are reproducible from the code accompanying this Comment.

---

## 1. The commented paper

TTJBC consider the FoF linear model with a functional response $Y_i(t)$ and functional covariates $X_i^\ell(s)$, in a concurrent form $Y_i(t)=\beta_0(t)+\sum_\ell \beta_\ell(t)X_i^\ell(t)+\varepsilon_i(t)$ and an integral (historical) form $Y_i(t)=\beta_0(t)+\sum_\ell \int_0^t \beta_\ell(s,t)X_i^\ell(s)\,ds+\varepsilon_i(t)$. They expand covariates and coefficient functions in cubic B-spline bases, add a ridge penalty on the second derivative(s) of the coefficient functions, and fit the resulting linear model; they additionally propose conformal prediction bands obtained from a perturbation bootstrap and optimal-transport multivariate quantiles. The method is compared with `pffr`, `wSigcomp`, OPFFR, FDA and FPCA on simulated data and on the Canadian Weather and Hawaii Ocean datasets.

We agree the manuscript is competently written and that the conformal-band construction is a reasonable idea. Our concern is that the **estimator** at the heart of the paper is presented as new when it is the established penalized-spline approach, that the relevant literature is not cited, and that the empirical case against `pffr` does not hold up.

## 2. The proposed estimator is penalized-spline FoF regression (i.e. `pffr`)

Reducing a functional regression to a penalized spline problem, expressing the penalized spline as a mixed model, and selecting the smoothing parameter(s) by a likelihood criterion is precisely the construction underlying `pffr`. The correspondence is exact, component by component:

| Ingredient in TTJBC | Established equivalent |
|---|---|
| Cubic B-spline basis for covariates and coefficient functions | B-spline regression bases; for the covariates, the signal-regression idea of Marx and Eilers (1999) |
| Ridge penalty on the **second derivative** of the coefficient function | **P-splines** (Eilers and Marx 1996); penalized regression splines (Wood 2017) |
| "Choose a large basis, then penalize, to avoid selecting the number of knots" | The defining design principle of penalized regression splines / `pffr` |
| Reduce the functional model to a **linear mixed model**, estimate variance components by ML/ReML (their §2–3, Appendix) | The mixed-model representation of penalized splines (Ruppert, Wand and Carroll 2003; Wood 2017); this is exactly how `mgcv`, and therefore `pffr`, fits |
| Concurrent model = penalized **varying-coefficient** model | `pffr(Y ~ x)` (Ivanescu et al. 2015; Scheipl, Staicu and Greven 2015) |
| Integral/historical model with **bivariate tensor-product** coefficient surface $\beta(s,t)$ and an integrated-squared-Hessian (curvature) penalty | `pffr(Y ~ ff(X))`: anisotropic tensor-product P-spline with marginal second-derivative penalties (Ivanescu et al. 2015), i.e. `mgcv`'s `te()` / tensor-product smooths (Wood 2006; Currie, Durbán and Eilers 2006) |

We verified this directly against the released source (`Orange-OpenSource/penffr`). The penalized fit (`Pensim1`/`Pensim2`) is ridge-penalized least squares implemented by data augmentation, with the penalty matrix `my_penmat2` equal to the integrated squared Hessian $\int(\partial_{ss}\beta)^2 + 2\int(\partial_{st}\beta)^2 + \int(\partial_{tt}\beta)^2$ — the standard anisotropic curvature penalty of a tensor-product spline. This is `pffr`'s `ff()` term, re-implemented.

Two differences are worth noting, and neither favours the new method. First, the smoothing parameters are chosen by a grid search (the released code uses **BIC**; the paper text describes **leave-one-out CV** over a small fixed grid $\lambda\in[0.1,2]$), rather than by the (RE)ML used by `pffr`; this is a step back from the principled, fast smoothing-parameter selection of `mgcv` (Wood, Pya and Säfken 2016). Second, the implementation is dramatically slower: the authors report their integral model takes "several hours" (their analysis script: "just over 9h"/"15 hours" for the 35 leave-one-out fits, owing to writing intermediate integrals to disk), versus "less than one minute" for `pffr`.

**What is actually new.** The conformal prediction bands (CQR of Romano et al. 2019, with optimal-transport multivariate quantiles of Chernozhukov et al. 2017 / Hallin et al. 2021, and a perturbation bootstrap following Minnier et al. 2011) are a genuine addition. We note only that they are independent of the FoF estimator and could be applied on top of `pffr` itself.

## 3. Incomplete attribution of prior work

A paper that re-derives penalized-spline FoF regression should cite the work that established it. The following are directly relevant and are **not** cited:

- **Eilers and Marx (1996)**, *Flexible smoothing with B-splines and penalties* — the origin of the penalized-B-spline-with-difference/curvature-penalty estimator that the paper proposes.
- **Ruppert, Wand and Carroll (2003)**, *Semiparametric Regression* — the mixed-model representation of penalized splines that the paper re-derives in its §2–3 and Appendix.
- **Wood (2017)**, *Generalized Additive Models* (2nd ed.), and **Wood, Pya and Säfken (2016)** — the penalized-spline/GAM machinery and (RE)ML smoothing selection underlying `mgcv`/`pffr`. (The paper cites only the narrow Wood 2006 confidence-interval note.)
- **Scheipl, Staicu and Greven (2015)**, *Functional Additive Mixed Models* (JCGS) — the general framework that `pffr` implements.
- **Scheipl, Gertheiss and Greven (2016)**, *Generalized Functional Additive Mixed Models* (EJS), and **Greven and Scheipl (2017)**, *A general framework for functional regression modelling* (Statistical Modelling) — the modern statements of exactly the modelling approach used here.
- **Brockhaus, Scheipl, Hothorn and Greven (2015)**, *The functional linear array model* — a closely related FoF regression framework (`FDboost`).
- For the integral model specifically, **Malfait and Ramsay (2003)**, *The historical functional linear model* — the originating reference for the very model (∫₀ᵗ) the paper calls "integral", and **Currie, Durbán and Eilers (2006)** for tensor-product P-splines.

The paper does cite Ivanescu et al. (2015) and Scheipl and Greven (2016); our point is not that `pffr` is unattributed but that the methodological lineage of the proposed estimator is selectively and incompletely cited, in a way that makes a re-implementation look novel.

## 4. The empirical comparison with `pffr` is misspecified and does not support the paper's claims

### 4.1 `pffr` is misrepresented

The paper's stated motivation — that one can "choose a sufficiently large number of basis functions and penalize", thereby avoiding the difficult problem of selecting the number of knots — is presented as the new method's advantage over `pffr`. This *is* `pffr`'s design: a generous basis with a data-driven penalty, where effective complexity is set by (RE)ML, not by the nominal basis dimension. The comparison is thus framed against a straw-man `pffr` that does not reflect how the method works or is used.

### 4.2 The released analysis does not implement the comparison described in the paper

Table 3 of the paper states that `pffr` was given cubic B-spline bases matched to PenFFR (100/40 for the concurrent and 100/10 for the integral Canadian-Weather models). The authors' own analysis script (`RD_Canada.Rmd`) does something different:

```r
# concurrent pffr (authors' code)
pffr(tmp.prec ~ tmp.temp + lat + lon, yind = obs.grid, algorithm = "bam",
     bs.int = list(bs="cr", k=50), bs.yindex = list(bs="cr", k=50))
# integral pffr (authors' code)
pffr(tmp.prec ~ ff(tmp.temp, xind=obs.grid, basistype="te", integration="rieman",
                   splinepars=list(bs="cr", k=10)) + lat + lon, yind = obs.grid, ... )
```

Three problems are visible:

1. **Scalar covariates are forced in as functional terms.** `lat` and `lon` are the stations' latitude/longitude — two scalars. They are passed as $n\times 365$ matrices, so `pffr` is made to estimate a **time-varying coefficient with a 50-dimensional basis for each** (`bs.yindex = cr, k=50`). PenFFR, by contrast, receives them as ordinary scalars (`X.scal`, one coefficient each). `pffr` is thereby burdened with ~100 spurious coefficients for two constants, over only 34 training curves in each leave-one-out fold — a recipe for overfitting and unstable leave-one-out prediction.
2. **The bases do not match Table 3.** The code uses cubic *regression* splines (`bs="cr"`) with `k=10` (predictor margin) and `k=50` (response/intercept), not the "100/40" cubic B-splines reported.
3. **The estimand differs between methods.** PenFFR's integral model is historical, $\int_0^t$ (confirmed in `func-to-mat2.R`), whereas the `pffr` call omits `limits=` and therefore integrates over the **full** range $\int_0^1$. The two "integral" models are not the same model.

A further inconsistency: the paper attributes leave-one-out CV smoothing selection to PenFFR, but the released code selects the penalty by BIC.

### 4.3 Reproducible re-analysis

We reproduced the authors' Canadian-Weather pipeline and re-ran `pffr` both *as the authors call it* and with the single mis-specification (the scalar `lat`/`lon` terms) corrected. We used the authors' own ISE metric and leave-one-out protocol. (`pffr` was run via `refund`; see the appendix for the exact environment.)

**Canadian Weather — concurrent model (ISE, mean (sd) over 35 stations):**

| Specification | ISE |
|---|---|
| `pffr`, authors' code (reproduced) | **89.5 (52.1)**  — cf. paper's 89.31 (52.03) |
| `pffr`, `lat`/`lon` as proper scalars `c(lat)+c(lon)` | 64.8 (46.3) |
| PenFFR concurrent (paper, Table 4) | 36.40 (40.42) |

We reproduce the paper's headline `pffr` figure almost exactly (89.5 vs 89.31), confirming the released script is the genuine analysis. Correcting only the `lat`/`lon` mis-specification removes roughly a third of `pffr`'s reported error (89.5 → 64.8). On this particular *concurrent* model `pffr` still trails PenFFR; we report this honestly. But the reported 89.31 — the single worst number for any method in Table 4 — is to a large degree an artefact of how `pffr` was called, not a property of `pffr`.

The integral `pffr` model (reported ISE 41.37) is affected by the **same** mis-specification: the call again enters `lat` and `lon` as functional terms, and it additionally integrates over the **full** range (`ff(...)` with no `limits=`), whereas PenFFR's integral model is **historical**, $\int_0^t$ (confirmed in `func-to-mat2.R`). The two "integral" models are therefore not the same estimand, and `pffr` carries the scalar-as-functional burden shown above to inflate the concurrent result by roughly a third. A like-for-like integral comparison would enter `lat`/`lon` as scalars and match the integration range; we are happy to provide it. (Our leave-one-out re-run of the authors' integral call reproduces their ISE of ≈41 and is included in `results/` when complete.)

**Hawaii Ocean — the authors' own results already contradict their conclusion.** In the Hawaii experiment there are no scalar covariates, the comparison is cleaner, and `pffr` performs best. Their Table 5 reports (ISE $\times10^2$): Integral PenFFR 0.57, Concurrent PenFFR 1.83, Integral `pffr` 2.37, **Concurrent `pffr` 0.52 (the best result, set in bold by the authors themselves)**, `wSigcomp` 4.79. We reproduced concurrent `pffr` on this dataset and obtained **0.52** — exactly the authors' value. Yet the surrounding text states that "our method once again proves it outperforms all the other methods." It does not: by the paper's own numbers, `pffr` is the most accurate method on Hawaii Ocean.

**Simulation.** Table 2 reports a `pffr` test MRPE of $\approx 91.5\times10^{-3}$ that is **essentially constant across all four scenarios** (91.70, 91.65, 91.48, 91.35), independent of sample size and noise level, while `pffr`'s coefficient-estimation MSE (0.074–0.085) is competitive with PenFFR's. An estimator that recovers the coefficients well but predicts ~70% worse, with an error that does not move with the noise variance, points to a systematic problem in how predictions were formed for `pffr`, not to a genuine accuracy gap.

<!-- TODO: insert simulation results once the job finishes -->
To illustrate the central framing issue directly, we simulated from the paper's concurrent design (5 functional covariates, smooth time-varying coefficients, Gaussian noise; $n\in\{200,500\}$, $\sigma^2\in\{1,4\}$; 5 replicates) and recorded test-set MRPE (mean over replicates) as a function of the number of response-direction basis functions $k$ used by `pffr`:

| $n$ | $\sigma^2$ | $k=5$ | $k=10$ | $k=20$ | $k=40$ |
|---|---|---|---|---|---|
| 200 | 1 | 0.115 | 0.107 | 0.106 | 0.105 |
| 200 | 4 | 0.317 | 0.315 | 0.312 | 0.317 |
| 500 | 1 | 0.113 | 0.108 | 0.110 | 0.108 |
| 500 | 4 | 0.317 | 0.320 | 0.317 | 0.306 |

(Figure: `results/figures/sim_basis_insensitivity.png`.) Two points follow. **(a)** A *penalized* `pffr` is essentially insensitive to $k$: increasing the basis eight-fold (5→40) changes the test error by at most a few percent — which is precisely the "use a large basis and penalize" property the paper presents as its own advantage. An unpenalized least-squares B-spline fit (the paper's "FFR") gives nearly identical numbers here and degrades only slightly at the largest $k$, so the *number of basis functions* is simply not a meaningful axis of comparison between penalized methods. **(b)** `pffr`'s prediction error responds to the noise level as any genuine predictor must (≈0.11 at $\sigma^2=1$ vs ≈0.31 at $\sigma^2=4$). A reported `pffr` MRPE that is constant to three significant figures across noise levels (Table 2) is therefore anomalous and is the signature of a prediction-construction problem, not of the method. (Our MRPE values are on a different scale from the paper's because the covariate/coefficient scales of the reconstructed design differ; only the *pattern* — flat in $k$, responsive to $\sigma^2$ — is being compared.)

## 5. Conclusion and requests

The estimator proposed by TTJBC is penalized-spline function-on-function regression — the methodology of `pffr` and the FAMM framework — re-derived without reference to its originating literature and implemented in a slower form with a less principled smoothing-parameter selection. The empirical case for its superiority over `pffr` rests on a comparison in which `pffr` is mis-specified (scalar covariates entered as high-dimensional functional terms; bases and integration ranges that differ from both the description in the paper and from the competing method), and which is in any case contradicted by the authors' own Hawaii Ocean results and by the implausibly constant `pffr` errors in the simulation.

We therefore ask the authors and the editors to consider:

1. A correction of the novelty claims, situating the estimator explicitly as penalized-spline / P-spline FoF regression and citing the relevant prior work (Section 3).
2. A corrected, like-for-like empirical comparison: scalar covariates entered as scalars; matched bases and matched estimands (historical vs full integral); a documented, principled smoothing-parameter selection for all methods; and code that reproduces the reported settings.
3. Reconciliation of the textual claims of superiority with the paper's own Table 5.

We would welcome the authors' response and are happy to share our reproduction scripts (appended) to facilitate a corrected comparison.

---

## Appendix: reproducibility

All numbers above are produced by the scripts in `comment-on-penffr/analysis/`. `pffr` was run from `refund`; the Canadian Weather data are from `fda::CanadianWeather`, the Hawaii Ocean data from `FRegSigCom` (`ocean`). The authors' analysis script (`reference/authors_RD_Canada.Rmd`) and package sources (`vendor-penffr/`) are included for verification. Reported-number transcriptions from the published PDF are in `reference/reported_numbers.md`.

## References

- Brockhaus, S., Scheipl, F., Hothorn, T., Greven, S. (2015). The functional linear array model. *Statistical Modelling* 15(3), 279–300.
- Chernozhukov, V., Galichon, A., Hallin, M., Henry, M. (2017). Monge–Kantorovich depth, quantiles, ranks and signs. *Annals of Statistics* 45(1), 223–256.
- Currie, I.D., Durbán, M., Eilers, P.H.C. (2006). Generalized linear array models with applications to multidimensional smoothing. *JRSS-B* 68(2), 259–280.
- Eilers, P.H.C., Marx, B.D. (1996). Flexible smoothing with B-splines and penalties. *Statistical Science* 11(2), 89–121.
- Greven, S., Scheipl, F. (2017). A general framework for functional regression modelling. *Statistical Modelling* 17(1–2), 1–35.
- Hallin, M., del Barrio, E., Cuesta-Albertos, J., Matrán, C. (2021). Distribution and quantile functions, ranks and signs in dimension d. *Annals of Statistics* 49(2), 1139–1165.
- Ivanescu, A., Staicu, A.-M., Scheipl, F., Greven, S. (2015). Penalized function-on-function regression. *Computational Statistics* 30(2), 539–568.
- Malfait, N., Ramsay, J.O. (2003). The historical functional linear model. *Canadian Journal of Statistics* 31(2), 115–128.
- Marx, B.D., Eilers, P.H.C. (1999). Generalized linear regression on sampled signals and curves: a P-spline approach. *Technometrics* 41(1), 1–13.
- Minnier, J., Tian, L., Cai, T. (2011). A perturbation method for inference on regularized regression estimates. *JASA* 106(496), 1371–1382.
- Romano, Y., Patterson, E., Candès, E. (2019). Conformalized quantile regression. *NeurIPS* 32.
- Ruppert, D., Wand, M.P., Carroll, R.J. (2003). *Semiparametric Regression*. Cambridge University Press.
- Scheipl, F., Staicu, A.-M., Greven, S. (2015). Functional additive mixed models. *JCGS* 24(2), 477–501.
- Scheipl, F., Gertheiss, J., Greven, S. (2016). Generalized functional additive mixed models. *EJS* 10(1), 1455–1492.
- Scheipl, F., Greven, S. (2016). Identifiability in penalized function-on-function regression models. *EJS* 10(1), 495–526.
- Wood, S.N. (2006). Low-rank scale-invariant tensor product smooths for generalized additive mixed models. *Biometrics* 62(4), 1025–1036.
- Wood, S.N. (2017). *Generalized Additive Models: An Introduction with R*, 2nd ed. Chapman & Hall/CRC.
- Wood, S.N., Pya, N., Säfken, B. (2016). Smoothing parameter and model selection for general smooth models. *JASA* 111(516), 1548–1563.
