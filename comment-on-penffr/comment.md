---
title: "Comment on *A penalized spline estimator for functional linear regression with functional response* (Tamo Tchomgui, Jacques, Fraysse, Barriac & Chretien, *ADAC*, 2026)"
subtitle: "A re-implementation of penalized-spline function-on-function regression, presented as new and supported by an invalid comparison"
author: "Fabian Scheipl (corresponding) — [co-authors to be confirmed]"
date: "2026"
---

## Abstract

Tamo Tchomgui et al. (2026, *Advances in Data Analysis and Classification*; hereafter TTJBC) present "FFR/PenFFR" — spline expansions of covariates and coefficient functions with a curvature penalty, fitted as a linear (mixed) model — as a new method for function-on-function (FoF) regression, and report that it is more accurate and more interpretable than `pffr` (Ivanescu et al. 2015; R package `refund`). Neither claim holds. First, the estimator is penalized-spline FoF regression as implemented in `pffr` and the functional additive mixed model framework since 2015, built from standard components (Eilers and Marx 1996; Ruppert, Wand and Carroll 2003; Wood 2017) whose originating references the paper does not cite; it is in fact a restriction of that framework to homoscedastic Gaussian responses without random effects. Second, the empirical case does not hold up: the authors' code enters two scalar covariates into `pffr` as 50-dimensional functional terms; it implements neither `pffr` nor PenFFR itself with the settings the paper describes — not the stated estimation, smoothing-parameter selection, per-curve random effect, or basis dimensions; and three competitors in the headline table were not run but copied from another paper. Reproducing the authors' analysis, we show their headline `pffr` error is largely an artefact of the mis-specification. Third, the paper's own Table 5 ranks `pffr` the most accurate method on the Hawaii Ocean data, contradicting its conclusion, and the one new ingredient — conformal prediction bands — is neither implemented in the released package nor attains its nominal coverage where evaluated. Each of these defects runs in the authors' favour. We document every point below and ask for a correction of the record: proper attribution of the prior work and a corrected, like-for-like comparison. Complete reproduction materials accompany this Comment.

---

## 1. The commented paper

TTJBC consider the FoF linear model with a functional response $Y_i(t)$ and functional covariates $X_i^\ell(s)$, in a concurrent form $Y_i(t)=\beta_0(t)+\sum_\ell \beta_\ell(t)X_i^\ell(t)+\varepsilon_i(t)$ and an integral/historical form $Y_i(t)=\beta_0(t)+\sum_\ell \int_0^t \beta_\ell(s,t)X_i^\ell(s)\,ds+\varepsilon_i(t)$. They expand covariates and coefficient functions in cubic B-spline bases, penalize the second derivative(s) of the coefficient functions, fit the resulting linear model, and describe conformal prediction bands. The method is compared with `pffr`, `wSigcomp`, OPFFR, FDA and FPCA on simulated data and on the Canadian Weather and Hawaii Ocean datasets, and is claimed to be the most accurate and most interpretable throughout.

We have checked the paper against its released package (`Orange-OpenSource/penffr`) and the authors' own analysis script, and re-run the experiments. This Comment documents that the estimator is not new (§2), that the work it re-derives is not cited (§3), and that the comparison supporting the empirical claims is invalid (§4).

## 2. The proposed estimator is penalized-spline FoF regression — i.e. `pffr`

Expanding covariates and coefficients in a spline basis, penalizing the curvature of the coefficient functions, writing the penalized problem as a mixed model, and selecting the smoothing parameters by a likelihood criterion *is* the construction underlying `pffr`. The correspondence is exact, component by component:

| Ingredient in TTJBC | Established equivalent |
|---|---|
| Cubic B-spline bases for covariates and coefficient functions | B-spline regression; for covariates, signal regression (Marx and Eilers 1999) |
| Ridge penalty on the **second derivative** of the coefficient function | **P-splines** (Eilers and Marx 1996); penalized regression splines (Wood 2017) |
| "Choose a large basis, then penalize, to avoid selecting the number of knots" | The defining design principle of penalized regression splines / `pffr` |
| Reduce the functional model to a **linear mixed model**, estimate variance components by ML/ReML (their §2–3 and Appendix) | The mixed-model representation of penalized splines (Ruppert, Wand and Carroll 2003; Wood 2017) — exactly how `mgcv`, and hence `pffr`, fits |
| Concurrent model = penalized **varying-coefficient** model | `pffr(Y ~ x)` (Ivanescu et al. 2015; Scheipl, Staicu and Greven 2015) |
| Integral/historical model: **bivariate tensor-product** $\beta(s,t)$ with an integrated-squared-Hessian (curvature) penalty | `pffr(Y ~ ff(X))`: anisotropic tensor-product P-spline with marginal second-derivative penalties (Ivanescu et al. 2015); `mgcv` tensor smooths (Wood 2006; Currie, Durbán and Eilers 2006) |

We confirmed this against the source. The penalized fit (`Pensim1`/`Pensim2`) is ridge-penalized least squares by data augmentation; the penalty matrix `my_penmat2` is the integrated squared Hessian $\int(\partial_{ss}\beta)^2 + 2\int(\partial_{st}\beta)^2 + \int(\partial_{tt}\beta)^2$ — the standard anisotropic curvature penalty of a tensor-product spline. This is `pffr`'s `ff()` term, re-implemented.

Where the implementation differs from `pffr`, it differs for the worse. Smoothing parameters are selected by a coarse grid search rather than by the fast, principled (RE)ML of `mgcv` (Wood, Pya and Säfken 2016) — and, as §4.5 documents, not even by the criterion the paper states. The implementation is also far slower: the authors report "several hours" for their integral model (their own script logs 9–15 hours for 35 leave-one-out fits), against "less than one minute" for `pffr`.

The proposal is also strictly less general than the method it re-implements — a step backwards, not forwards. PenFFR fits a single homoscedastic-Gaussian model by least squares (`Pensim1`/`Pensim2` call `lm()`): only a conditional mean is modelled. `pffr` inherits the generality of `mgcv` and fits the same models for non-Gaussian and heteroscedastic responses — Beta, Negative Binomial, Tweedie, scaled-$t$, Gaussian location–scale, zero-inflated and ordinal families, among others — with scalar or functional random effects for longitudinal, hierarchical or otherwise correlated functional data. That is the content of the generalized functional additive mixed model (GFAMM) framework of Scheipl, Gertheiss and Greven (2016), of which PenFFR re-derives the Gaussian special case (§3 gives the correspondence in detail). A method offering a strict subset of an existing method's capabilities, more slowly and without working uncertainty quantification (§4.4), is not a contribution.

## 3. The paper does not attribute the methodology it re-derives

A paper that re-derives penalized-spline FoF regression must cite the work that established it. The following are directly on point and absent from the reference list:

- **Eilers and Marx (1996)**, *Flexible smoothing with B-splines and penalties* — the origin of the penalized-B-spline estimator the paper proposes.
- **Ruppert, Wand and Carroll (2003)**, *Semiparametric Regression* — the mixed-model representation of penalized splines that the paper re-derives in its §2–3 and Appendix.
- **Wood (2017)**, *Generalized Additive Models* (2nd ed.), and **Wood, Pya and Säfken (2016)** — the penalized-spline machinery and (RE)ML smoothing selection underlying `mgcv`/`pffr`. (Only the narrow Wood 2006 confidence-interval note is cited.)
- **Scheipl, Staicu and Greven (2015)**, *Functional Additive Mixed Models* (JCGS) — the framework `pffr` implements. **Scheipl, Gertheiss and Greven (2016)**, *Generalized Functional Additive Mixed Models* (EJS) — contains the framework PenFFR occupies: its equations (1)–(2) and Table 1 give the general response distribution and link, the concurrent effect $v(t)\beta(t)$, the historical/integral effect $\int_{l(t)}^{u(t)} x(s)\beta(s,t)ds$, and functional random effects; its equations (3)–(4) give the penalized tensor-product-spline estimator with the anisotropic (Kronecker-sum) curvature penalty that PenFFR re-derives; its §§4.1–4.2 demonstrate non-Gaussian responses. **Greven and Scheipl (2017)** restate this generality.
- **Malfait and Ramsay (2003)**, *The historical functional linear model* — the originating reference for the model the paper calls "integral" ($\int_0^t$); **Currie, Durbán and Eilers (2006)** for tensor-product P-splines.

That `refund` itself cites the Scheipl/Greven papers as the basis of `pffr` makes the omission harder to understand: the authors cite Ivanescu et al. (2015) and Scheipl and Greven (2016), but not the methodological lineage that would make clear their estimator is not new.

## 4. The empirical comparison is invalid

### 4.1 `pffr` is misrepresented and then mis-specified

The paper motivates its method as being able to "choose a sufficiently large number of basis functions and penalize", avoiding the need to select the number of knots. That is `pffr`'s design, not a point of difference: a generous basis with a data-driven penalty. The comparison is framed against a `pffr` that does not exist.

The actual `pffr` calls are worse than the framing. From the authors' Canadian-Weather script:

```r
# concurrent pffr (authors' code)
pffr(tmp.prec ~ tmp.temp + lat + lon, yind = obs.grid, algorithm = "bam",
     bs.int = list(bs="cr", k=50), bs.yindex = list(bs="cr", k=50))
# integral pffr (authors' code)
pffr(tmp.prec ~ ff(tmp.temp, xind=obs.grid, basistype="te", integration="rieman",
                   splinepars=list(bs="cr", k=10)) + lat + lon, yind = obs.grid, ...)
```

Three defects, each inflating `pffr`'s error:

1. **Two scalar covariates are entered as functional terms.** `lat` and `lon` are the stations' coordinates — one number per station. The code passes them as $n\times 365$ matrices, so `pffr` treats each as a functional covariate and fits a time-varying coefficient with a 50-basis smooth, consuming roughly 26 (lat) and 20 (lon) effective degrees of freedom in our refit: two constants fitted as seasonal curves. PenFFR enters the same coordinates as scalars with one constant coefficient each (`penffr1` appends them as plain design columns). `pffr` fits constant scalar effects via its `c()` notation, which the matched specification below uses. Comparing PenFFR's 2 coefficients with `pffr`'s ~100 is the dominant cause of the reported Canadian-Weather gap, and it is an artefact of the call, not of `pffr`.
2. **The bases contradict the paper's own Table 3.** Table 3 states `pffr` used "100/40" cubic B-splines; the code uses cubic *regression* splines (`bs="cr"`) with `k=10`/`k=50`. The reported settings were not the ones run.
3. **The estimand differs between methods.** PenFFR's integral model is historical, $\int_0^t$ (confirmed in `func-to-mat2.R`); the `pffr` call omits `limits=` and integrates over the full range $\int_0^1$. The two "integral" models are not the same model.

### 4.2 Three of the competitors were not run at all

For the Canadian-Weather comparison (Table 4), the paper states (their §6.1) that "[d]ue to the unavailability of code for the OPFFR approach, we simply use the published results as presented in [Sun et al. 2018]." The analysis script hard-codes the OPFFR, FDA *and* FPCA entries as literal constants — `40.28 (45.76)`, `44.16 (56.95)`, `45.51 (45.78)`: three of the nine rows in the headline table are numbers copied from a different paper, obtained under a different protocol, not a comparison the authors performed. (For Hawaii Ocean the same three methods are dropped "due to the unavailability of the code.") A table mixing the authors' leave-one-out runs with figures lifted from elsewhere is not a controlled benchmark. `wSigcomp`, by contrast, *was* run, with generous settings (80/80/40 basis functions) — the problem is not that every competitor was handicapped, but that the comparison is uncontrolled exactly where the paper draws its conclusions.

### 4.3 Re-analysis: the reported deficits of `pffr` are artefacts

We reproduced the Canadian-Weather pipeline and re-ran `pffr` both as the authors call it and matched to PenFFR as closely as possible: cubic B-spline bases with a second-derivative penalty, PenFFR's coefficient-basis dimensions (40 concurrent, 10×10 integral), scalars entered as scalars, historical integration `limits="s<=t"`, fitted with `mgcv`'s exact-REML `gam` engine. The as-is rows use the authors' own `algorithm="bam"` call; metric and leave-one-out protocol are theirs.

**Canadian Weather (ISE, mean (sd) over 35 stations):**

| Specification | concurrent | integral |
|---|---|---|
| `pffr`, authors' code (reproduced) | **89.5 (52.1)** — cf. paper 89.31 (52.03) | cf. paper 41.37 (48.91) |
| `pffr`, matched to PenFFR (scalars as scalars; historical integral) | 64.80 (46.31) | 37.98 (36.65) |
| PenFFR (paper, Table 4) | 36.40 (40.42) | 33.66 (22.99) |

We reproduce the paper's headline `pffr` figure to within rounding (89.5 vs 89.31), confirming the released script is the analysis behind the paper. Correcting only the scalar-covariate mis-specification removes roughly a third of the reported `pffr` error. A well-specified `pffr` still trails PenFFR on this concurrent model (64.8 vs 36.4) — we do not dispute that PenFFR fits this dataset better. But the reported 89.31, the worst number for any method in Table 4, is substantially an artefact of how `pffr` was called. On the integral model the matched, historical `pffr` (37.98) likewise improves on the reported 41.37, with lower variability (sd 36.7 vs 48.9), while still trailing PenFFR's 33.66. What the paper's numbers overstate is the *size* of the gap; the blanket claim of superiority, the Hawaii data below contradict outright.

**Only the scalar-covariate correction moves these numbers.** With the coordinates entered correctly, the matched `pffr` is insensitive to every remaining modelling choice. Replacing cubic regression splines (`bs="cr"`) with cubic B-splines (`bs="ps"`, second-derivative penalty) changes the concurrent ISE by less than one unit at the authors' settings (89.9 vs 89.5) and moves the matched historical integral from 37.98 to 39.58. GCV in place of exact REML gives 64.9 / 38.7 (concurrent / integral) against 64.8 / 38.0. `mgcv`'s exact-REML `gam` and fast-REML `bam` engines agree to the precision reported. Pre-smoothing each temperature curve into a 100-dimensional B-spline representation before modelling, as PenFFR does, moves the concurrent ISE from 64.80 to 65.02. The one choice that matters is whether the two coordinates enter as scalars or as 50-basis functional terms (89.5 → 64.8). The gap that remains to PenFFR is a genuine, if modest, difference on this dataset — not an artefact of how `pffr` was tuned, and not grounds for the paper's conclusions.

**Hawaii Ocean: `pffr` is the best method in the authors' own table.** Here there are no scalar covariates and the comparison is cleaner. Table 5 reports (ISE $\times10^2$): Integral PenFFR 0.57, Concurrent PenFFR 1.83, Integral `pffr` 2.37, Concurrent `pffr` **0.52** — the best result, in the authors' own bold — and `wSigcomp` 4.79. We reproduced concurrent `pffr` and obtained 0.52, the reported value to two figures. Yet the text states that "our method once again proves it outperforms all the other methods." By the paper's own numbers, it does not.

**Simulation.** Table 2 reports a `pffr` MRPE of about $91.5\times10^{-3}$, constant to three figures across all four scenarios (91.70, 91.65, 91.48, 91.35) — unmoved by sample size and by the noise variance — while `pffr`'s coefficient-estimation MSE (0.074–0.085) is competitive with PenFFR's. An estimator that recovers the coefficients well but predicts ~70% worse, with a prediction error that does not respond to the noise level, points to a defect in how the predictions were formed, not to an accuracy gap. Simulating from the paper's concurrent design with `pffr` correctly specified, its MRPE is flat in the basis dimension (0.105–0.115 across $k=5,\dots,40$ at $\sigma^2=1$) — the insensitivity the paper claims as its distinguishing advantage — and scales with the noise level (about 0.11 at $\sigma^2=1$, 0.31 at $\sigma^2=4$), as any working predictor must (figure `results/figures/sim_basis_insensitivity.png`).

### 4.4 The one new element — uncertainty quantification — is neither delivered nor validated

The paper's one ingredient with a claim to novelty is its functional prediction band: a perturbation bootstrap (cf. Minnier et al. 2011) combined with optimal-transport multivariate quantiles (Chernozhukov et al. 2017; Hallin et al. 2021) and conformal calibration (Romano et al. 2019). Three problems:

- **It is not implemented in the released software.** In `Orange-OpenSource/penffr`, the helpers for the procedure (`generate_Grid`, `mvt.quant`, `OR` in `utils.R`) are never called by any model, prediction or band-construction function; the code contains no resampling loop; and every prediction function (`pred.PenFFR`, `pred.penffr1/2`) returns point predictions only. As released, PenFFR provides no uncertainty quantification, and the procedure of the paper's Section 4 cannot be reproduced from its code.
- **Where it is evaluated, it fails.** Coverage (CovP) is reported for FFR/PenFFR only, never for a competitor. For a nominal 95% band the reported coverage is 2–13% (Table 2) and 2% and 54% for the two illustrative curves (Fig. 5). A conformal procedure whose point is a coverage guarantee has not been shown to deliver one.
- **It is redundant.** `pffr`/`mgcv` already provide calibrated uncertainty quantification: standard errors and pointwise and simultaneous credible bands for all terms and predictions.

### 4.5 The released code does not implement the method the paper describes

The defects in §4.1 concern how `pffr` was configured. They are compounded by discrepancies between how PenFFR itself is described and what the released package computes: the numbers reported for PenFFR were not produced by the procedure the paper sets out.

- **Estimation.** The paper reduces the functional model to a linear mixed model estimated by ReML (their Appendix). The code does no such thing: `Pensim1`/`Pensim2` compute ridge-penalized least squares by data augmentation — appending $\sqrt{\lambda}\,\text{chol}(\text{penalty})$ rows to the design matrix and calling `lm()`. No mixed model is fitted and no variance component is estimated.
- **Smoothing-parameter selection.** The main text selects $\lambda$ by leave-one-out cross-validation over ten equally spaced values in $[0.1, 2.0]$ (their §3.1 and §6); the Appendix instead describes ReML. The code uses neither: $\lambda$ is chosen by BIC over `seq(0.01, 5)` (concurrent) and `seq(0, 5)` (integral), with `n.lam = min(10, floor(100^(1/(d+1))))` grid points per coefficient function for $d$ functional covariates. For the four-covariate Hawaii Ocean model this evaluates to `n.lam = 2`: the "grid" consists of its two endpoints, no interior value of $\lambda$ is ever tried, and the amount of regularization is effectively not tuned at all. Neither the criterion, nor the range, nor the number of grid points matches the paper, and a two-point endpoint search is not a credible data-driven choice of a smoothing parameter.
- **The per-curve random intercept is not fitted.** The paper's model includes a per-curve random intercept to capture dependence along the response curve (their eqs. 6–7) — the mechanism it offers for within-curve autocorrelation. In the code, `penffr1()`/`penffr2()` fit a fixed-effects `lm()`: the curve identifier is dropped from the model matrix and no random effect or variance component is estimated. The helper that would build the random-effect design (`Rdeff` in `utils.R`) is defined but never called. (A random intercept appears only in the separate mixture-of-experts routines, via `FLXMRlmm(random = ~1)`, which are not the models reported in the tables.)

Together with §4.1 — cubic regression splines rather than the cubic B-splines named in Table 3, scalar coordinates entered as 50-basis functional terms, full-range rather than historical integration — this means **neither the proposed method nor its principal competitor was run as the paper describes.** Tables 2, 4 and 5 were generated by code that departs from the paper's account of it in the estimation method, the smoothing-parameter selection, the random-effect structure, the basis type and dimensions, and the integration range. That is an independent reason the comparison cannot be relied upon.

## 5. A pattern, and what we ask

These defects are not independent accidents that happen to cancel out. Each runs in the same direction, against `pffr` and for PenFFR: scalar covariates inflated into 50-dimensional functional terms for `pffr` only; bases and an integration range matching neither the paper's Table 3 nor PenFFR's own model; three competitors copied from another paper rather than run; a simulated `pffr` prediction error frozen at a constant while its estimation accuracy is fine; coverage reported only for the proposed method. A comparison in which every uncontrolled choice favours the authors' method is not evidence that the method is superior. Nothing here requires a finding of intent — this is what results when authors optimise their own method and run competitors with ill-chosen settings — and the effect on the conclusions is the same either way.

In sum: the paper's central methodological claim, a new estimator, is false — it is penalized-spline FoF regression, uncited, and a strict restriction of it (Gaussian-only, no random effects, no working uncertainty quantification). Its central empirical claim, superiority over `pffr`, is unsupported by its own results and contradicted by its own Table 5, and rests on code that implements neither method as described. Its one novel ingredient is absent from the released software and misses its stated coverage guarantee where evaluated at all. These are not matters of wording. They call for a substantial correction: proper attribution of the prior work, and a corrected, like-for-like comparison or withdrawal of the superiority claims. We are not asking for retraction; we are asking that the record be set right.

Every point above is elementary to verify — the released code can be compared line by line with the paper's Table 3 and Appendix, the missing references are canonical, and the headline conclusion is contradicted by the paper's own Table 5 — so a correction should be straightforward to adjudicate. We raise these points to correct the literature, not to litigate priority; we recognise the work invested in the paper, we provide complete reproduction materials to the authors, reviewers and editors, and we support the authors' right of reply.

---

## Appendix: reproducibility

All numbers are produced by the scripts in `comment-on-penffr/analysis/`. `pffr` is run from `refund`; Canadian Weather is `fda::CanadianWeather`, Hawaii Ocean is `FRegSigCom`'s `ocean`. The authors' analysis script (`reference/authors_RD_Canada.Rmd`) and package sources (`vendor-penffr/`) are included for verification; the transcribed reported numbers are in `reference/reported_numbers.md`.

## References

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
- Sun, X., Du, P., Wang, X., Ma, P. (2018). Optimal penalized function-on-function regression under a reproducing kernel Hilbert space framework. *JASA* 113(524), 1601–1611.
- Wood, S.N. (2006). Low-rank scale-invariant tensor product smooths for generalized additive mixed models. *Biometrics* 62(4), 1025–1036.
- Wood, S.N. (2017). *Generalized Additive Models: An Introduction with R*, 2nd ed. Chapman & Hall/CRC.
- Wood, S.N., Pya, N., Säfken, B. (2016). Smoothing parameter and model selection for general smooth models. *JASA* 111(516), 1548–1563.
