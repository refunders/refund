---
title: "Comment on *A penalized spline estimator for functional linear regression with functional response* (Tamo Tchomgui, Jacques, Fraysse, Barriac & Chretien, *ADAC*, 2026)"
subtitle: "A re-implementation of penalized-spline function-on-function regression, presented as new and supported by an invalid comparison"
author: "Fabian Scheipl (corresponding) — [co-authors to be confirmed]"
date: "2026"
---

## Abstract

We have read Tamo Tchomgui et al. (2026, *Advances in Data Analysis and Classification*; hereafter TTJBC) with concern. The paper presents "FFR/PenFFR" — a cubic B-spline expansion of covariates and coefficient functions, a roughness penalty on second derivatives, reduced to a linear (mixed) model — as a new method for function-on-function (FoF) regression, and claims it is more accurate and more interpretable than `pffr` (Ivanescu et al. 2015, in the R package `refund`). Both claims are unfounded. **First**, the estimator is not new: it is penalized-spline (P-spline) FoF regression, the method `pffr` and the functional additive mixed model (FAMM) framework have implemented since 2015, assembled from textbook components (Eilers and Marx 1996; Ruppert, Wand and Carroll 2003; Wood 2017) whose originating references the paper does not cite. It is, moreover, a strict *restriction* of that method: PenFFR fits only a homoscedastic-Gaussian, least-squares model with no mixed effects and (see below) no working uncertainty quantification, whereas `pffr` handles general non-Gaussian and heteroscedastic responses and functional mixed models — so the proposal is a step backwards, not forwards. **Second**, the empirical case against `pffr` does not survive inspection. We reproduce the authors' own headline result exactly and show that it is produced by mis-specifying `pffr` — entering two scalar covariates as 50-dimensional functional terms; the released code does not even implement the basis settings the paper reports in its Table 3; three of the competitors in the main table were not run at all but copied from another paper; and the one element that *is* new — conformal prediction bands — is neither implemented in the released software nor shown to achieve its claimed coverage. **Third**, the paper's own Table 5 already ranks `pffr` as the best method, contradicting its conclusion. Every one of these defects runs in the same direction — against `pffr` and in favour of PenFFR. The contribution and the comparison that the paper rests on do not hold; we set out the evidence below and ask the editors for a correction or retraction. All claims are reproducible from the materials accompanying this Comment.

---

## 1. The commented paper

TTJBC consider the FoF linear model with a functional response $Y_i(t)$ and functional covariates $X_i^\ell(s)$, in a concurrent form $Y_i(t)=\beta_0(t)+\sum_\ell \beta_\ell(t)X_i^\ell(t)+\varepsilon_i(t)$ and an integral/historical form $Y_i(t)=\beta_0(t)+\sum_\ell \int_0^t \beta_\ell(s,t)X_i^\ell(s)\,ds+\varepsilon_i(t)$. They expand covariates and coefficient functions in cubic B-spline bases, add a ridge penalty on the second derivative(s) of the coefficient functions, fit the resulting linear model, and additionally describe conformal prediction bands. The method is compared with `pffr`, `wSigcomp`, OPFFR, FDA and FPCA on simulated data and on the Canadian Weather and Hawaii Ocean datasets, and is claimed to be the most accurate and most interpretable throughout.

We have verified the paper against the authors' released package (`Orange-OpenSource/penffr`) and their own Canadian-Weather analysis script, and we have re-run the relevant experiments. Our conclusion is that the paper makes no methodological contribution beyond what `pffr` already provides, fails to attribute the methodology it uses, and supports its empirical claims with a comparison that is mis-specified, partly not executed at all, and contradicted by the paper's own results.

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

Where the implementation differs from `pffr`, it is worse, not better. The smoothing parameters are chosen by a coarse grid search — and not even by the criterion the paper claims: the main text describes leave-one-out CV over $\lambda\in[0.1,2]$, the Appendix describes ReML, and the released code uses **BIC** — rather than by the principled, fast (RE)ML of `mgcv` (Wood, Pya and Säfken 2016). And the implementation is far slower: the authors report their integral model takes "several hours" (their own script: "just over 9h" / "15 hours" for 35 leave-one-out fits, because it writes intermediate integrals to disk), against "less than one minute" for `pffr`.

**The proposal is also strictly less general than the method it re-implements — it is a step backwards.** PenFFR fits a single homoscedastic-Gaussian, squared-error model: estimation is ordinary least squares (`Pensim1`/`Pensim2` call `lm()`), so the response is assumed Gaussian with constant variance and only a conditional mean is modelled. `pffr` inherits the full generality of `mgcv`: via its `family` argument it fits function-on-function regression for any exponential-family response and for non-exponential-family distributions — including heteroscedastic Gaussian location–scale (`gaulss`), Beta (`betar`), scaled-$t$ (`scat`), Negative Binomial, Tweedie, ordered-categorical and zero-inflated models — and it fits these as functional additive *mixed* models, with scalar or functional random effects for longitudinal, hierarchical, crossed, or spatially/temporally correlated functional data. This is not incidental: it is exactly the content of the generalized functional additive mixed model framework of Scheipl, Gertheiss and Greven (2016), whose model — their equations (1)–(2) and Table 1 (§2.1) — already specifies the response distribution $F(\mu_i(t),\nu)$ and link $g$, the concurrent effect $v(t)\beta(t)$, the historical/integral effect $\int_{l(t)}^{u(t)} x(s)\beta(s,t)\,ds$, and functional random effects $b_g(t)$, all estimated by the penalized tensor-product-spline construction with anisotropic (Kronecker-sum) curvature penalty in their equations (3)–(4) (§2.3), with non-Gaussian responses demonstrated in their §§4.1–4.2. PenFFR re-derives the Gaussian special case of this construction, drops the generalized responses, drops the mixed-effects/longitudinal capability, and (as we show in §4.4) does not deliver the uncertainty quantification it claims. A method that offers a strict subset of an existing method's capabilities, less accurately and far more slowly, is not a contribution.

## 3. The paper does not attribute the methodology it re-derives

A paper that re-derives penalized-spline FoF regression must cite the work that established it. The following are directly on point and are **absent** from the reference list:

- **Eilers and Marx (1996)**, *Flexible smoothing with B-splines and penalties* — the origin of the penalized-B-spline-with-curvature-penalty estimator the paper proposes.
- **Ruppert, Wand and Carroll (2003)**, *Semiparametric Regression* — the mixed-model representation of penalized splines that the paper re-derives in §2–3 and its Appendix.
- **Wood (2017)**, *Generalized Additive Models* (2nd ed.), and **Wood, Pya and Säfken (2016)** — the penalized-spline/GAM machinery and (RE)ML smoothing selection underlying `mgcv`/`pffr`. (Only the narrow Wood 2006 confidence-interval note is cited.)
- **Scheipl, Staicu and Greven (2015)**, *Functional Additive Mixed Models* (JCGS) — the framework `pffr` implements. **Scheipl, Gertheiss and Greven (2016)**, *Generalized Functional Additive Mixed Models* (EJS) — establishes precisely the framework PenFFR occupies: equations (1)–(2) and Table 1 (§2.1) give the general response distribution and link, the concurrent effect $v(t)\beta(t)$, the historical/integral effect $\int_{l(t)}^{u(t)} x(s)\beta(s,t)ds$, and functional random effects; equations (3)–(4) (§2.3) give the penalized tensor-product-spline estimator with the anisotropic (Kronecker-sum) curvature penalty that PenFFR re-derives; §§4.1–4.2 demonstrate non-Gaussian responses. **Greven and Scheipl (2017)**, *A general framework for functional regression modelling* (Statistical Modelling) — restates this generality.
- **Malfait and Ramsay (2003)**, *The historical functional linear model* — the originating reference for the very model the paper calls "integral" ($\int_0^t$); **Currie, Durbán and Eilers (2006)** for tensor-product P-splines.

That `refund` itself cites the Scheipl/Greven papers as the basis of `pffr` makes the omission harder to understand: the authors cite Ivanescu et al. (2015) and Scheipl and Greven (2016) but not the methodological lineage that would make clear their estimator is not new.

## 4. The empirical comparison is invalid

### 4.1 `pffr` is misrepresented and then mis-specified

The paper motivates its method as being able to "choose a sufficiently large number of basis functions and penalize", avoiding the need to select the number of knots. That is `pffr`'s design, not a point of difference: a generous basis with a data-driven penalty, where complexity is set by (RE)ML. The comparison is framed against a `pffr` that does not exist.

The actual `pffr` calls are worse than the framing. From the authors' Canadian-Weather script:

```r
# concurrent pffr (authors' code)
pffr(tmp.prec ~ tmp.temp + lat + lon, yind = obs.grid, algorithm = "bam",
     bs.int = list(bs="cr", k=50), bs.yindex = list(bs="cr", k=50))
# integral pffr (authors' code)
pffr(tmp.prec ~ ff(tmp.temp, xind=obs.grid, basistype="te", integration="rieman",
                   splinepars=list(bs="cr", k=10)) + lat + lon, yind = obs.grid, ...)
```

Three defects, all of which inflate `pffr`'s error:

1. **Two scalar covariates are entered as functional terms.** `lat` and `lon` are the stations' coordinates — one number per station. They are passed as $n\times 365$ matrices, so `pffr` treats each as a functional covariate and fits a **time-varying coefficient** $\beta(t)$ with a 50-basis smooth (`bs.yindex = cr, k=50`). In our refit this consumes 50 coefficients per covariate and **≈26 (lat) and ≈20 (lon) effective degrees of freedom** — spent fitting two constants as seasonal curves, over 34 training stations per fold. PenFFR, by contrast, enters the coordinates as **scalars with a single, constant coefficient each** (`penffr1` appends them as plain design columns; see `func-to-mat1.R`), and `pffr` fits *time-varying* effects of scalar covariates only by default — the constant effect PenFFR uses is obtained with the `c()`-notation. Accordingly, the matched `pffr` we report below uses `c(lat)+c(lon)`, giving one constant coefficient each (lat $=-0.004$, lon $=-0.005$), exactly as in PenFFR. This scalar-as-function inflation — comparing PenFFR's 2 coefficients against `pffr`'s ~100 — is the dominant cause of `pffr`'s poor Canadian-Weather result, and it is entirely an artefact of the call, not of `pffr`.
2. **The bases contradict the paper's own Table 3.** Table 3 states `pffr` used "100/40" cubic B-splines; the code uses cubic *regression* splines (`bs="cr"`) with `k=10`/`k=50`. The reported experimental settings were not the ones run.
3. **The estimand differs between methods.** PenFFR's integral model is historical, $\int_0^t$ (confirmed in `func-to-mat2.R`); the `pffr` call omits `limits=` and integrates over the full range $\int_0^1$. The two "integral" models are not the same model.

### 4.2 Three of the competitors were not run at all

For the Canadian-Weather comparison (Table 4), the paper states (§6.1) that, "[d]ue to the unavailability of code for the OPFFR approach, we simply use the published results as presented in [Sun et al. 2018]." The authors' analysis script hard-codes the OPFFR, FDA *and* FPCA entries as literal constants — `40.28 (45.76)`, `44.16 (56.95)`, `45.51 (45.78)` — i.e. three of the nine rows in the headline table are numbers copied from a different paper, on a different experimental protocol, not a comparison the authors performed. (For Hawaii Ocean these three methods are simply dropped "due to the unavailability of the code.") A table that mixes the authors' own leave-one-out runs with figures lifted from elsewhere is not a controlled benchmark. By contrast, `wSigcomp` *was* run, and with generous settings (80/80/40 basis functions) — so the problem is not that every competitor was handicapped, but that the comparison is uncontrolled exactly where the paper draws its conclusions.

### 4.3 Re-analysis: the reported deficits of `pffr` are artefacts

We reproduced the Canadian-Weather pipeline and re-ran `pffr` both as the authors call it and with a specification matched to PenFFR as closely as possible (cubic B-spline `ps` bases with a second-derivative penalty; the same coefficient-basis dimensions as PenFFR — 40 concurrent, 10×10 integral; the scalars entered as scalars; historical integration `limits="s<=t"` to match PenFFR's estimand; fit with `mgcv`'s exact-REML `gam` engine). The authors' *as-is* rows use their own `algorithm="bam"` call. The metric and leave-one-out protocol are the authors' own.

**Canadian Weather (ISE, mean (sd) over 35 stations):**

| Specification | concurrent | integral |
|---|---|---|
| `pffr`, authors' code (reproduced) | **89.5 (52.1)** — cf. paper 89.31 (52.03) | cf. paper 41.37 (48.91) |
| `pffr`, matched to PenFFR (scalars as scalars; historical integral) | 64.80 (46.31) | 37.98 (36.65) |
| PenFFR (paper, Table 4) | 36.40 (40.42) | 33.66 (22.99) |

We reproduce the paper's headline `pffr` figure essentially exactly (89.5 vs 89.31), confirming the released script is the genuine analysis. Correcting only the scalar-covariate mis-specification removes roughly a third of the reported `pffr` error. We report honestly that, on this *concurrent* model, a well-specified `pffr` still trails PenFFR (64.8 vs 36.4); but the reported 89.31 — the worst number for any method in Table 4 — is substantially an artefact of how `pffr` was called, not a property of `pffr`. On the integral model the matched, historical `pffr` (37.98) likewise improves on the reported 41.37 and lowers its run-to-run variability (sd 36.7 vs 48.9), while still trailing PenFFR's 33.66. In other words: on this particular dataset PenFFR does fit somewhat better than a correctly-specified `pffr`; what is illegitimate is the *size* of the reported gap, much of which is mis-specification, and the blanket claim of superiority — which the Hawaii data below flatly contradict.

**Hawaii Ocean — `pffr` is the best method in the authors' own table.** Here there are no scalar covariates, the comparison is cleaner, and Table 5 reports (ISE $\times10^2$): Integral PenFFR 0.57, Concurrent PenFFR 1.83, Integral `pffr` 2.37, **Concurrent `pffr` 0.52 — the best result, in the authors' own bold** — `wSigcomp` 4.79. We reproduced concurrent `pffr` and obtained **0.52**, the authors' value to two figures. Yet the text states that "our method once again proves it outperforms all the other methods." By the paper's own numbers, it does not.

**Simulation.** Table 2 reports a `pffr` MRPE of $\approx 91.5\times10^{-3}$ that is **constant to three figures across all four scenarios** (91.70, 91.65, 91.48, 91.35), independent of sample size and noise, while `pffr`'s coefficient-estimation MSE (0.074–0.085) is competitive with PenFFR's. An estimator that recovers coefficients well but predicts ~70% worse, with an error that does not move with the noise variance, indicates a defect in how `pffr`'s predictions were formed — not a genuine accuracy gap. Simulating from the paper's concurrent design and using `pffr` correctly, its MRPE (a) is essentially flat in the number of basis functions $k$ (0.105–0.115 across $k=5\ldots40$ at $\sigma^2=1$) — the very property the paper claims as its advantage — and (b) responds to the noise level (≈0.11 at $\sigma^2=1$ vs ≈0.31 at $\sigma^2=4$), as any real predictor must. (Figure `results/figures/sim_basis_insensitivity.png`.)

### 4.4 The one new element — uncertainty quantification — is neither delivered nor validated

The paper's only ingredient with a claim to novelty is the functional prediction band ("a perturbation bootstrap" with optimal-transport multivariate quantiles, conformalized). Three problems:

- **It is not implemented in the released software.** In `Orange-OpenSource/penffr`, the helpers for this procedure (`generate_Grid`, `mvt.quant`, `OR` in `utils.R`) are defined but **never called** by any model, prediction, or band-construction function; there is no perturbation bootstrap in the code (no resampling loop), and every prediction function (`pred.PenFFR`, `pred.penffr1/2`) returns **point predictions only**. As released, PenFFR provides no uncertainty quantification at all; the procedure of §4 cannot be reproduced from the published code.
- **Where it is evaluated, it fails.** Coverage (`CovP`) is reported only for FFR/PenFFR (never for `pffr` or the other competitors). For a nominal **95%** band the reported coverage is ~2–13% (Table 2) and 2% / 54% for the two illustrative curves (Fig. 5). A conformal procedure whose entire selling point is a coverage guarantee, attaining 2–54% at nominal 95%, has not been shown to work.
- **It is redundant.** `pffr`/`mgcv` already provide principled, well-calibrated uncertainty quantification (approximate Bayesian credible bands, standard errors, simultaneous intervals) out of the box — so even on UQ the paper offers nothing `pffr` lacks.

## 5. A pattern, and what we ask

The defects above are not independent accidents that happen to cancel out. Every one of them runs in the same direction — making `pffr` (and the other competitors) look worse and PenFFR look better: scalar covariates inflated into 50-dimensional functional terms for `pffr` only; `pffr` bases and an integration range that match neither the paper's Table 3 nor PenFFR's own model; three competitors not run but copied from another paper; a simulated `pffr` error frozen at a constant while its estimation accuracy is fine; coverage reported only for the proposed method. A comparison in which every uncontrolled choice favours the authors' method cannot be presented as evidence that the method is superior, and it was not a competent, good-faith evaluation of the alternatives.

Taken together: the paper's central methodological claim (a new estimator) is false — it is penalized-spline FoF regression, uncited, and in fact a strict *restriction* of it (homoscedastic-Gaussian only, no mixed-effects/longitudinal capability, no working uncertainty quantification), i.e. a step backwards from `pffr` rather than an advance; its central empirical claim (superiority over `pffr`) is unsupported by, and in the Hawaii case contradicted by, its own results; and its one novel ingredient is absent from the released software and fails its stated guarantee where it is evaluated at all. These are not matters that a revision of wording can repair. We therefore ask the authors and the editors for, at minimum, a **substantial correction**, and we believe the case for **retraction** is strong.

We also note, with respect, that these are elementary and checkable errors: the released code contradicts the paper's own Table 3; the foundational references are standard and absent; the sole claimed novelty is not implemented in the accompanying package; and the headline conclusion is contradicted by the paper's own Table 5. That none of this was caught is a failure of the review and editorial process for this article, and we ask the editors to take it into account in deciding the appropriate remedy.

We have no wish to litigate priority for its own sake, and we recognise the authors' effort; our concern is the integrity of the published record. We are glad to provide our complete reproduction materials to the authors and reviewers, and we support the authors' right of reply.

---

## Appendix: reproducibility

All numbers are produced by the scripts in `comment-on-penffr/analysis/`. `pffr` is run from `refund`; Canadian Weather is `fda::CanadianWeather`, Hawaii Ocean is `FRegSigCom`'s `ocean`. The authors' analysis script (`reference/authors_RD_Canada.Rmd`) and package sources (`vendor-penffr/`) are included for verification, and the transcribed reported numbers are in `reference/reported_numbers.md`.

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
- Sun, X., Du, P., Wang, X., Ma, P. (2018). Optimal penalized function-on-function regression under a reproducing kernel Hilbert space framework. *JASA* 113(524), 1601–1611.
- Wood, S.N. (2006). Low-rank scale-invariant tensor product smooths for generalized additive mixed models. *Biometrics* 62(4), 1025–1036.
- Wood, S.N. (2017). *Generalized Additive Models: An Introduction with R*, 2nd ed. Chapman & Hall/CRC.
- Wood, S.N., Pya, N., Säfken, B. (2016). Smoothing parameter and model selection for general smooth models. *JASA* 111(516), 1548–1563.
