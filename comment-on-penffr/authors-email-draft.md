# Email to the authors — draft (for review before sending)

**To:** jean-steve.tamo-tchomgui@univ-lyon2.fr, julien.jacques@univ-lyon2.fr, guillaume.fraysse@orange.com, vincent.barriac@orange.com, stephane.chretien@univ-lyon2.fr
**Cc:** (co-authors of this note)
**From:** fabian.scheipl@googlemail.com
**Subject:** Reproducing the pffr comparison in your ADAC paper — a few concerns, and an offer to help

Dear Dr Tamo Tchomgui and colleagues,

I read your recent paper, "A penalized spline estimator for functional linear regression with functional response" (ADAC, 2026), with interest — the conformal prediction-band idea is a nice one, and on the Canadian Weather data your penalized method does fit better than a well-specified `pffr`, which I am happy to acknowledge. I am one of the authors of `pffr`/`refund`, so I read the comparisons closely, reproduced them from your released code, and I wanted to raise a few concerns with you directly before considering any more formal step. My aim is genuinely constructive: I think most of these are straightforward to put right.

**1. Attribution / positioning.** The estimator — a B-spline expansion of covariates and coefficient functions with a second-derivative (curvature) penalty, reduced to a mixed model — is the penalized-spline approach that `pffr` and the functional additive mixed model framework implement. It would help readers to position it that way and to cite the originating literature: Eilers & Marx (1996, P-splines); Ruppert, Wand & Carroll (2003) and Wood (2017) for the mixed-model representation of penalized splines; Scheipl, Staicu & Greven (2015) and Scheipl, Gertheiss & Greven (2016) for the (generalized) functional additive mixed model; and Malfait & Ramsay (2003) for the historical/integral model. In particular, Scheipl, Gertheiss & Greven (2016, EJS) already contains your concurrent and integral/historical models, the tensor-product spline estimator with anisotropic curvature penalty (their eqs. (1)–(4), Table 1), generalized responses, and functional random effects.

**2. The `pffr` comparison.** Reproducing your `RD_Canada.Rmd`, I obtain your reported concurrent `pffr` ISE (≈89) — but the call enters the two scalar station coordinates (`lat`, `lon`) as `n×365` matrices, so `pffr` fits a *time-varying* coefficient with a 50-basis smooth for each (≈26 and ≈20 effective df), whereas your own method enters them as constant scalars. Entering them as scalars in `pffr` (`c(lat)+c(lon)`) drops the error substantially. A few related points: the bases in the code (`bs="cr"`, k=10/50) do not match the "100/40 cubic B-splines" stated in your Table 3; the `pffr` integral term integrates over the full range while PenFFR's integral model is historical (∫₀ᵗ), so the two are not the same estimand; and the OPFFR/FDA/FPCA figures in Table 4 are taken from Sun et al. (2018) rather than run on your split. I would be glad to share a corrected, like-for-like `pffr` script.

**3. Hawaii Ocean.** Your Table 5 already reports concurrent `pffr` as the most accurate method (0.52×10⁻², ahead of all PenFFR variants), which I reproduce exactly. The surrounding text ("our method once again … outperforms all the other methods") seems to overstate this and could be adjusted.

**4. The prediction bands.** The conformal/optimal-transport procedure of §4 does not appear to be exercised in the released package — `pred.PenFFR`/`pred.penffr*` return point predictions, and the OT-quantile helpers are not called — and the coverage reported (Table 2, Fig. 5) is well below the nominal level. If I have missed where this is implemented, please point me to it; otherwise it may be worth either including a working implementation or softening the claims. (`pffr`/`mgcv`, for what it's worth, already provide calibrated interval estimates and non-Gaussian/mixed-model fits, if those are useful to you.)

I have a fully reproducible set of scripts (data, your own analysis file, and the corrected `pffr` runs) that I am happy to send so you can check all of the above independently. I think a correction addressing the attribution and the comparison would resolve my concerns. I did want to be transparent that, absent a correction, I would likely submit a short Comment to ADAC — but I would much rather sort this out collegially and first, and I am very open to being corrected myself if I have misread anything.

With best regards,
Fabian Scheipl
[affiliation, contact]
