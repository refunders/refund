# Cover letter to the Editors, *Advances in Data Analysis and Classification*

Dear Editors,

I am writing, with my co-authors, to submit a Comment on

> J. S. Tamo Tchomgui, J. Jacques, G. Fraysse, V. Barriac and S. Chretien (2026),
> *A penalized spline estimator for functional linear regression with functional response*,
> Advances in Data Analysis and Classification. DOI 10.1007/s11634-026-00681-w.

I am a co-author of the `pffr` method (Ivanescu, Staicu, Scheipl and Greven, 2015) and of the functional additive mixed model framework it implements, against which this article compares its proposal. I have a declared interest, which is precisely why I am raising the matter formally, in a Comment, with the authors given a full right of reply.

I do not write this lightly. Having reproduced the authors' own analysis, I have concluded that the article should not have been published in its present form, and that the case for a substantial correction — and arguably a retraction — is strong. The reasons are documented and fully reproducible:

1. **The method is not new — and is in fact a step backwards.** "FFR/PenFFR" is penalized-spline (P-spline) function-on-function regression — the approach `pffr` and the functional additive mixed model framework have implemented since 2015 — assembled from textbook components (Eilers and Marx 1996; Ruppert, Wand and Carroll 2003; Wood 2017), whose originating references are standard and are not cited. It is, moreover, a strict restriction of that method: PenFFR fits only a homoscedastic-Gaussian, least-squares model with no mixed effects and no working uncertainty quantification, whereas `pffr` handles general non-Gaussian and heteroscedastic responses and functional mixed models. The generalized functional additive mixed model framework (Scheipl, Gertheiss and Greven 2016, EJS — eqs. (1)–(4) and Table 1) already contains PenFFR's estimator, its integral/historical model, its generalized responses and functional random effects as special cases.

2. **The comparison against `pffr` is invalid.** The released code enters two scalar covariates into `pffr` as 50-dimensional functional terms; it does not implement the basis settings the paper reports in its own Table 3; the integral models being compared are different estimands; and three of the competitors in the headline Canadian-Weather table (OPFFR, FDA, FPCA) were not run at all but copied from another paper. I reproduce the authors' worst-case `pffr` number exactly and show it is an artefact of this mis-specification. Every uncontrolled choice in the comparison favours the authors' method.

3. **The conclusion is contradicted by the authors' own results.** Their Table 5 already ranks `pffr` as the most accurate method on the Hawaii Ocean data, which I reproduce.

4. **The one novel element does not exist in the software and fails where evaluated.** The conformal prediction bands are not implemented in the released package (which returns point predictions only), and the only coverage numbers reported reach 2–54% for nominal 95% intervals.

These are elementary, checkable points — the released code contradicts the paper's own Table 3, the missing references are canonical, and the claimed novelty is absent from the accompanying package. I therefore also ask that the handling of this submission be reviewed.

I attach the Comment and our complete reproduction materials (data, code, and the authors' own analysis script), which I am glad to share with the authors and reviewers. I support the authors' right of reply and would welcome their response.

Yours sincerely,

Fabian Scheipl
[affiliation, contact]
on behalf of the co-authors
