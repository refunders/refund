# Editor email — draft (for review before sending)

**To:** Editorial office, *Advances in Data Analysis and Classification* (Editors: M. Vichi, A. Cerioli, H. A. Kestler, A. Okada, C. Weihs) — set recipient before sending
**From:** fabian.scheipl@googlemail.com
**Subject:** Request for correction/retraction + Comment: Tamo Tchomgui et al. (2026), DOI 10.1007/s11634-026-00681-w

Dear Editors,

I am submitting the attached Comment on Tamo Tchomgui, Jacques, Fraysse, Barriac & Chretien (2026), "A penalized spline estimator for functional linear regression with functional response" (DOI 10.1007/s11634-026-00681-w), and asking you to consider a correction or retraction.

I am a co-author of the `pffr` method against which the article compares its proposal, so I have a declared interest; that is exactly why I am raising this formally and supporting the authors' right of reply. Having reproduced the authors' own analysis, my conclusion is that the paper's two central claims are both unfounded:

- The proposed estimator is not new — it is penalized-spline function-on-function regression (the `pffr` / functional additive mixed model approach since 2015), and the foundational references are not cited. It is in fact a strict restriction of `pffr` (homoscedastic-Gaussian only, no mixed effects, no working uncertainty quantification), i.e. a step backwards; the generalized FAMM framework (Scheipl, Gertheiss & Greven 2016) already contains it as a special case.
- The empirical case against `pffr` does not hold: the released code mis-specifies `pffr` (two scalar covariates entered as 50-dimensional functional terms), does not implement the settings stated in the paper's own Table 3, compares different integral estimands, and copies three competitors' numbers from another paper rather than running them. I reproduce the authors' worst `pffr` result exactly and show it is an artefact of this. The paper's own Table 5 already ranks `pffr` as the best method, contradicting its conclusion, and the one novel element — conformal prediction bands — is not implemented in the released package and attains 2–54% coverage at nominal 95% where evaluated.

Every uncontrolled choice in the comparison favours the authors' method. These are elementary, checkable issues, and I would also ask that the editorial handling of the submission be reviewed.

I attach the Comment and a cover letter; complete reproduction materials (data, code, the authors' own analysis script) are available to the authors and reviewers.

Thank you for your consideration.

Yours sincerely,
Fabian Scheipl
