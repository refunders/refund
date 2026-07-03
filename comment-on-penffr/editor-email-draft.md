Editor email - draft (for review before sending)

To: Editorial office, Advances in Data Analysis and Classification (Editors: M. Vichi, A. Cerioli, H. A. Kestler, A. Okada, C. Weihs) - set recipient before sending
From: fabian.scheipl@googlemail.com
Subject: Comment + request for correction: Tamo Tchomgui et al. (2026), DOI 10.1007/s11634-026-00681-w

Dear Editors,

I am submitting the attached Comment on Tamo Tchomgui, Jacques, Fraysse, Barriac and Chretien (2026), "A penalized spline estimator for functional linear regression with functional response" (DOI 10.1007/s11634-026-00681-w), and asking you to consider a correction. I first raised these concerns with the authors directly and in private, offering our reproduction materials so they could correct the record themselves. The corresponding author replied only that the paper had "gone through the standard review process and is now published" and that I was "free to contact the ADAC editorial board", without engaging with any of the specific points; I am therefore doing exactly that.

I am a co-author of the pffr method against which the article compares its proposal, so I have a declared interest; that is exactly why I am raising this formally and supporting the authors' right of reply. Having reproduced the authors' own analysis, my conclusion is that the paper's central claims do not hold - and that even the results it reports for its own method were not produced by the algorithm it describes:

- The proposed estimator is not new - it is penalized-spline function-on-function regression (the pffr / functional additive mixed model approach since 2015), and the foundational references are not cited. It is in fact a strict restriction of pffr (homoscedastic-Gaussian only, no mixed effects, no working uncertainty quantification), i.e. a step backwards; the generalized FAMM framework (Scheipl, Gertheiss and Greven 2016) already contains it as a special case.
- The empirical case against pffr does not hold: the released code mis-specifies pffr (two scalar covariates entered as 50-dimensional functional terms), does not implement the settings stated in the paper's own Table 3, compares different integral estimands, and copies three competitors' numbers from another paper rather than running them. I reproduce the authors' worst pffr result exactly and show it is an artefact of this. The paper's own Table 5 already ranks pffr as the best method, contradicting its conclusion, and the one novel element - conformal prediction bands - is not implemented in the released package and attains 2-54% coverage at nominal 95% where evaluated.
- Even the results reported for the authors' own method are not produced by the algorithm the paper describes. The released code estimates PenFFR by ridge-penalized least squares (an lm() call), not the linear mixed model estimated by ReML that the paper sets out; it selects the penalty by BIC over a degenerate grid - for the four-covariate Hawaii model, only the two endpoints of the range - rather than by the cross-validation the paper states; and it omits the per-curve random intercept entirely (the curve identifier is dropped from the design, and the helper that would build the random effect is never called). The PenFFR numbers in Tables 2, 4 and 5 therefore do not correspond to the method as described.

Every uncontrolled choice in the comparison favours the authors' method (a pattern common in such studies; I do not impute intent). These are elementary, checkable issues, so a correction should be straightforward.

One further point, raised strictly as a matter of process and without any implication of impropriety on anyone's part: I note that one of the authors, Prof. Jacques, is an Associate Editor of ADAC. I would be grateful if you could confirm that the original submission was handled entirely independently of him, in line with the journal's and COPE's guidance on editor-authored submissions, and I ask that this Comment likewise be handled by editors with no connection to the authors. It was Prof. Jacques who, replying for the authors, directed me to the editorial board; given his position on it, I am writing to the Editors-in-Chief directly to avoid any conflict.

I attach the Comment and a cover letter; complete reproduction materials (data, code, the authors' own analysis script) are available to the authors and reviewers, and the full correspondence with the authors is available to you on request.

Thank you for your consideration.

Yours sincerely,
Fabian Scheipl
