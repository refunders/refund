# Cover letter to the Editors, *Advances in Data Analysis and Classification*

Dear Editors,

I am writing, with my co-authors, to submit a Comment on the article

> J. S. Tamo Tchomgui, J. Jacques, G. Fraysse, V. Barriac and S. Chretien (2026),
> *A penalized spline estimator for functional linear regression with functional response*,
> Advances in Data Analysis and Classification. DOI 10.1007/s11634-026-00681-w.

I am a co-author of the `pffr` method (Ivanescu, Staicu, Scheipl and Greven, 2015) and of the functional additive mixed model framework it implements, against which the article compares its proposal. I have a clear interest in the matter, which is precisely why I believe a Comment, with the authors given the right of reply, is the appropriate venue.

Our Comment raises three substantive points, each documented and fully reproducible:

1. **Novelty.** The estimator proposed in the article — a cubic B-spline expansion with a roughness penalty on second derivatives, reduced to a (mixed) linear model — is the established penalized-spline approach to function-on-function regression, i.e. the methodology that `pffr` and the functional additive mixed model framework have implemented since 2015, built on P-splines (Eilers and Marx, 1996) and the mixed-model representation of penalized splines (Ruppert, Wand and Carroll, 2003; Wood, 2017). We document the correspondence component by component.

2. **Attribution.** The article omits the foundational references for the methodology it re-derives. We list the specific works that should be cited.

3. **Empirical comparison.** The comparison with `pffr` is mis-specified and internally inconsistent, and the released analysis code does not implement the settings described in the paper. We reproduce the authors' headline result for `pffr` exactly, show that a large part of its reported error is an artefact of mis-specification (scalar covariates entered as high-dimensional functional terms), and note that the authors' own Hawaii Ocean table already ranks `pffr` as the best method — contradicting the article's stated conclusion.

Our intent is corrective rather than adversarial: we ask for proper attribution of prior work and a corrected, like-for-like comparison. We would welcome the authors' response and are glad to provide our complete reproduction materials (data, code and the authors' own analysis script) to the authors and to the reviewers.

Thank you for considering this Comment.

Yours sincerely,

Fabian Scheipl
[affiliation, contact]
on behalf of the co-authors
