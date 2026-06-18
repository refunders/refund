# E4 — fastFMM / FUI comparison: findings & decision to ABANDON

**Decision (2026-06-18, user):** Do **not** include a fastFMM/FUI comparison in the
paper. This note records what was investigated and why we dropped it, plus the
reusable byproducts worth keeping.

---

## 1. Why E4 was proposed
`fastFMM` (FUI = Fast Univariate Inference; Cui, Leroux, Smirnova & Crainiceanu
2022, *JCGS* 31(1):219–230) is a popular package for functional regression with
analytic CIs. We wanted it as a competitor CI method alongside the paper's
default / HC / cluster / cl2 sandwich intervals for `pffr()`.

## 2. Finding — FUI cannot be fit on the paper's DGPs
The paper's Study 1/2 have **one functional curve per subject** (independent
curves). FUI is a functional **mixed-model** method and requires a subject random
intercept `(1|id)`; with one curve per subject that RE is unidentifiable (n levels
= n observations) and `lme4` rejects it before variance estimation
(`override_zero_var=TRUE` does not help). FUI's own example uses `refund::DTI`
(multiple visits per subject). FUI also only does function-on-**scalar**, never the
`ff` surface that is the paper's focus. → A fair comparison would require a
**different**, repeated-measures sub-study.

## 3. The repeated-measures sub-study we built (and its problems)
We built a longitudinal DGP where FUI applies:

    Y_ij(t) = β0(t) + zlin_ij·β1(t) + b_i(t) + ε_ij(t),   i=1..N subjects, j=1..J visits,

with a functional random intercept `b_i(t)=Σ_{k=1}^3 ξ_ik ψ_k(t)`, iid pointwise
noise ε, and zero-mean smooth β0/β1. Fit by `pffr(Y ~ zlin + s(id, bs="re"))` and
`fui(Y ~ zlin + (1|id))`. Two issues surfaced (first via review, then a council):

- An initial driver (a) silently swapped the planned functional RE for a *scalar*
  one to dodge a (non-existent) pffr fitting issue, and (b) clustered the sandwich
  **by curve**. Both were wrong (see §4).
- `pffr(Y ~ zlin + s(id, bs="re"))` fits a functional RE fine (~13 s); the
  "centering problem" was a non-issue.

## 4. Council verdict (Codex + Gemini + Claude, UNANIMOUS)

### (a) Cluster by SUBJECT, not by curve
Clustering the sandwich by curve is **wrong**; cluster by subject. The tempting
argument — "the scores use *conditional* residuals `y − Xβ̂ − Zb̂`, and given the
true `b_i` the J curves are independent, so curve-wise is fine" — fails because the
sandwich plugs in the **estimated** BLUP `b̂_i`, which is computed from *all J
curves* of the subject. The conditional residual therefore carries the shared term
`−Z(b̂_i − b_i)`, correlating the curve-level scores **within** a subject. Cluster-
robust theory requires independence *across* clusters at the level of the
independent sampling unit = the subject; the frequentist CI target resamples
subjects, not curves. (Liang & Zeger 1986; Cameron & Miller 2015; MacKinnon,
Nielsen & Webb 2022.) Conditional-vs-marginal residual choice does **not** change
the unit of independence.

Empirical fingerprint: **cluster-by-curve is numerically identical to the
model-based default** — curve-level clustering captures none of the subject-level
correlation, so it just reproduces the model SE.

### (b) The binding constraint is the PENALIZED BREAD, not the cluster level
Even **subject-clustered** pffr undercovers badly relative to FUI. All three legs
attribute this to `V_p = (XᵀWX + S)⁻¹`: the RE penalty `S` (REML estimate of the
RE precision) **over-shrinks the between-subject variance** — especially with the
Kronecker-spline `s(id,bs="re")` RE not matching the rank-3 truth and only N≈50
subjects — so the fixed-effect SEs collapse. The sandwich is applied on top of an
already-shrunk estimator; no choice of *meat* (cluster level) can recover variance
the *bread* has thrown away. Gemini also notes the likelihood over-counts
information (L grid points treated as independent). Possible fixes the council
floated: unpenalized-bread `(XᵀWX)⁻¹` for the fixed-effect block, a marginal-
residual GEE-style subject-clustered sandwich, a Kenward–Roger variance-component
correction, or an FPCA random effect that can represent the true covariance.

### (c) FUI vs pffr is not apples-to-apples
FUI targets the **marginal** Var(β̂(t)) directly via per-point LMEs (+ subject-level
inference / joint bootstrap bands); the pffr sandwich is a pointwise, penalized-
*joint* object. Different variance functionals and interval types (pointwise vs
joint), so a raw coverage/width comparison conflates several things. FUI being
near-nominal at *smaller* width than subject-clustered pffr means it is better
**calibrated** for this DGP, not merely more conservative.

### (d) Monte-Carlo size
All exploratory numbers below are R = 8–12 reps (MC-SE on a 0.2 coverage ≈ 0.14):
useful as smoke tests, **not citable**. A real study needs R ≥ 500.

## 5. Empirical evidence (exploratory; functional-RE DGP, N=50, J=4, L=50)

β1(t) pointwise coverage, nominal 0.90:

| method | coverage | mean width |
|--------|----------|------------|
| pffr default (Bayesian) | 0.18 | 0.036 |
| cluster, **by curve** | 0.18 | 0.036 |  ← identical to default |
| cluster, **by subject** | 0.34 | 0.079 |
| cl2, by curve | 0.36 | 0.078 |
| cl2, **by subject** | 0.59 | 0.179 |
| fastFMM (fui) | 0.86 | 0.065 |

A partial 10-rep/4-cell pilot (before the run was stopped) corroborated the
ordering across cells: pffr_default ≈ 0.06–0.28, pffr_cluster ≈ 0.19–0.60,
pffr_cl2 ≈ 0.49–0.86 (higher at larger RE variance), fastfmm ≈ 0.82–0.90.

## 6. Why we abandoned it
- The comparison is tangential to the paper's core contribution (sandwich CIs for
  **function-on-function** regression on **independent** curves).
- Done honestly it requires a whole separate repeated-measures study, an
  apples-to-apples redesign (matched estimands/interval types), and it opens a
  methodological investigation into pffr's penalized-RE variance under-propagation
  (the "bread" problem) that is a paper of its own, not a competitor row.
- FUI addresses a *different data structure* (longitudinal/repeated functional
  measures); it simply isn't a competitor on the paper's setting.

## 7. Reusable byproducts (keep)
- **`coef.pffr(..., cluster=)`** — new argument for subject/nested-level sandwich
  clustering (commit `b0cee3d8`, with regression test; suite green). General,
  backward-compatible (NULL = old by-curve behaviour), and correct per the council.
  Now unused by any paper study, but a genuine, tested improvement — recommend
  **keep** (revert only if you want the package surface minimal).
- **LRZ `Rfast`/`fastFMM` install** fixed by Codex (`.Rprofile` stdout banner
  removed + `<numeric>` patch + `PKG_CPPFLAGS +=`); patched tarball at
  `~/rfast-src/Rfast-patched2.tar.gz` on LRZ. Harmless to leave.

## 8. Abandoned artifacts (in `ci-experiments`, history only)
- Drivers `ci-benchmark/sim-study-fastfmm-extension.R`,
  `sim-study-fastfmm-longitudinal.R` (the functional-RE revision was left
  uncommitted in the working tree when the run was stopped).
- Data dirs `ci-benchmark/study4-longitudinal/`, `study2-fastfmm/`.
- Helper script `~/fda/check-e4-dgp.R` (+ `.pdf`).
These can stay as history or be cleaned up; nothing in the paper depends on them.

## References (council-cited)
- Liang & Zeger (1986), *Biometrika* — cluster-robust GEE.
- Cameron & Miller (2015), *J. Human Resources* — practitioner's guide to cluster-robust inference.
- MacKinnon, Nielsen & Webb (2022), arXiv:2205.03285 — cluster-robust inference.
- Cui, Leroux, Smirnova & Crainiceanu (2022), *JCGS* 31(1):219–230 — FUI / fastFMM.
- Wood, Pya & Säfken (2016), *JASA* — smoothing-parameter uncertainty (`Vc`).
- Nychka (1988); Marra & Wood (2012) — Bayesian CI coverage for penalized splines.
- Nobre & Singer (2007), *Biometrical J.* — conditional residuals in LMMs.
