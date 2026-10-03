#' Penalized flexible functional regression
#'
#' Implements additive regression for functional and scalar covariates and
#' functional responses. This function is a wrapper for \code{mgcv}'s
#' \code{\link[mgcv]{gam}} and its siblings to fit models of the general form
#' \cr \eqn{E(Y_i(t)) = g(\mu(t) + \int X_i(s)\beta(s,t)ds + f(z_{1i}, t) +
#' f(z_{2i}) + z_{3i} \beta_3(t) + \dots )}\cr with a functional (but not
#' necessarily continuous) response \eqn{Y(t)}, response function \eqn{g},
#' (optional) smooth intercept \eqn{\mu(t)}, (multiple) functional covariates
#' \eqn{X(t)} and scalar covariates \eqn{z_1}, \eqn{z_2}, etc.
#'
#' @section Details: The routine can estimate \enumerate{ \item linear
#'   functional effects of scalar (numeric or factor) covariates that vary
#'   smoothly over \eqn{t} (e.g. \eqn{z_{1i} \beta_1(t)}, specified as
#'   \code{~z1}), \item nonlinear, and possibly multivariate functional effects
#'   of (one or multiple) scalar covariates \eqn{z} that vary smoothly over the
#'   index \eqn{t} of \eqn{Y(t)} (e.g. \eqn{f(z_{2i}, t)}, specified in the
#'   \code{formula} simply as \code{~s(z2)}) \item (nonlinear) effects of scalar
#'   covariates that are constant over \eqn{t} (e.g. \eqn{f(z_{3i})}, specified
#'   as \code{~c(s(z3))}, or \eqn{\beta_3 z_{3i}}, specified as \code{~c(z3)}),
#'   \item function-on-function regression terms (e.g. \eqn{\int
#'   X_i(s)\beta(s,t)ds}, specified as \code{~ff(X, yindex=t, xindex=s)}, see
#'   \code{\link{ff}}). Terms given by \code{\link{sff}} and \code{\link{ffpc}}
#'   provide nonlinear and FPC-based effects of functional covariates,
#'   respectively. \item concurrent effects of functional covariates \code{X}
#'   measured on the same grid as the response  are specified as follows:
#'   \code{~s(x)} for a smooth, index-varying effect \eqn{f(X(t),t)}, \code{~x}
#'   for a linear index-varying effect \eqn{X(t)\beta(t)}, \code{~c(s(x))} for a
#'   constant nonlinear effect \eqn{f(X(t))}, \code{~c(x)} for a constant linear
#'   effect \eqn{X(t)\beta}. \item Smooth functional random intercepts
#'   \eqn{b_{0g(i)}(t)} for a grouping variable \code{g} with levels \eqn{g(i)}
#'   can be specified via \code{~s(g, bs="re")}), functional random slopes
#'   \eqn{u_i b_{1g(i)}(t)} in a numeric variable \code{u} via \code{~s(g, u,
#'   bs="re")}). Scheipl, Staicu, Greven (2013) contains code examples for
#'   modeling correlated functional random intercepts using
#'   \code{\link[mgcv]{mrf}}-terms. } Use the \code{c()}-notation to denote
#'   model terms that are constant over the index of the functional response.\cr
#'
#'   Internally, univariate smooth terms without a \code{c()}-wrapper are
#'   expanded into bivariate smooth terms in the original covariate and the
#'   index of the functional response. Bivariate smooth terms (\code{s(), te()}
#'   or \code{t2()}) without a \code{c()}-wrapper are expanded into trivariate
#'   smooth terms in the original covariates and the index of the functional
#'   response. Linear terms for scalar covariates or categorical covariates are
#'   expanded into varying coefficient terms, varying smoothly over the index of
#'   the functional response. For factor variables, a separate smooth function
#'   with its own smoothing parameter is estimated for each level of the
#'   factor.\cr \cr The marginal spline basis used for the index of the the
#'   functional response is specified via the \emph{global} argument
#'   \code{bs.yindex}. If necessary, this can be overriden for any specific term
#'   by supplying a \code{bs.yindex}-argument to that term in the formula, e.g.
#'   \code{~s(x, bs.yindex=list(bs="tp", k=7))} would yield a tensor product
#'   spline over \code{x} and the index of the response in which the marginal
#'   basis for the index of the response are 7 cubic thin-plate spline functions
#'   (overriding the global default for the basis and penalty on the index of
#'   the response given by the \emph{global} \code{bs.yindex}-argument).\cr Use
#'   \code{~-1 + c(1) + ...} to specify a model with only a constant and no
#'   functional intercept. \cr
#'
#'   The functional covariates have to be supplied as a \eqn{n} by <no. of
#'   evaluations> matrices, i.e. each row is one functional observation. For
#'   data on a regular grid, the functional response is supplied in the same
#'   format, i.e. as a matrix-valued entry in \code{data},  which can contain
#'   missing values.\cr
#'
#'   If the functional responses are \emph{sparse or irregular} (i.e., not
#'   evaluated on the same evaluation points across all observations), the
#'   \code{ydata}-argument can be used to specify the responses: \code{ydata}
#'   must be a \code{data.frame} with 3 columns called \code{'.obs', '.index',
#'   '.value'} which specify which curve the point belongs to
#'   (\code{'.obs'}=\eqn{i}), at which \eqn{t} it was observed
#'   (\code{'.index'}=\eqn{t}), and the observed value
#'   (\code{'.value'}=\eqn{Y_i(t)}). Note that the vector of unique sorted
#'   entries in \code{ydata$.obs} must be equal to \code{rownames(data)} to
#'   ensure the correct association of entries in \code{ydata} to the
#'   corresponding rows of \code{data}. For both regular and irregular
#'   functional responses, the model is then fitted with the data in long
#'   format, i.e., for data on a grid the rows of the matrix of the functional
#'   response evaluations \eqn{Y_i(t)} are stacked into one long vector and the
#'   covariates are expanded/repeated correspondingly. This means the models get
#'   quite big fairly fast, since the effective number of rows in the design
#'   matrix is number of observations times number of evaluations of \eqn{Y(t)}
#'   per observation.\cr
#'
#'   When a \code{rho} argument (as defined in \code{\link[mgcv]{bam}}) is supplied,
#'   \code{pffr} automatically constructs the required \code{AR.start} indicator
#'   so that residuals can follow an AR(1) process along the functional response.
#'   This facility is available for \code{algorithm = "bam"} fits that use
#'   \code{method = "fREML"}. For non-Gaussian families, \code{discrete = TRUE}
#'   is required; \code{pffr} will set this automatically if needed, mirroring
#'   the constraints documented in \code{\link[mgcv]{bam}}.\cr
#'
#'   The following automatic defaults apply when \code{rho} is supplied:
#'   \itemize{
#'     \item \code{algorithm} auto-switches to \code{"bam"} (errors if the user
#'       explicitly requests a non-\code{bam} algorithm)
#'     \item \code{method} auto-switches to \code{"fREML"} (errors if the user
#'       explicitly supplies a different method)
#'     \item \code{discrete} is auto-set to \code{TRUE} for non-Gaussian families
#'       (errors if the user explicitly supplies \code{discrete = FALSE})
#'     \item \code{sandwich} resolves to \code{FALSE} (model-based
#'       intervals, with a message) when left at its default: the sandwich
#'       assumes working-independence scores and is not valid for the
#'       AR(1)-whitened fit. An explicit \code{sandwich = TRUE} errors.
#'   }
#'   Explicit user overrides that conflict with these constraints will produce
#'   an informative error.\cr
#'
#'   Note that \code{pffr} does not use \code{mgcv}'s default identifiability
#'   constraints (i.e., \eqn{\sum_{i,t} \hat f(z_i, x_i, t) = 0} or
#'   \eqn{\sum_{i,t} \hat f(x_i, t) = 0}) for tensor product terms whose
#'   marginals include the index \eqn{t} of the functional response.  Instead,
#'   \eqn{\sum_i \hat f(z_i, x_i, t) = 0} for all \eqn{t} is enforced, so that
#'   effects varying over \eqn{t} can be interpreted as local deviations from
#'   the global functional intercept. This is achieved by using
#'   \code{\link[mgcv]{ti}}-terms with a suitably modified \code{mc}-argument.
#'   Note that this is not possible if \code{algorithm='gamm4'} since only
#'   \code{t2}-type terms can then be used and these modified constraints are
#'   not available for \code{t2}. We recommend using centered scalar covariates
#'   for terms like \eqn{z \beta(t)} (\code{~z}) and centered functional
#'   covariates with \eqn{\sum_i X_i(t) = 0} for all \eqn{t} in \code{ff}-terms
#'   so that the global functional intercept can be interpreted as the global
#'   mean function.
#'
#'   The \code{family}-argument can be used to specify all of the response
#'   distributions and link functions described in
#'   \code{\link[mgcv]{family.mgcv}}. Note that  \code{family = "gaulss"} is
#'   treated in a special way: Users can supply the formula for the variance by
#'   supplying a special argument \code{varformula}, but this is not modified in
#'   the way that the \code{formula}-argument is but handed over to the fitter
#'   directly, so this is for expert use only. If \code{varformula} is not
#'   given, \code{pffr} will use the parameters from argument \code{bs.int} to
#'   define a spline basis along the index of the response, i.e., a smooth
#'   variance function over $t$ for responses $Y(t)$.
#'
#' @param formula a formula with special terms as for \code{\link[mgcv]{gam}},
#'   with additional special terms \code{\link{ff}(), \link{sff}(),
#'   \link{ffpc}(), \link{pcre}()} and \code{c()}.
#' @param yind a vector with length equal to the number of columns of the matrix
#'   of functional responses giving the vector of evaluation points \eqn{(t_1,
#'   \dots ,t_{G})}. If not supplied, \code{yind} is set to
#'   \code{1:ncol(<response>)}.
#' @param algorithm the name of the function used to estimate the model.
#'   Defaults to \code{\link[mgcv]{gam}} if the matrix of functional responses
#'   has less than \code{2e5} data points and to \code{\link[mgcv]{bam}} if not.
#'   \code{'\link[mgcv]{gamm}'}, \code{'\link[gamm4]{gamm4}'} and
#'   \code{'\link[mgcv]{jagam}'} are valid options as well. See Details for
#'   \code{'\link[gamm4]{gamm4}'} and \code{'\link[mgcv]{jagam}'}.
#' @param data an (optional) \code{data.frame} containing the data. Can also be
#'   a named list for regular data. Functional covariates have to be supplied as
#'   <no. of observations> by <no. of evaluations> matrices, i.e. each row is
#'   one functional observation.
#' @param ydata an (optional) \code{data.frame} supplying functional responses
#'   that are not observed on a regular grid. See Details.
#' @param method Defaults to \code{"REML"}-estimation, including of unknown
#'   scale. If \code{algorithm="bam"}, the default is switched to
#'   \code{"fREML"}. See \code{\link[mgcv]{gam}} and \code{\link[mgcv]{bam}} for
#'   details.
#' @param bs.yindex a named (!) list giving the parameters for spline bases on
#'   the index of the functional response. Defaults to \code{list(bs="ps", k=5,
#'   m=c(2, 1))}, i.e. 5 cubic B-splines bases with first order difference
#'   penalty.
#' @param bs.int a named (!) list giving the parameters for the spline basis for
#'   the global functional intercept. Defaults to \code{list(bs="ps", k=20,
#'   m=c(2, 1))}, i.e. 20 cubic B-splines bases with first order difference
#'   penalty.
#' @param tensortype which typ of tensor product splines to use. One of
#'   "\code{\link[mgcv]{ti}}" or "\code{\link[mgcv]{t2}}", defaults to
#'   \code{ti}. \code{t2}-type terms do not enforce the more suitable special
#'   constraints for functional regression, see Details.
#' @param sandwich Covariance for standard errors and confidence intervals.
#'   \code{TRUE} (default): the curve-clustered CL2 sandwich, computed at fit
#'   time; see the section \sQuote{Inference}. \code{FALSE}: the model-based
#'   (Bayesian posterior) covariance of the working-independence fit with
#'   Gaussian critical values. Model-based intervals assume independent errors
#'   within each curve and are too narrow under within-curve dependence, more
#'   so on dense grids; use them only if that dependence is known to be
#'   absent. The character values of refund 0.1-40 are deprecated:
#'   \code{"cl2"} means \code{TRUE}, \code{"none"} means \code{FALSE}, and
#'   \code{"cluster"} and \code{"hc"} are no longer available and give the CL2
#'   sandwich.
#' @param cluster Optional grouping with one nonmissing entry per curve,
#'   evaluated in \code{data}. With several curves per subject, cluster by
#'   subject: the CL2 sandwich (and, for \code{method = "NCV"}, the
#'   neighbourhoods) then treat subjects, not curves, as the independent units.
#'   This leaves fewer clusters and does not address dependence between
#'   subjects. Supported for densely observed responses only.
#' @param ncv_blocks For \code{method = "NCV"}, \code{"cluster"} (default)
#'   leaves out each whole curve, or each group of curves defined by
#'   \code{cluster}. \code{"point"} uses leave-one-point-out NCV for comparisons.
#'   A supplied \code{nei} in \code{...} takes precedence; its indices
#'   must refer to retained model-frame rows, after response omissions.
#' @section Inference:
#' By default (\code{sandwich = TRUE}) standard errors and intervals use the
#' curve-clustered CL2 sandwich: per-curve (or per-\code{cluster}) score sums
#' of the working-independence fit, each adjusted by the full block of the
#' penalized hat matrix \eqn{H},
#' \eqn{A_g = \{(I - H)^2\}_{gg}^{-1/2}}{A_g = ((I - H)^2)_gg^(-1/2)}
#' (Bell and McCaffrey), in the Bayesian form that adds \eqn{V_p - V_e} to the
#' sandwich. Pointwise intervals from \code{\link{coef.pffr}},
#' \code{\link{predict.pffr}} and \code{\link{pffr_predict_ci}} use
#' Satterthwaite critical values, i.e. \eqn{t_\nu}{t_nu} quantiles with a
#' separate working-model moment df \eqn{\nu}{nu} for every evaluation point.
#' Fit by REML (the default) for inference. In simulations with Gaussian, count
#' and binary responses, dependent and independent errors and 20 to 100
#' curves, these intervals covered the mean, coefficient surfaces and
#' coefficient functions of scalar covariates close to the nominal level; they
#' are wide for coefficient surfaces. Only pointwise intervals were evaluated.
#' The following are exceptions or caveats:
#' \itemize{
#'   \item With fewer than 40 curves or clusters \code{pffr()} warns once at
#'     fit time (class \code{"pffr_small_G_warning"}): the procedure was
#'     evaluated down to 20 curves; treat intervals as approximate.
#'   \item For binary responses a message (once per session) states that
#'     intervals for the functional intercept and for coefficient functions of
#'     scalar covariates can undercover, also with Satterthwaite critical
#'     values. Intervals for count responses whose curves are misregistered
#'     (observed on their own clocks) undercover as well.
#'   \item Families: exponential-family responses (including quasi families)
#'     use their exact scores, \code{gaulss} and \code{scat} exact
#'     (two-block) scores. Other extended families (\code{nb}, \code{tw},
#'     \code{betar}, \code{ocat}, ...) use the exponential-family
#'     working-residual approximation to the score, with a warning once per
#'     session. Families without a cluster score path (those defining their
#'     own \code{family$sandwich}, e.g. \code{multinom}), the
#'     \code{"gamm"}/\code{"gamm4"} algorithms and fits with a single curve
#'     fall back to model-based intervals with a message (a warning if
#'     \code{sandwich = TRUE} was given explicitly).
#'   \item The full-block adjustment grows with the number of curves, the
#'     curve length and the basis dimension, but stays well below the cost of
#'     the fit (a few seconds for 100 or 300 curves of 60 points with 104
#'     coefficients, against 15 and 60 seconds for the REML fit); there is no
#'     upper limit on the number of curves.
#'     \code{\link{predict.pffr}} and \code{\link{pffr_predict_ci}} compute
#'     Satterthwaite df for at most
#'     \code{getOption("refund.pffr_satterthwaite_max_points", 1e4)}
#'     prediction points and use Gaussian critical values beyond that, with a
#'     message. Where a df is undefined (a contrast with zero sampling
#'     variance) the Gaussian critical value is used, with a message.
#'   \item \code{\link{summary.pffr}} reports mgcv's model-based tests;
#'     \code{\link{plot.pffr}} draws mgcv's bands of \eqn{\pm 2}{+/- 2} CL2
#'     standard errors. Use \code{\link{coef.pffr}} with
#'     \code{ci = "pointwise"} for the intervals described here.
#' }
#' Storage: the fit's model-based (Bayesian posterior) covariance matrices
#' \code{$Vp}, \code{$Vc} and \code{$Ve} are left exactly as
#' \code{\link[mgcv]{gam}} produced them; the CL2 covariance is stored in
#' \code{$pffr$Vsandwich} (with metadata in \code{$pffr$sandwich_info}) and
#' returned by \code{\link{vcov.pffr}}. mgcv's generics applied to the
#' underlying gam report model-based uncertainty.
#' @section Neighbourhood cross-validation:
#' With mgcv >= 1.9, \code{method = "NCV"} targets prediction of whole new
#' curves (or clusters). Leaving out single points can undersmooth severely
#' when errors within curves are dependent. Under within-curve dependence REML
#' undersmooths and its estimate of a coefficient surface is noisy, while
#' curve-blocked NCV smooths more and estimates the surface's shape better.
#' Use the NCV fit's estimate to describe the shape of coefficient surfaces,
#' and a REML fit of the same model with the default CL2 intervals for
#' inference: intervals centred at the NCV estimate undercover, so do not use
#' the NCV fit's intervals for inference (\code{\link{coef.pffr}} and the
#' prediction functions say so once per session). For the intercept, effects
#' of scalar covariates and fitted means, the REML fit gives both estimate and
#' intervals.
#' Blocked NCV typically costs roughly
#' 3--5 times as much as REML, depending on the model and block sizes.
#' Supported backends are \code{algorithm = "gam"} and \code{algorithm = "bam"}.
#' An unspecified algorithm selects gam even for large data. For bam, pffr sets
#' \code{discrete = TRUE}; explicitly supplying \code{discrete = FALSE} errors.
#' Discrete NCV inverts a matrix of side equal to the neighbourhood size for
#' each neighbourhood. It is usually slower than gam for whole-curve blocks;
#' \code{nei$sample} can reduce the cost (see \code{\link[mgcv]{bam}}).
#' Nonzero \code{rho} (AR1 residuals) is unavailable with NCV and errors.
#' With bam, all neighbourhoods must have the same size: mgcv's discrete NCV
#' (checked for 1.9-3 and 1.9-4) corrupts memory when they differ, so curves
#' with missing responses or irregular grids stop with an error; use
#' \code{algorithm = "gam"} for those.
#' gamm and gamm4 do not provide this NCV selector. \code{subset} is unsupported.
#' Missing dense response values are allowed with \code{na.action = na.omit};
#' \code{na.exclude} padding and additional model-frame omissions cause an
#' error. Sparse responses use the supplied \code{ydata} row order; custom
#' \code{cluster} grouping for sparse responses remains unsupported.
#' For both backends, user-supplied \code{nei} indices refer to retained
#' model-frame rows. Bam requires the \code{a/ma/d/md} neighbourhood names
#' (not the older gam names \code{k/m/i/mi}); its documented defaults for
#' omitted prediction indices apply. Pffr removes missing response rows before
#' bam's fresh setup, including corresponding weights and offsets, so bam
#' receives these retained-row indices without a second NA translation.
#' Smoothing parameters supplied via \code{sp} must be either all fixed or all
#' free: mgcv (up to at least 1.9-4) fails in its NCV optimizer when only some
#' are fixed.
#' Generated neighbourhoods use \code{jackknife = FALSE}, retaining model-based
#' covariance for the sandwich machinery. The first NCV call checks that the
#' installed mgcv honours blocks, cached separately for each backend. The fit
#' records \code{blocks}, \code{n_blocks}
#' and \code{nei} in \code{fit$pffr$ncv}. Model-based covariance falls back to
#' \code{Vp} when smoothing-parameter uncertainty covariance is unavailable.
#' @param ... additional arguments that are valid for \code{\link[mgcv]{gam}},
#'   \code{\link[mgcv]{bam}}, \code{'\link[gamm4]{gamm4}'} or
#'   \code{'\link[mgcv]{jagam}'}. \code{subset} is not implemented.
#' @return A fitted \code{pffr}-object, which is a
#'   \code{\link[mgcv]{gam}}-object with some additional information in an
#'   \code{pffr}-entry. If \code{algorithm} is \code{"gamm"} or \code{"gamm4"},
#'   only the \code{$gam} part of the returned list is modified in this way.\cr
#'   Available methods/functions to postprocess fitted models:
#'   \code{\link{summary.pffr}}, \code{\link{plot.pffr}},
#'   \code{\link{coef.pffr}}, \code{\link{fitted.pffr}},
#'   \code{\link{residuals.pffr}}, \code{\link{predict.pffr}},
#'   \code{\link{model.matrix.pffr}},  \code{\link{qq.pffr}},
#'   \code{\link{pffr.check}}.\cr If \code{algorithm} is \code{"jagam"}, only
#'   the location of the model file and the usual
#'   \code{\link[mgcv]{jagam}}-object are returned, you have to run the sampler
#'   yourself.\cr
#' @author Fabian Scheipl, Sonja Greven
#' @seealso \code{\link[mgcv]{smooth.terms}} for details of \code{mgcv} syntax
#'   and available spline bases and penalties.
#' @references Ivanescu, A., Staicu, A.-M., Scheipl, F. and Greven, S. (2015).
#'   Penalized function-on-function regression. Computational Statistics,
#'   30(2):539--568. \url{https://biostats.bepress.com/jhubiostat/paper254/}
#'
#'   Scheipl, F., Staicu, A.-M. and Greven, S. (2015). Functional Additive Mixed
#'   Models. Journal of Computational & Graphical Statistics, 24(2): 477--501.
#'   \url{ https://arxiv.org/abs/1207.5947}
#'
#'   F. Scheipl, J. Gertheiss, S. Greven (2016):  Generalized Functional Additive Mixed Models,
#'   Electronic Journal of Statistics, 10(1), 1455--1492.
#'   \url{https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-10/issue-1/Generalized-functional-additive-mixed-models/10.1214/16-EJS1145.full}
#' @export
#' @importFrom mgcv ti jagam gam gam.fit3 bam gamm
#' @importFrom gamm4 gamm4
#' @importFrom lme4 lmer
#' @examples
#' ###############################################################################
#' # univariate model:
#' # Y(t) = f(t)  + \int X1(s)\beta(s,t)ds + eps
#' set.seed(2121)
#' data1 <- pffr_simulate(Y ~ ff(X1), n=40)
#' t <- attr(data1, "yindex")
#' s <- attr(data1, "xindex")
#' m1 <- pffr(Y ~ ff(X1, xind=s), yind=t, data=data1)
#' summary(m1)
#' plot(m1, pages=1)
#'
#' \dontrun{
#' ###############################################################################
#' # multivariate model:
#' # E(Y(t)) = \beta_0(t)  + \int X1(s)\beta_1(s,t)ds + xlin \beta_3(t) +
#' #        f_1(xte1, xte2) + f_2(xsmoo, t) + \beta_4 xconst
#' data2 <- pffr_simulate(Y ~ ff(X1) + xlin + c(te(xte1, xte2)) + s(xsmoo) + c(xconst), n=200)
#' t <- attr(data2, "yindex")
#' s <- attr(data2, "xindex")
#' m2 <- pffr(Y ~  ff(X1, xind=s) + #linear function-on-function
#'                 xlin  +  #varying coefficient term
#'                 c(te(xte1, xte2)) + #bivariate smooth term in xte1 & xte2, const. over Y-index
#'                 s(xsmoo) + #smooth effect of xsmoo varying over Y-index
#'                 c(xconst), # linear effect of xconst constant over Y-index
#'         yind=t,
#'         data=data2)
#' summary(m2)
#' plot(m2)
#' str(coef(m2))
#' # convenience functions:
#' preddata <- pffr_simulate(Y ~ ff(X1) + xlin + c(te(xte1, xte2)) + s(xsmoo) + c(xconst), n=20)
#' str(predict(m2, newdata=preddata))
#' str(predict(m2, type="terms"))
#' cm2 <- coef(m2)
#' cm2$pterms
#' str(cm2$smterms, 2)
#' str(cm2$smterms[["s(xsmoo)"]]$coef)
#'
#' #############################################################################
#' # sparse data (80% missing on a regular grid):
#' set.seed(88182004)
#' data3 <- pffr_simulate(Y ~ 1 + s(xsmoo), n=100, propmissing=0.8)
#' t <- attr(data3, "yindex")
#' m3.sparse <- pffr(Y ~ s(xsmoo), data=data3$data, ydata=data3$ydata, yind=t)
#' summary(m3.sparse)
#' plot(m3.sparse,pages=1)
#' }
pffr <- function(
  formula,
  yind,
  data = NULL,
  ydata = NULL,
  algorithm = NA,
  method = "REML",
  tensortype = c("ti", "t2"),
  bs.yindex = list(bs = "ps", k = 5, m = c(2, 1)),
  bs.int = list(bs = "ps", k = 20, m = c(2, 1)),
  sandwich = TRUE,
  cluster = NULL,
  ncv_blocks = c("cluster", "point"),
  ...
) {
  call <- match.call()
  pffr_check_removed_args(list(...), "pffr")
  ncv_blocks <- match.arg(ncv_blocks)
  if (identical(method, "NCV")) {
    dots <- list(...)
    if (is.na(algorithm)) algorithm <- "gam"
    if (!algorithm %in% c("gam", "bam")) {
      stop('pffr NCV supports only algorithm = "gam" or "bam".', call. = FALSE)
    }
    if (!is.null(dots$rho) && !identical(as.numeric(dots$rho), 0)) {
      stop(
        "Nonzero `rho` (AR1 residuals) is unavailable with NCV.",
        call. = FALSE
      )
    }
    if (algorithm == "bam") {
      if ("discrete" %in% names(dots) && !isTRUE(dots$discrete)) {
        stop(
          'pffr NCV with algorithm = "bam" requires discrete = TRUE.',
          call. = FALSE
        )
      }
      call$discrete <- TRUE
    } else if (!is.null(dots$discrete) && !identical(dots$discrete, FALSE)) {
      stop(
        'pffr NCV with discrete = TRUE requires explicit algorithm = "bam".',
        call. = FALSE
      )
    }
  }
  cluster_value <- if (missing(cluster)) NULL else
    eval(
      substitute(cluster),
      envir = data %||% parent.frame(),
      enclos = parent.frame()
    )
  tensortype <- as.symbol(match.arg(tensortype))
  sandwich_missing <- missing(sandwich)
  use_sandwich <- pffr_sandwich_arg(sandwich)
  yind_missing <- missing(yind)

  prep <- pffr_prepare(
    call = call,
    formula = formula,
    yind = if (yind_missing) NULL else yind,
    yind_missing = yind_missing,
    yind_expr = if (yind_missing) NULL else substitute(yind),
    data = data,
    ydata = ydata,
    algorithm = algorithm,
    method = method,
    tensortype = tensortype,
    bs_yindex = bs.yindex,
    bs_int = bs.int,
    dots = list(...)
  )
  algorithm_chr <- as.character(prep$algorithm)

  # AR(1) working correlation: the sandwich assumes working-independence
  # scores, so the default resolves to model-based intervals and an explicit
  # sandwich request errors before anything is fitted.
  if (prep$use_ar && use_sandwich) {
    if (sandwich_missing) {
      message(
        "pffr inference: model-based intervals (sandwich = FALSE) because ",
        "an AR(1) working correlation (rho) is used; the CL2 sandwich ",
        "assumes working-independence scores."
      )
      use_sandwich <- FALSE
    } else {
      pffr_check_sandwich_ar1(NULL, TRUE, rho = prep$dots$rho)
    }
  }

  meta_dims <- list(
    nobs = prep$nobs,
    nyindex = prep$nyindex,
    is_sparse = prep$is_sparse,
    missing_indices = prep$missing_indices,
    ydata = prep$ydata
  )
  # Validates a user-supplied grouping before anything is fitted.
  if (!is.null(cluster_value)) build_cluster_id(meta_dims, cluster_value)

  # Fit the model
  ncv <- NULL
  if (identical(method, "NCV")) {
    ncv <- pffr_ncv_setup(prep, cluster_value, ncv_blocks, environment())
    if (algorithm_chr == "gam") {
      ncv_setup <- ncv$setup
      prep$new_call$G <- quote(ncv_setup)
      prep$new_call$sp <- NULL # fixed sp already incorporated by gam(fit = FALSE)
    } else {
      prep$new_call <- ncv$fit_call
    }
    prep$new_call$nei <- ncv$info$nei
  }
  m <- eval(prep$new_call)
  if (algorithm_chr == "jagam") {
    m$modelfile <- prep$new_call$file
    message("JAGS/BUGS model code written to \n", m$modelfile, ",\n see ?jagam")
    return(m)
  }

  # Families without a cluster score path, gamm/gamm4 and single-curve fits
  # fall back to model-based intervals.
  if (use_sandwich) {
    reason <- pffr_cl2_unavailable(
      m,
      algorithm_chr,
      build_cluster_id(meta_dims, cluster_value)
    )
    if (!is.null(reason)) {
      msg <- paste0(
        "pffr inference: model-based intervals (sandwich = FALSE) because ",
        reason,
        "."
      )
      if (sandwich_missing) message(msg) else warning(msg, call. = FALSE)
      use_sandwich <- FALSE
    }
  }

  # Post-processing
  m_smooth <- if (algorithm_chr %in% c("gamm4", "gamm")) m$gam$smooth else
    m$smooth

  label_map_result <- pffr_build_label_map(
    new_term_strings = prep$new_term_strings,
    terms = prep$terms,
    add_f_int = prep$add_f_int,
    int_string = prep$int_string,
    m_smooth = m_smooth,
    where_specials = prep$where_specials,
    ffpc_terms = prep$ffpc_terms,
    formula_env = prep$formula_env,
    yindex_vec_name = prep$yindex_vec_name
  )
  label_map <- label_map_result$label_map
  term_map <- label_map_result$term_map

  names(m_smooth) <- sapply(m_smooth, \(x) x$label)
  if (algorithm_chr %in% c("gamm4", "gamm")) {
    m$gam$smooth <- m_smooth
  } else {
    m$smooth <- m_smooth
  }

  short_labels <- create_shortlabels(
    label_map = label_map,
    m_smooth = m_smooth,
    yind_name = prep$yind_name,
    where_specials = prep$where_specials,
    family = m$family
  )

  ret <- pffr_build_metadata(
    call = call,
    formula = formula,
    term_map = term_map,
    label_map = label_map,
    short_labels = short_labels,
    response_name = prep$response_name,
    nobs = prep$nobs,
    nyindex = prep$nyindex,
    yind_name = prep$yind_name,
    yind = prep$yind,
    where_specials = prep$where_specials,
    ff_terms = prep$ff_terms,
    ffpc_terms = prep$ffpc_terms,
    pcre_terms = prep$pcre_terms,
    missing_indices = prep$missing_indices,
    is_sparse = prep$is_sparse,
    ydata = prep$ydata,
    sandwich = if (use_sandwich) "cl2" else "none"
  )
  ret$cluster <- cluster_value
  ret$ncv <- ncv$info
  m <- pffr_attach_metadata(m, prep$algorithm, ret)

  if (!use_sandwich) {
    return(m)
  }
  # $Vp/$Vc/$Ve remain the model-based matrices mgcv produced; the CL2
  # covariance is stored separately in $pffr$Vsandwich. Only gam/bam fits
  # reach this point.
  apply_sandwich_correction(m)
}
