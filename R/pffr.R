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
#' @param sandwich logical; if \code{TRUE}, apply sandwich correction to
#'   covariance matrices to obtain robust standard errors. This can improve
#'   inference when residuals are heteroscedastic or correlated. Uses
#'   \code{\link[mgcv]{vcov.gam}} with \code{sandwich = TRUE}. Default is
#'   \code{FALSE}.
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
  bs.yindex = list(bs = "ps", k = 5, m = c(2, 1)), # only bs, k, m are propagated...
  bs.int = list(bs = "ps", k = 20, m = c(2, 1)), # only bs, k, m are propagated...
  sandwich = FALSE,
  ...
) {
  # TODO: subset args!

  call <- match.call()
  tensortype <- as.symbol(match.arg(tensortype))
  # make sure we use values for the args that were defined as close to the
  #  actual function call as possible:
  lapply(
    names(head(call, -1))[-1],
    function(nm) try(assign(nm, eval(nm, parent.frame())))
  )
  ## TODO: does this make sense? useful if pffr is called from a function that
  ## supplies args as variables that are also defined differently in GlobalEnv:
  ## this should then ensure that the args as defined in the calling function,
  ## and not in GlobalEnv get used....

  ## warn if entries in ... aren't arguments for gam/gam.fit3/jagam or gamm4/lmer
  ## check for special case of gaulss family
  dots <- list(...)
  rhoArg <- dots[["rho"]]
  useAR <- FALSE
  gaulss <- FALSE
  if (length(dots)) {
    if ("AR.start" %in% names(dots)) {
      stop(
        "Please do not supply `AR.start` directly; pffr constructs it automatically when `rho` is specified."
      )
    }
    validDots <- if (!is.na(algorithm) && algorithm == "gamm4") {
      c(names(formals(gamm4)), names(formals(lmer)))
    } else {
      c(
        names(formals(gam)),
        names(formals(bam)),
        names(formals(gam.fit3)),
        names(formals(jagam))
      )
    }
    if (!is.null(dots$family)) {
      if (
        (is.character(dots$family) && dots$family == "gaulss") |
          (is.list(dots$family) && dots$family$family == "gaulss")
      ) {
        validDots <- c(validDots, "varformula")
        gaulss <- TRUE
      }
    }

    notUsed <- names(dots)[!(names(dots) %in% validDots)]
    if (length(notUsed))
      warning(
        "Arguments <",
        paste(notUsed, collapse = ", "),
        "> supplied but not used."
      )
  }

  # Validate rho argument for AR(1) errors
  if (!is.null(rhoArg)) {
    if (!(is.numeric(rhoArg) && length(rhoArg) == 1 && !is.na(rhoArg))) {
      stop("`rho` must be a single numeric value.")
    }
    if (abs(rhoArg) >= 1) {
      stop("`rho` must have absolute value strictly less than 1.")
    }
    useAR <- abs(rhoArg) > 0
  }
  if (useAR) {
    # Validate family constraints for AR(1) errors
    familyObj <- gaussian()
    if ("family" %in% names(dots) && !is.null(dots$family)) {
      fam <- dots$family
      if (is.character(fam)) {
        if (length(fam) != 1) {
          stop("Character `family` specifications must have length 1.")
        }
        fam <- match.fun(fam)
      }
      if (is.function(fam)) {
        fam <- fam()
      }
      if (!is.list(fam) || is.null(fam$family) || is.null(fam$link)) {
        stop("Unable to interpret `family` argument when `rho` is supplied.")
      }
      familyObj <- fam
    }
    gaussianIdentity <- identical(familyObj$family, "gaussian") &&
      identical(familyObj$link, "identity")
    discreteRequested <- isTRUE(dots$discrete)
    if (!gaussianIdentity && !discreteRequested) {
      # Auto-enable discrete=TRUE for non-Gaussian families with rho
      dots$discrete <- TRUE
      message(
        "Note: Setting `discrete = TRUE` automatically because `rho` is specified with a non-Gaussian family (see ?mgcv::bam)."
      )
    }
  }

  is_sparse <- !is.null(ydata)
  if (is_sparse) {
    stopifnot(ncol(ydata) == 3)
    stopifnot(c(".obs", ".index", ".value") == colnames(ydata))
  }

  # Parse formula and classify terms using modular parser
  parsed <- parse_pffr_model_formula(formula, data, ydata)
  tf <- parsed$tf
  trmstrings <- parsed$trmstrings
  terms <- parsed$terms
  frmlenv <- parsed$frmlenv
  where.specials <- parsed$where.specials
  responsename <- parsed$responsename

  #start new formula
  newfrml <- paste(responsename, "~", sep = "")
  formula_env <- new.env()
  evalenv <- if ("data" %in% names(call)) eval.parent(call$data) else NULL

  if (is_sparse) {
    nobs <- length(unique(ydata$.obs))
    stopifnot(all(ydata$.obs %in% rownames(data)))
    # FIXME: allow for non-1:nobs .obs-formats!
    stopifnot(all(ydata$.obs %in% 1:nobs))

    #works for data-lists or matrix-valued covariates as well:
    nobs.data <- nrow(as.matrix(data[[1]]))
    stopifnot(nobs == nobs.data)
    ntotal <- nrow(ydata)

    #generate yind for estimates/predictions etc
    yind <- if (length(unique(ydata$.index)) > 100) {
      seq(min(ydata$.index), max(ydata$.index), l = 100)
    } else {
      sort(unique(ydata$.index))
    }
    nyindex <- length(yind)
  } else {
    nobs <- nrow(eval(responsename, envir = evalenv, enclos = frmlenv))
    nyindex <- ncol(eval(responsename, envir = evalenv, enclos = frmlenv))
    ntotal <- nobs * nyindex
  }

  if (missing(algorithm) || is.na(algorithm)) {
    algorithm <- ifelse(ntotal > 1e5, "bam", "gam")
  }
  if (algorithm == "bam" & missing(method)) {
    call$method <- "fREML"
  }
  algorithm <- as.symbol(algorithm)
  if (as.character(algorithm) == "bam" && !("chunk.size" %in% names(call))) {
    call$chunk.size <- 10000
    #same default as in bam
  }
  ## no te-terms possible in gamm4:
  if (as.character(algorithm) == "gamm4") {
    stopifnot(length(unlist(where.specials[c("te", "ti")])) < 1)
  }
  ## AR(1) errors only supported for bam
  if (useAR && as.character(algorithm) != "bam") {
    stop(
      "Autocorrelated errors via `rho` are currently supported only when `algorithm = \"bam\"`."
    )
  }

  if (!is_sparse) {
    #if missing, define y-index or get it from first ff/sff-term, then assign expanded versions to formula_env
    if (missing(yind)) {
      if (length(c(where.specials$ff, where.specials$sff))) {
        if (length(where.specials$ff)) {
          ffcall <- expand.call(ff, as.call(terms[where.specials$ff][1])[[1]])
        } else
          ffcall <- expand.call(sff, as.call(terms[where.specials$sff][1])[[1]])
        if (!is.null(ffcall$yind)) {
          yind <- eval(ffcall$yind, envir = evalenv, enclos = frmlenv)
          yindname <- deparse(ffcall$yind)
        } else {
          yind <- 1:nyindex
          yindname <- "yindex"
        }
      } else {
        yind <- 1:nyindex
        yindname <- "yindex"
      }
    } else {
      if (is.symbol(substitute(yind)) | is.character(yind)) {
        yindname <- deparse(substitute(yind))
        if (!is.null(data) && !is.null(data[[yindname]])) {
          yind <- data[[yindname]]
        }
      } else {
        yindname <- "yindex"
      }
      stopifnot(is.vector(yind), is.numeric(yind), length(yind) == nyindex)
    }
    #make sure it's a valid name
    if (length(yindname) > 1) yindname <- "yindex"
    # make sure yind is sorted
    stopifnot(all.equal(order(yind), 1:nyindex))

    yindvec <- rep(yind, times = nobs)
    yindex_vec_name <- as.symbol(paste(yindname, ".vec", sep = ""))
    assign(x = deparse(yindex_vec_name), value = yindvec, envir = formula_env)

    #assign response in _long_ format to formula_env
    assign(
      x = deparse(responsename),
      value = as.vector(t(eval(
        responsename,
        envir = evalenv,
        enclos = frmlenv
      ))),
      envir = formula_env
    )

    missing_indices <- if (
      any(is.na(get(as.character(responsename), formula_env)))
    ) {
      which(is.na(get(as.character(responsename), formula_env)))
    } else NULL

    # repeat which row in <data> how many times
    obs_indices <- rep(1:nobs, each = nyindex)
  } else {
    # is_sparse:
    yindname <- "yindex"
    yindvec <- ydata$.index
    yindex_vec_name <- as.symbol(paste(yindname, ".vec", sep = ""))
    assign(
      x = deparse(yindex_vec_name),
      value = ydata$.index,
      envir = formula_env
    )

    #assign response in _long_ format to formula_env
    assign(x = deparse(responsename), value = ydata$.value, envir = formula_env)

    missing_indices <- NULL

    # repeat which row in <data> how many times:
    obs_indices <- ydata$.obs
  }

  ##################################################################################
  #modify formula terms....
  newtrmstrings <- attr(tf, "term.labels")

  #if intercept, add \mu(yindex)
  if (parsed$has_intercept) {
    # Use modular intercept transformer
    int_result <- transform_intercept_term(yindex_vec_name, bs.int, yindname)
    intstring <- int_result$term_string

    newfrml <- paste(newfrml, intstring, sep = " ")
    addFint <- TRUE
  } else {
    newfrml <- paste(newfrml, "0", sep = "")
    addFint <- FALSE
    intstring <- NULL
  }

  #transform: c(foo) --> foo using modular transformer
  if (length(where.specials$c)) {
    newtrmstrings[where.specials$c] <- sapply(
      trmstrings[where.specials$c],
      transform_c_term
    )
  }

  #prep function-on-function-terms
  if (length(c(where.specials$ff, where.specials$sff))) {
    ffterms <- lapply(
      terms[c(where.specials$ff, where.specials$sff)],
      function(x) {
        eval(x, envir = evalenv, enclos = frmlenv)
      }
    )

    newtrmstrings[c(where.specials$ff, where.specials$sff)] <- sapply(
      ffterms,
      function(x) {
        safeDeparse(x$call)
      }
    )

    #apply limits function and assign stacked data to formula_env
    makeff <- function(x) {
      tmat <- matrix(yindvec, nrow = length(yindvec), ncol = length(x$xind))
      smat <- matrix(
        x$xind,
        nrow = length(yindvec),
        ncol = length(x$xind),
        byrow = TRUE
      )
      if (!is.null(x[["LX"]])) {
        # for ff: stack weights * covariate
        LStacked <- x$LX[obs_indices, ]
      } else {
        # for sff: stack weights, X separately
        LStacked <- x$L[obs_indices, ]
        XStacked <- x$X[obs_indices, ]
      }
      if (!is.null(x$limits)) {
        # find int-limits and set weights to 0 outside
        use <- x$limits(smat, tmat)
        LStacked <- LStacked * use

        # find indices for row-wise int-range & maximal width
        windows <- compute_integration_windows(use)
        max_width <- max(windows[, 3])

        # reduce size of matrix-covariates if possible
        if (max_width < ncol(smat)) {
          eff_windows <- expand_windows_to_maxwidth(windows, ncol(smat))
          smat <- shift_and_shorten_matrix(smat, eff_windows)
          tmat <- shift_and_shorten_matrix(tmat, eff_windows)
          LStacked <- shift_and_shorten_matrix(LStacked, eff_windows)
          if (is.null(x$LX)) {
            # sff
            XStacked <- shift_and_shorten_matrix(XStacked, eff_windows)
          }
        }
      }
      assign(x = x$yindname, value = tmat, envir = formula_env)
      assign(x = x$xindname, value = smat, envir = formula_env)

      assign(x = x$LXname, value = LStacked, envir = formula_env)
      if (is.null(x[["LX"]])) {
        # sff
        assign(x = x$xname, value = XStacked, envir = formula_env)
      }
      invisible(NULL)
    }
    lapply(ffterms, makeff)
  } else ffterms <- NULL

  if (length(where.specials$ffpc)) {
    ##TODO for sparse
    ffpcterms <- lapply(terms[where.specials$ffpc], function(x) {
      eval(x, envir = evalenv, enclos = frmlenv)
    })

    lapply(ffpcterms, function(trm) {
      lapply(colnames(trm$data), function(nm) {
        assign(x = nm, value = trm$data[obs_indices, nm], envir = formula_env)
        invisible(NULL)
      })
      invisible(NULL)
    })

    getFfpcFormula <- function(trm) {
      frmls <- lapply(colnames(trm$data), function(pc) {
        arglist <- c(
          name = "s",
          x = as.symbol(yindex_vec_name),
          by = as.symbol(pc),
          id = trm$id,
          trm$splinepars
        )
        call <- do.call("call", arglist, envir = formula_env)
        call$x <- as.symbol(yindex_vec_name)
        call$by <- as.symbol(pc)
        safeDeparse(call)
      })
      return(paste(unlist(frmls), collapse = " + "))
    }
    newtrmstrings[where.specials$ffpc] <- sapply(ffpcterms, getFfpcFormula)

    ffpcterms <- lapply(ffpcterms, function(x) x[names(x) != "data"])
  } else ffpcterms <- NULL

  #prep PC-based random effects
  if (length(where.specials$pcre)) {
    pcreterms <- lapply(terms[where.specials$pcre], function(x) {
      eval(x, envir = evalenv, enclos = frmlenv)
    })
    #assign newly created data to formula_env
    lapply(pcreterms, function(trm) {
      if (!is_sparse && all(trm$yind == yind)) {
        lapply(colnames(trm$efunctions), function(nm) {
          assign(
            x = nm,
            value = trm$efunctions[rep(1:nyindex, times = nobs), nm],
            envir = formula_env
          )
          invisible(NULL)
        })
      } else {
        # don't ever extrapolate eigenfunctions:
        stopifnot(min(trm$yind) <= min(yind))
        stopifnot(max(trm$yind) >= max(yind))

        # interpolate given eigenfunctions to observed index values:
        lapply(colnames(trm$efunctions), function(nm) {
          tmp <- approx(
            x = trm$yind,
            y = trm$efunctions[, nm],
            xout = yindvec,
            method = "linear"
          )$y
          assign(x = nm, value = tmp, envir = formula_env)
          invisible(NULL)
        })
      }
      assign(x = trm$idname, value = trm$id[obs_indices], envir = formula_env)
      invisible(NULL)
    })

    newtrmstrings[where.specials$pcre] <- sapply(pcreterms, function(x) {
      safeDeparse(x$call)
    })
  } else pcreterms <- NULL

  #transform: s(x, ...), te(x, z,...), t2(x, z, ...) --> <ti|t2>(x, <z,> yindex, ..., <bs.yindex>)
  # Using modular transformer
  if (length(c(where.specials$s, where.specials$te, where.specials$t2))) {
    newtrmstrings[c(where.specials$s, where.specials$te, where.specials$t2)] <-
      sapply(
        terms[c(where.specials$s, where.specials$te, where.specials$t2)],
        function(x)
          transform_smooth_term(
            x,
            yindex_vec_name,
            bs.yindex,
            tensortype,
            algorithm
          )
      )
  }

  #transform: x --> s(YINDEX, by=x) using modular transformer
  if (length(where.specials$par)) {
    newtrmstrings[where.specials$par] <- sapply(
      terms[where.specials$par],
      function(x) transform_par_term(x, yindex_vec_name, bs.yindex)
    )
  }

  #... & assign expanded/additional variables to formula_env
  where.specials$notff <- c(
    where.specials$c,
    where.specials$par,
    where.specials$s,
    where.specials$te,
    where.specials$t2
  )
  if (length(where.specials$notff)) {
    # evalenv below used to be list2env(eval.parent(call$data)), frmlenv),
    # but that assigned everything in <data> to the global workspace if frmlenv was the global
    # workspace.
    evalenv <- if ("data" %in% names(call)) {
      list2env(eval.parent(call$data))
    } else frmlenv
    lapply(terms[where.specials$notff], function(x) {
      #nms <- all.vars(x)
      isC <- safeDeparse(x) %in% sapply(terms[where.specials$c], safeDeparse)
      if (isC) {
        # drop c()
        # FIXME: FUGLY!
        x <- formula(paste(
          "~",
          gsub("\\)$", "", gsub("^c\\(", "", deparse(x)))
        ))[[2]]
      }
      ## remove names in xt, k, bs,  information (such as variable names for MRF penalties etc)
      nms <- if (!is.null(names(x))) {
        all.vars(x[names(x) %in% c("", "by")])
      } else all.vars(x)

      sapply(nms, function(nm) {
        var <- get(nm, envir = evalenv)
        if (is.matrix(var)) {
          stopifnot(!is_sparse || ncol(var) == nyindex)
          assign(x = nm, value = as.vector(t(var)), envir = formula_env)
        } else {
          stopifnot(length(var) == nobs)
          assign(x = nm, value = var[obs_indices], envir = formula_env)
        }
        invisible(NULL)
      })
      invisible(NULL)
    })
  }

  # Build mgcv formula using modular function
  newfrml <- build_mgcv_formula(
    responsename = responsename,
    intercept_string = if (addFint) intstring else NULL,
    term_strings = newtrmstrings,
    has_intercept = addFint,
    formula_env = formula_env
  )

  # variance formula for gaulss
  if (gaulss) {
    if (is.null(dots$varformula)) {
      dots$varformula <- formula(paste(
        "~",
        safeDeparse(
          as.call(c(as.name("s"), x = as.symbol(yindex_vec_name), bs.int))
        )
      ))
    }
    environment(dots$varformula) <- formula_env
    newfrml <- list(newfrml, dots$varformula)
  }

  # Build AR.start indicator if AR(1) errors requested
  if (useAR) {
    resp_long <- get(as.character(responsename), envir = formula_env)
    obs_ids <- obs_indices
    valid_idx <- which(!is.na(resp_long))
    if (!length(valid_idx)) {
      stop("Cannot build AR.start because all responses are missing.")
    }
    start_idx <- rep(FALSE, length(resp_long))
    splits <- split(valid_idx, obs_ids[valid_idx])
    start_positions <- as.integer(vapply(
      splits,
      function(idx) idx[1],
      integer(1)
    ))
    start_idx[start_positions] <- TRUE
    assign("AR.start", value = start_idx, envir = formula_env)
  }

  # Build mgcv data using modular function
  pffrdata <- build_mgcv_data(formula_env)

  newcall <- expand.call(pffr, call)
  newcall$yind <- newcall$tensortype <- newcall$bs.int <-
    newcall$bs.yindex <- newcall$algorithm <- newcall$ydata <- NULL
  newcall$sandwich <- NULL
  newcall$formula <- newfrml
  newcall$data <- quote(pffrdata)
  newcall[[1]] <- algorithm
  if (useAR) {
    newcall$AR.start <- quote(pffrdata$AR.start)
    # Ensure discrete=TRUE is passed to bam if it was auto-set
    if (isTRUE(dots$discrete)) {
      newcall$discrete <- TRUE
    }
  }
  # make sure ...-args are taken from ..., not GlobalEnv:
  dotargs <- names(newcall)[names(newcall) %in% names(dots)]
  newcall[dotargs] <- dots[dotargs]
  if ("subset" %in% dotargs) {
    stop("<subset>-argument is not supported.")
  }
  if ("weights" %in% dotargs) {
    wtsdone <- FALSE
    if (length(dots$weights) == nobs) {
      newcall$weights <- dots$weights[obs_indices]
      wtsdone <- TRUE
    }
    if (
      !is.null(dim(dots$weights)) &&
        all(dim(dots$weights) == c(nobs, nyindex))
    ) {
      newcall$weights <- as.vector(t(dots$weights))
      wtsdone <- TRUE
    }
    if (!wtsdone) {
      stop(
        "weights have to be supplied as a vector with length=rows(data) or
               a matrix with the same dimensions as the response."
      )
    }
  }
  if ("offset" %in% dotargs) {
    ofstdone <- FALSE
    if (length(dots$offset) == nobs) {
      newcall$offset <- dots$offset[obs_indices]
      ofstdone <- TRUE
    }
    if (
      !is.null(dim(dots$offset)) &&
        all(dim(dots$offset) == c(nobs, nyindex))
    ) {
      newcall$offset <- as.vector(t(dots$offset))
      ofstdone <- TRUE
    }
    if (!ofstdone) {
      stop(
        "offsets have to be supplied as a vector with length=rows(data) or
             a matrix with the same dimensions as the response."
      )
    }
  }
  if (as.character(algorithm) == "jagam") {
    newcall <- newcall[names(newcall) %in% c("", names(formals(jagam)))]
    if (is.null(newcall$file)) {
      newcall$file <- tempfile(
        "pffr2jagam",
        tmpdir = getwd(),
        fileext = ".jags"
      )
    }
  }
  # call algorithm to estimate model
  m <- eval(newcall)
  if (as.character(algorithm) == "jagam") {
    m$modelfile <- newcall$file
    message("JAGS/BUGS model code written to \n", m$modelfile, ",\n see ?jagam")
    return(m)
  }

  m.smooth <- if (as.character(algorithm) %in% c("gamm4", "gamm")) {
    m$gam$smooth
  } else m$smooth

  #return some more info s.t. custom predict/plot/summary will work
  trmmap <- newtrmstrings
  names(trmmap) <- names(terms)
  if (addFint) trmmap <- c(trmmap, intstring)

  # map labels to terms --
  # ffpc are associated with multiple smooths
  # parametric are associated with multiple smooths if covariate is a factor
  labelmap <- as.list(trmmap)
  lbls <- sapply(m.smooth, function(x) x$label)
  if (length(c(where.specials$par, where.specials$ffpc))) {
    if (length(where.specials$par)) {
      for (w in where.specials$par) {
        # only combine if <by>-variable is a factor!
        if (is.factor(get(names(labelmap)[w], envir = formula_env))) {
          labelmap[[w]] <- {
            #covariates for parametric terms become by-variables:
            where <- sapply(m.smooth, function(x) x$by) == names(labelmap)[w]
            sapply(m.smooth[where], function(x) x$label)
          }
        } else {
          labelmap[[w]] <- paste0(
            "s(",
            yindex_vec_name,
            "):",
            names(labelmap)[w]
          )
        }
      }
    }
    if (length(where.specials$ffpc)) {
      ind <- 1
      for (w in where.specials$ffpc) {
        labelmap[[w]] <- {
          #PCs for X become by-variables:
          where <- sapply(m.smooth, function(x) x$id) == ffpcterms[[ind]]$id
          sapply(m.smooth[where], function(x) x$label)
        }
        ind <- ind + 1
      }
    }
    labelmap[-c(where.specials$par, where.specials$ffpc)] <- lbls[pmatch(
      sapply(
        labelmap[-c(where.specials$par, where.specials$ffpc)],
        function(x) {
          ## FUGLY: check whether x is a function call of some sort
          ##  or simply a variable name.
          if (length(parse(text = x)[[1]]) != 1) {
            tmp <- eval(parse(text = x))
            return(tmp$label)
          } else {
            return(x)
          }
        }
      ),
      lbls
    )]
  } else {
    labelmap[1:length(labelmap)] <- lbls[pmatch(
      sapply(labelmap[1:length(labelmap)], function(x) {
        ## FUGLY: check whether x is a function call of some sort
        ##  or simply a variable name.
        if (length(parse(text = x)[[1]]) != 1) {
          tmp <- eval(parse(text = x))
          return(tmp$label)
        } else {
          return(x)
        }
      }),
      lbls
    )]
  }
  # check whether any parametric terms were left out & add them
  nalbls <- sapply(labelmap, function(x) {
    any(is.null(x)) | any(is.na(x[!is.null(x)]))
  })
  if (any(nalbls)) {
    labelmap[nalbls] <- trmmap[nalbls]
  }

  names(m.smooth) <- lbls
  if (as.character(algorithm) %in% c("gamm4", "gamm")) {
    m$gam$smooth <- m.smooth
  } else {
    m$smooth <- m.smooth
  }

  # Create shortlabels mapping: mgcv label -> human-readable pffr label
  shortlabels <- create_shortlabels(
    labelmap = labelmap,
    m.smooth = m.smooth,
    yindname = yindname,
    where.specials = where.specials,
    family = m$family
  )

  ret <- list(
    call = call,
    formula = formula,
    termmap = trmmap,
    labelmap = labelmap,
    shortlabels = shortlabels,
    responsename = responsename,
    nobs = nobs,
    nyindex = nyindex,
    yindname = yindname,
    yind = yind,
    where = where.specials,
    ff = ffterms,
    ffpc = ffpcterms,
    pcreterms = pcreterms,
    missing_indices = missing_indices,
    is_sparse = is_sparse,
    ydata = ydata,
    sandwich = sandwich
  )

  if (as.character(algorithm) %in% c("gamm4", "gamm")) {
    m$gam$pffr <- ret
    class(m$gam) <- c("pffr", class(m$gam))
  } else {
    m$pffr <- ret
    class(m) <- c("pffr", class(m))
  }

  # Apply sandwich correction if requested
  if (sandwich) {
    gam_obj <- if (as.character(algorithm) %in% c("gamm4", "gamm")) m$gam else m
    # Strip pffr class for vcov calls to avoid predict.pffr warnings
    gam_obj_stripped <- gam_obj
    class(gam_obj_stripped) <- setdiff(class(gam_obj_stripped), "pffr")
    # Overwrite both the frequentist and the Bayesian covariance matrix
    gam_obj$Vp <- gam_obj$Vc <- stats::vcov(gam_obj_stripped, sandwich = TRUE)
    gam_obj$Ve <- stats::vcov(gam_obj_stripped, sandwich = TRUE, freq = TRUE)
    if (as.character(algorithm) %in% c("gamm4", "gamm")) {
      m$gam <- gam_obj
    } else {
      m <- gam_obj
    }
  }

  return(m)
} # end pffr()
