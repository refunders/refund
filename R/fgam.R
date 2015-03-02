#' Functional Generalized Additive Models
#'
#' Implements functional generalized additive models for functional and scalar covariates and scalar responses.
#' Additionally implements functional linear models.  This function is a wrapper for mgcv's \code{\link{gam}}
#' and its siblings to fit models of the general form
#' \deqn{g(E(Y_i)) = \beta_0 + \int_{T_1} F(X_{i1},t)dt+ \int_{T_2} \beta(t)X_{i2}dt + f(z_{i1}) + f(z_{i2}, z_{i3}) + \ldots}
#' with a scalar (but not necessarily continuous) response Y, and link function g
#' @param formula a formula with special terms as for gam, with additional special terms
#' \code{\link{af}}() and \code{\link{lf}}().
#' @param fitter the name of the function used to estimate the model. Defaults to \code{\link{gam}}
#' if the matrix of functional responses has less than 2e5 data points and to
#' \code{\link{bam}} if not. "gamm" (see \code{\link{gamm}}) and "gamm4"
#' (see \code{\link{gamm4}}) are valid options as well.
#' @param ... additional arguments that are valid for \code{\link{gam}} or \code{\link{bam}}; for example,
#' specify a \code{gamma} > 1 to increase amount of smoothing when using GCV to choose smoothing
#' parameters or \code{method="REML"} to change to REML for estimation of smoothing parameters
#' (default is GCV).
#' @section Warning:
#' Binomial responses should be specified as a numeric vector rather than as a matrix or a factor.
#' @return a fitted fgam-object, which is a \code{\link{gam}}-object with some additional information
#' in a fgam-entry. If fitter is "gamm" or "gamm4", only the $gam part of the
#' returned list is modified in this way.
#' @references McLean, M. W., Hooker, G., Staicu, A.-M., Scheipl, F., and Ruppert, D. (2014). Functional
#' generalized additive models. \emph{Journal of Computational and Graphical Statistics}, \bold{23 (1)},
#' pp. 249-269.  Available at \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3982924}.
#' @author Mathew W. McLean \email{mathew.w.mclean@@gmail.com} and Fabian Scheipl
#' @seealso \code{\link{af}}, \code{\link{lf}}, \code{\link{predict.fgam}}, \code{\link{vis.fgam}}
#' @importFrom mgcv gam gam.fit gamm4 bam s te t2
#' @importFrom gamm4 gamm4
#' @importFrom nlme lme4
#' @importFrom stats terms.formula
#' @examples
#' data(DTI)
#' ## only consider first visit and cases (no PASAT scores for controls)
#' y <- DTI$pasat[DTI$visit==1 & DTI$case==1]
#' X <- DTI$cca[DTI$visit==1 & DTI$case==1, ]
#' X_2 <- DTI$rcst[DTI$visit==1 & DTI$case==1, ]
#'
#' ## remove samples containing missing data
#' ind <- rowSums(is.na(X)) > 0
#' ind2 <- rowSums(is.na(X_2)) > 0
#'
#' y <- y[!(ind | ind2)]
#' X <- X[!(ind | ind2), ]
#' X_2 <- X_2[!(ind | ind2), ]
#'
#' N <- length(y)
#'
#' ## fit fgam using FA measurements along corpus callosum
#' ## as functional predictor with PASAT as response
#' ## using 8 cubic B-splines for marginal bases with third
#' ## order marginal difference penalties
#' ## specifying gamma > 1 enforces more smoothing when using
#' ## GCV to choose smoothing parameters
#' fit <- fgam(y ~ af(X, splinepars = list(k = c(8, 8), m = list(c(2, 3), c(2, 3)))), gamma = 1.2)
#' vis.fgam(fit)
#'
#'
#' ## fgam term for the cca measurements plus an flm term for the rcst measurements
#' ## leave out 10 samples for prediction
#' test <- sample(N, 10)
#' fit <- fgam(y ~ af(X, splinepars = list(k = c(7, 7), m = list(c(2, 2), c(2, 2)))) +
#'   lf(X_2, splinepars = list(k=7, m = c(2, 2))), subset=(1:N)[-test])
#' ## plot the fits
#' plot(fit)
#' ## vis.fgam(fit, af.term = "X")
#' vis.fgam(fit, af.term="X", xval=.6)
#' ## predict the ten left outs samples
#' pred <- predict(fit, newdata = list(X=X[test, ], X_2 = X_2[test, ]), type='response',
#'                 PredOutOfRange = TRUE)
#' sqrt(mean((y[test] - pred)^2))

#' ## Try to predict the binary response disease status (case or control)
#' ##   using the quantile transformed measurements from the rcst tract
#' ##   with a smooth component for a scalar covariate that is pure noise
#' y <- DTI$case[DTI$visit==1]
#' X <- DTI$cca[DTI$visit==1, ]
#' X_2 <- DTI$rcst[DTI$visit==1, ]
#'
#' ind <- rowSums(is.na(X)) > 0
#' ind2 <- rowSums(is.na(X_2)) > 0
#'
#' y <- y[!(ind | ind2)]
#' X <- X[!(ind | ind2), ]
#' X_2 <- X_2[!(ind | ind2), ]
#' z1 <- rnorm(length(y))
#'
#' ## select=TRUE allows terms to be zeroed out of model completely
#' fit <- fgam(y ~ s(z1, k = 10) + af(X_2, splinepars = list(k=c(7,7), m = list(c(2, 1), c(2, 1))),
#'             Qtransform=TRUE), family=binomial(), select=TRUE)
#' plot(fit)
#' vis.fgam(fit,af.term='X_2',plot.type='contour')
fgam <- function(formula, fitter=NA, ...){

  call <- match.call()
  dots <- list(...)
  if (length(dots)) {
    validDots <- if (!is.na(fitter) && fitter == "gamm4") {
      c(names(formals(gamm4)), names(formals(lmer)))
    }
    else {
      c(names(formals(gam)), names(formals(gam.fit)))
    }
    notUsed <- names(dots)[!(names(dots) %in% validDots)]
    if (length(notUsed))
      warning("Arguments <", paste(notUsed, collapse = ", "),
              "> supplied but not used.")
  }
  tf <- terms.formula(formula, specials = c("s", "te", "t2", "lf", "af"))
  trmstrings <- attr(tf, "term.labels")
  terms <- sapply(trmstrings, function(trm) as.call(parse(text = trm))[[1]],
                  simplify = FALSE)
  frmlenv <- environment(formula)
  where.af <- attr(tf, "specials")$af - 1
  where.lf <- attr(tf, "specials")$lf - 1
  where.s <- attr(tf, "specials")$s - 1
  where.te <- attr(tf, "specials")$te - 1
  where.t2 <- attr(tf, "specials")$t2 - 1


  if (length(trmstrings)) {
    where.par <- which(!(1:length(trmstrings) %in% c(where.af, where.lf, where.s,where.te, where.t2)))
  }else where.par <- numeric(0)

  responsename <- attr(tf, "variables")[2][[1]]
  newfrml <- paste(responsename, "~", sep = "")
  newfrmlenv <- new.env()
  evalenv <- if ("data" %in% names(call))
    eval(call$data)
  else NULL
  nobs <- length(eval(responsename, envir = evalenv, enclos = frmlenv))

  if (missing(fitter) || is.na(fitter)) {
    fitter <- ifelse(nobs > 1e+05, "bam", "gam")
  }

  fitter <- as.symbol(fitter)
  if (as.character(fitter) == "bam" && !("chunk.size" %in%
    names(call))) {
    call$chunk.size <- max(nobs/5, 10000)
  }
  if (as.character(fitter) == "gamm4")
    stopifnot(length(where.te) < 1)

  assign(x = deparse(responsename), value = as.vector(t(eval(responsename,
                                      envir = evalenv, enclos = frmlenv))), envir = newfrmlenv)

  newtrmstrings <- attr(tf, "term.labels")
  if (!attr(tf, "intercept")) {
    newfrml <- paste(newfrml, "0", sep = "")
  }

  if (length(c(where.af, where.lf))) {
    fterms <- lapply(terms[c(where.af, where.lf)], function(x) {
      eval(x, envir = evalenv, enclos = frmlenv)
    })
    newtrmstrings[c(where.af, where.lf)] <- sapply(fterms,
                                                    function(x) {
                                                      safeDeparse(x$call)
                                                    })
    lapply(fterms, function(x) {
      lapply(names(x$data), function(nm) {
        assign(x = nm, value = x$data[[nm]], envir = newfrmlenv)
        invisible(NULL)
      })
      invisible(NULL)
    })
    fterms <- lapply(fterms, function(x) x[names(x) !=
      "data"])
  }
  else fterms <- NULL

  where.notf <- c(where.par,where.s,where.te,where.t2)
  if (length(where.notf)) {
    if ("data" %in% names(call))
      frmlenv <- list2env(eval(call$data), frmlenv)
    lapply(terms[where.notf], function(x) {

      nms <- if (!is.null(names(x))) {
        all.vars(x[names(x) == ""])
      }
      else all.vars(x)
      sapply(nms, function(nm) {
        stopifnot(length(get(nm, envir = frmlenv)) ==
          nobs)
        assign(x = nm, value = get(nm, envir = frmlenv), envir = newfrmlenv)
        invisible(NULL)
      })
      invisible(NULL)
    })
  }

  newfrml <- formula(paste(c(newfrml, newtrmstrings), collapse = "+"))
  environment(newfrml) <- newfrmlenv
  fgamdata <- list2df(as.list(newfrmlenv))
  datameans <- sapply(as.list(newfrmlenv),mean)
  newcall <- expand.call(fgam, call)
  newcall$fitter <- type <- newcall$bs.int <- newcall$bs.yindex <- newcall$fitter <- NULL
  newcall$formula <- newfrml
  newcall$data <- quote(fgamdata)
  newcall$fitter <- NULL
  newcall[[1]] <- fitter


  res <- eval(newcall)

  res.smooth <- if (as.character(fitter) %in% c("gamm4", "gamm")) {
    res$gam$smooth
  }else res$smooth

  trmmap <- newtrmstrings
  names(trmmap) <- names(terms)
  labelmap <- as.list(trmmap)
  names(trmmap) <- names(terms)


  lbls <- sapply(res.smooth, function(x) x$label)
  if (length(where.par)) {
    for (w in where.par) labelmap[[w]] <- {
      where <- sapply(res.smooth, function(x) x$by) ==
        names(labelmap)[w]
      sapply(res.smooth[where], function(x) x$label)
    }
  labelmap[-c(where.par)] <- lbls[pmatch(sapply(labelmap[-c(where.par)], function(x) {
                                                                tmp <- eval(parse(text = x))
                                                                if (is.list(tmp)) {
                                                                  return(tmp$label)
                                                                }else {
                                                                  return(x)
                                                                }
                                                              }), lbls)]
  }else{
    labelmap[1:length(labelmap)] <- lbls[pmatch(sapply(labelmap,
                                                       function(x) {
                                                         tmp <- eval(parse(text = x))
                                                         if (is.list(tmp)) {
                                                           return(tmp$label)
                                                         }
                                                         else {
                                                           return(x)
                                                         }
                                                       }), lbls)]
  }
  if (any(nalbls <- sapply(labelmap, function(x) any(is.na(x))))) {
    labelmap[nalbls] <- trmmap[nalbls]
  }
  names(res.smooth) <- lbls
  if (as.character(fitter) %in% c("gamm4", "gamm")) {
    res$gam$smooth <- res.smooth
  }else {
    res$smooth <- res.smooth
  }
  ret <- list(formula = formula, termmap = trmmap, labelmap = labelmap,
              responsename = responsename, nobs = nobs,
              where = list(where.af = where.af, where.lf = where.lf,
                              where.s = where.s, where.te = where.te, where.t2 = where.t2,
                               where.par = where.par), datameans=datameans,ft = fterms)
  if (as.character(fitter) %in% c("gamm4", "gamm")) {
    res$gam$fgam <- ret
    class(res$gam) <- c("fgam", class(res$gam))
  }
  else {
    res$fgam <- ret
    class(res) <- c("fgam", class(res))
  }

	return(res)
}
