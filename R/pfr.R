#' Penalized Functional Regression
#'
#' Implements various approaches to penalized scalar-on-function regression.
#' These techniques include
#' Penalized Functional Regression (Goldsmith et al., 2011),
#' Longitudinal Penalized Functional Regression (Goldsmith, et al., 2012),
#' Functional Principal Component Regression (Reiss and Ogden, 2007),
#' Partially Empirical Eigenvectors for Regression (Randolph et al., 2012),
#' Functional Generalized Additive Models (McLean et al., 2013),
#' and
#' Variable-Domain Functional Regression (Gellar et al., 2014).
#' This function is a wrapper for mgcv's \code{\link{gam}} and its siblings
#' to fit models with a scalar (but not necessarily continuous) response.
#'
#' @param formula a formula that could contain any of the following special terms:
#'   \code{\link{lf}()}, \code{\link{af}()}, \code{\link{lf.vd}()},
#'   \code{\link{peer}()},  \code{\link{fpc}()},
#'   or \code{\link{re}()}; also \code{mgcv}'s \code{\link{s}()},
#'   \code{\link{te}()}, or \code{\link{t2}()}.
#' @param fitter the name of the function used to estimate the model. Defaults
#'   to \code{\link{gam}} if the matrix of functional responses has less than 2e5
#'   data points and to \code{\link{bam}} if not. "gamm" (see \code{\link{gamm}})
#'   and "gamm4" (see \code{\link{gamm4}}) are valid options as well.
#' @param method The smoothing parameter estimation method. Default is
#'   \code{"REML"}. For options, see \code{\link{gam}}.
#' @param ... additional arguments that are valid for \code{\link{gam}} or
#'   \code{\link{bam}}. These include \code{data} and \code{family} to specify
#'   the input data and outcome family, as well as many options to control the
#'   estimation.
#'
#' @section Warning:
#'   Binomial responses should be specified as a numeric vector rather than as a
#'   matrix or a factor.
#' @return
#'   A fitted pfr-object, which is a \code{\link{gam}}-object with some
#'   additional information in a \code{$pfr}-element. If fitter is \code{"gamm"}
#'   or \code{"gamm4"}, only the \code{$gam} part of the returned list is
#'   modified in this way.
#'
#' @references
#' Goldsmith, J., Bobb, J., Crainiceanu, C., Caffo, B., and Reich, D. (2011).
#' Penalized functional regression. \emph{Journal of Computational and Graphical
#' Statistics}, 20(4), 830-851.
#'
#' Goldsmith, J., Crainiceanu, C., Caffo, B., and Reich, D. (2012). Longitudinal
#' penalized functional regression for cognitive outcomes on neuronal tract
#' measurements. \emph{Journal of the Royal Statistical Society: Series C},
#' 61(3), 453-469.
#'
#' Reiss, P. T., and Ogden, R. T. (2007). Functional principal component
#' regression and functional partial least squares. \emph{Journal of the
#' American Statistical Association}, 102, 984-996.
#'
#' Randolph, T. W., Harezlak, J, and Feng, Z. (2012). Structured penalties for
#' functional linear models - partially empirical eigenvectors for regression.
#' \emph{Electronic Journal of Statistics}, 6, 323-353.
#'
#' McLean, M. W., Hooker, G., Staicu, A.-M., Scheipl, F., and
#' Ruppert, D. (2014). Functional generalized additive models. \emph{Journal of
#' Computational and Graphical Statistics}, \bold{23 (1)}, pp. 249-269.
#' Available at \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3982924}.
#'
#' Gellar, J. E., Colantuoni, E., Needham, D. M., and Crainiceanu, C. M. (2014).
#' Variable-Domain Functional Regression for Modeling ICU Data. Journal of the
#' American Statistical Association, 109(508): 1425-1439.
#'
#' @author Jonathan Gellar \email{JGellar@@mathematica-mpr.com}, Mathew W. McLean,
#' Jeff Goldsmith, and Fabian Scheipl
#' @seealso \code{\link{af}}, \code{\link{lf}}, \code{\link{lf.vd}},
#'   \code{\link{fpc}}, \code{\link{peer}}, \code{\link{re}}.
#' @importFrom mgcv gam gam.fit bam s te t2
#' @importFrom gamm4 gamm4
#' @importFrom lme4 lmer
#' @importFrom stats terms.formula
#' @export
#'
#' @examples
#' # See lf(), lf.vd(), af(), fpc(), and peer() for additional examples
#'
#' data(DTI)
#' DTI1 <- DTI[DTI$visit==1 & complete.cases(DTI),]
#' par(mfrow=c(1,2))
#'
#' # Fit model with linear functional term for CCA
#' fit.lf <- pfr(pasat ~ lf(cca, k=30, bs="ps"), data=DTI1)
#' plot(fit.lf, ylab=expression(paste(beta(t))), xlab="t")
#' \dontrun{
#' # Alternative way to plot
#' bhat.lf <- coef(fit.lf, n=101)
#' bhat.lf$upper <- bhat.lf$value + 1.96*bhat.lf$se
#' bhat.lf$lower <- bhat.lf$value - 1.96*bhat.lf$se
#' matplot(bhat.lf$cca.argvals, bhat.lf[,c("value", "upper", "lower")],
#'         type="l", lty=c(1,2,2), col=1,
#'         ylab=expression(paste(beta(t))), xlab="t")
#'
#' # Fit model with additive functional term for CCA, using tensor product basis
#' fit.af <- pfr(pasat ~ af(cca, Qtransform=TRUE, k=c(7,7)), data=DTI1)
#' plot(fit.af, scheme=2, xlab="t", ylab="cca(t)", main="Tensor Product")
#' plot(fit.af, scheme=2, Qtransform=TRUE,
#'      xlab="t", ylab="cca(t)", main="Tensor Product")
#'
#' # Change basistype to thin-plate regression splines
#' fit.af.s <- pfr(pasat ~ af(cca, basistype="s", Qtransform=TRUE, k=50),
#'                 data=DTI1)
#' plot(fit.af.s, scheme=2, xlab="t", ylab="cca(t)", main="TPRS", rug=FALSE)
#' plot(fit.af.s, scheme=2, Qtransform=TRUE,
#'      xlab="t", ylab="cca(t)", main="TPRS", rug=FALSE)
#'
#' # Visualize bivariate function at various values of x
#' par(mfrow=c(2,2))
#' vis.pfr(fit.af, xval=.2)
#' vis.pfr(fit.af, xval=.4)
#' vis.pfr(fit.af, xval=.6)
#' vis.pfr(fit.af, xval=.8)
#'
#' # Include random intercept for subject
#' DTI.re <- DTI[complete.cases(DTI$cca),]
#' DTI.re$ID <- factor(DTI.re$ID)
#' fit.re <- pfr(pasat ~ lf(cca, k=30) + re(ID), data=DTI.re)
#' coef.re <- coef(fit.re)
#' par(mfrow=c(1,2))
#' plot(fit.re)
#'
#' # FPCR_R Model
#' fit.fpc <- pfr(pasat ~ fpc(cca), data=DTI.re)
#' plot(fit.fpc)
#'
#' # PEER Model with second order difference penalty
#' DTI.use <- DTI[DTI$case==1,]
#' DTI.use <- DTI.use[complete.cases(DTI.use$cca),]
#' fit.peer <- pfr(pasat ~ peer(cca, argvals=seq(0,1,length=93),
#'                              integration="riemann", pentype="D"), data=DTI.use)
#' plot(fit.peer)
#' }

pfr <- function(formula=NULL, fitter=NA, method="REML", ...){

  if (class(formula) != "formula") {
    warning(paste0("The interface for pfr() has changed to using a formula ",
                   "argument, with linear functional terms specified by lf(). ",
                   "See ?pfr for details. The old interface will be deprecated ",
                   "in the next refund release."))
    # Call pfr_old()
    call <- sys.call()
    call[[1]] <- as.symbol("pfr_old")
    fit <- eval(call)
    return(fit)
  }

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

  # Set up terms
  tf <- terms.formula(formula, specials = c("s", "te", "t2", "lf", "af",
                                            "lf.vd", "re", "peer", "fpc"))
  trmstrings <- attr(tf, "term.labels")
  terms <- sapply(trmstrings, function(trm) as.call(parse(text = trm))[[1]],
                  simplify = FALSE)
  frmlenv <- environment(formula)
  specials <- attr(tf, "specials")
  where.af <- specials$af - 1
  where.lf <- specials$lf - 1
  where.pr <- specials$peer - 1
  where.fp <- specials$fpc  - 1
  where.s  <- specials$s  - 1
  where.te <- specials$te - 1
  where.t2 <- specials$t2 - 1
  where.re <- specials$re - 1
  where.lf.vd <- specials$lf.vd - 1
  where.all <- c(where.af, where.lf, where.s, where.te, where.t2, where.re,
                 where.lf.vd, where.pr, where.fp)

  if (length(trmstrings)) {
    where.par <- which(!(1:length(trmstrings) %in% where.all))
  } else where.par <- numeric(0)

  # Set up new formula and response
  responsename <- attr(tf, "variables")[2][[1]]
  newfrml <- paste(responsename, "~", sep = "")
  newfrmlenv <- new.env()
  evalenv <- if ("data" %in% names(call))
    eval.parent(call$data)
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

  assign(x = deparse(responsename),
         value = as.vector(t(eval(responsename, envir = evalenv,
                                  enclos = frmlenv))),
         envir = newfrmlenv)

  newtrmstrings <- attr(tf, "term.labels")
  if (!attr(tf, "intercept")) {
    newfrml <- paste(newfrml, "0", sep = "")
  }

  # Process refund-type terms
  where.refund <- c(where.af, where.lf, where.lf.vd, where.pr, where.fp, where.re)
  if (length(where.refund)) {
    fterms <- lapply(terms[where.refund], function(x) {
      eval(x, envir = evalenv, enclos = frmlenv)
    })
    newtrmstrings[where.refund] <- sapply(fterms, function(x) {
      safeDeparse(x$call)
    })
    lapply(fterms, function(x) {
      lapply(names(x$data), function(nm) {
        assign(x = nm, value = x$data[[nm]], envir = newfrmlenv)
        invisible(NULL)
      })
      # allow for unevaluated variable names in xt arguments:
      if ("xt" %in% names(x$call)) {
        xtvars <- all.vars(x$call$xt)
        if (length(xtvars)) {
          sapply(xtvars, function(xtvar) {
            xtvarval <- eval(as.name(xtvar), envir = evalenv, enclos = frmlenv)
            #assign into parent of newfrmlenv because these are not
            #necessarily covariates so list2df(newfrmlenv) below would fail
            assign(x = xtvar, value = xtvarval, envir = parent.env(newfrmlenv))
            invisible(NULL)
          })
        }
      }
      invisible(NULL)
    })
    fterms <- lapply(fterms, function(x) x[names(x) != "data"])
  }
  else fterms <- NULL

  # Process mgcv-type terms
  where.mgcv <- c(where.par, where.s, where.te, where.t2)
  if (length(where.mgcv)) {
    if ("data" %in% names(call))
      frmlenv <- list2env(eval(call$data), frmlenv)
    lapply(terms[where.mgcv], function(x) {
      nms <- if (!is.null(names(x))) {
        all.vars(x[names(x) == ""])
      }
      else all.vars(x)
      sapply(nms, function(nm) {
        stopifnot(length(get(nm, envir = frmlenv)) == nobs)
        assign(x = nm, value = get(nm, envir = frmlenv), envir = newfrmlenv)
        invisible(NULL)
      })
      invisible(NULL)
    })
  }

  # Finalize call to fitter
  newfrml <- formula(paste(newfrml, paste(newtrmstrings, collapse="+")))
  environment(newfrml) <- newfrmlenv
  pfrdata <- list2df(as.list(newfrmlenv))
  datameans <- sapply(as.list(newfrmlenv), function(x){
    if (is.numeric(x) | is.logical(x)) {
      mean(x)
    }else if (is.factor(x)){
      names(which.max(table(x)))
    }else
        NA
  }, simplify = FALSE)
  newcall <- expand.call(pfr, call)
  newcall$fitter  <- newcall$bs.int <- newcall$bs.yindex <- NULL
  newcall$formula <- newfrml
  newcall$data <- quote(pfrdata)
  newcall$method <- method
  newcall[[1]] <- fitter

  # Evaluate call
  res <- eval(newcall)

  # Post-process fit
  res.smooth <- if (as.character(fitter) %in% c("gamm4", "gamm")) {
    res$gam$smooth
  } else res$smooth
  names(res.smooth) <- sapply(res.smooth, function(x) x$label)

  if (as.character(fitter) %in% c("gamm4", "gamm")) {
    res$gam$smooth <- res.smooth
  } else {
    res$smooth <- res.smooth
  }

  termtype <- rep("par", length(terms))
  for (i in 1:length(specials))
    termtype[specials[[i]]-1] <- names(specials)[i]

  ret <- list(formula = formula,
              #termmap = trmmap, labelmap = labelmap,
              responsename = responsename, nobs = nobs,
              termnames = names(terms),
              termtype = termtype, datameans=datameans, ft = fterms)
  if (as.character(fitter) %in% c("gamm4", "gamm")) {
    res$gam$pfr <- ret
    class(res$gam) <- c("pfr", class(res$gam))
  }
  else {
    res$pfr <- ret
    class(res) <- c("pfr", class(res))
  }

  return(res)
}
