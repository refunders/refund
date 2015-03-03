#' Penalized Functional Regression
#' 
#' Implements various approaches to penalized scalar-on-function regression.
#' These techniques include Penalized Functional Regression (Goldsmith et al.,
#' 2011), Functional Principal Component Regression (Reiss and Ogden, 2007),
#' Functional Generalized Additive Models (McLean et al., 2013), and Partially
#' Empirical Eigenvectors for Regression (Randolph et al., 2012).
#' This function is a wrapper for mgcv's \code{\link{gam}}
#' and its siblings to fit models with a scalar (but necessarily continuous)
#' response.
#' 
#' @param formula a formula with special terms as for \code{\link{gam}},
#' with additional special terms \code{\link{af}}() and \code{\link{lf}}().
#' @param fitter the name of the function used to estimate the model. Defaults
#' to \code{\link{gam}} if the matrix of functional responses has less than 2e5
#' data points and to \code{\link{bam}} if not. "gamm" (see \code{\link{gamm}})
#' and "gamm4" (see \code{\link{gamm4}}) are valid options as well.
#' @param ... additional arguments that are valid for \code{\link{gam}} or
#' \code{\link{bam}}; for example, specify a \code{gamma} > 1 to increase amount
#' of smoothing when using GCV to choose smoothing parameters or
#' \code{method="REML"} to change to REML for estimation of smoothing parameters
#' (default is GCV).
#' @section Warning:
#' Binomial responses should be specified as a numeric vector rather than as a
#' matrix or a factor.
#' @return A fitted pfr-object, which is a \code{\link{gam}}-object with some
#' additional information in a pfr-element If fitter is "gamm" or "gamm4", only
#' the $gam part of the returned list is modified in this way.
#' @references McLean, M. W., Hooker, G., Staicu, A.-M., Scheipl, F., and
#' Ruppert, D. (2014). Functional generalized additive models. \emph{Journal of
#' Computational and Graphical Statistics}, \bold{23 (1)}, pp. 249-269. 
#' Available at \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3982924}.
#' @author Mathew W. McLean \email{mathew.w.mclean@@gmail.com} and Fabian Scheipl
#' @seealso \code{\link{af}}, \code{\link{lf}}
#' @importFrom mgcv gam gam.fit gamm4 bam s te t2
#' @importFrom gamm4 gamm4
#' @importFrom nlme lme4
#' @importFrom stats terms.formula
#' @export
#' 
pfr <- function(formula=NULL, fitter=NA, ...){
  
  if (class(formula) != "formula") {
    # Call pfr_old()
    fit <- if (is.null(formula) & is.na(fitter))
      pfr_old(...)
    else if (is.null(formula))
      pfr_old(subj=fitter, ...)
    else if (is.na(fitter))
      pfr_old(Y=formula, ...)
    else
      pfr_old(Y=formula, subj=fitter, ...)
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
  tf <- terms.formula(formula, specials = c("s", "te", "t2", "lf", "af",
                                            "lf.vd", "re"))
  trmstrings <- attr(tf, "term.labels")
  terms <- sapply(trmstrings, function(trm) as.call(parse(text = trm))[[1]], 
                  simplify = FALSE)
  frmlenv <- environment(formula)
  specials <- attr(tf, "specials")
  where.af <- specials$af - 1
  where.lf <- specials$lf - 1
  where.s  <- specials$s  - 1
  where.te <- specials$te - 1
  where.t2 <- specials$t2 - 1
  where.re <- specials$re - 1
  where.lf.vd <- specials$lf.vd - 1
  where.all <- c(where.af, where.lf, where.s, where.te, where.t2, where.re,
                 where.lf.vd)
  
  if (length(trmstrings)) {
    where.par <- which(!(1:length(trmstrings) %in% where.all))
  } else where.par <- numeric(0)
  
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
  
  assign(x = deparse(responsename),
         value = as.vector(t(eval(responsename, envir = evalenv,
                                  enclos = frmlenv))),
         envir = newfrmlenv)
  
  newtrmstrings <- attr(tf, "term.labels")
  if (!attr(tf, "intercept")) {
    newfrml <- paste(newfrml, "0", sep = "")
  }
  
  where.special <- c(where.af, where.lf, where.lf.vd, where.re)
  if (length(where.special)) {
    fterms <- lapply(terms[where.special], function(x) {
      eval(x, envir = evalenv, enclos = frmlenv)
    })
    newtrmstrings[where.special] <- sapply(fterms, function(x) {
      safeDeparse(x$call)
    })
    lapply(fterms, function(x) {
      lapply(names(x$data), function(nm) {
        assign(x = nm, value = x$data[[nm]], envir = newfrmlenv)
        invisible(NULL)
      })
      invisible(NULL)
    })
    fterms <- lapply(fterms, function(x) x[names(x) != "data"])
  }
  else fterms <- NULL
  
  where.notf <- c(where.par, where.s, where.te, where.t2)
  if (length(where.notf)) {
    if ("data" %in% names(call)) 
      frmlenv <- list2env(eval(call$data), frmlenv)
    lapply(terms[where.notf], function(x) {
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
  
  newfrml <- formula(paste(c(newfrml, newtrmstrings), collapse = "+"))
  environment(newfrml) <- newfrmlenv
  psfrdata <- list2df(as.list(newfrmlenv))
  datameans <- sapply(as.list(newfrmlenv),mean)
  newcall <- expand.call(pfr, call)
  newcall$fitter  <- newcall$bs.int <- newcall$bs.yindex <- NULL
  newcall$formula <- newfrml
  newcall$data <- quote(psfrdata)
  newcall[[1]] <- fitter
  
  # Evaluate call to fitter
  res <- eval(newcall)
  
  res.smooth <- if (as.character(fitter) %in% c("gamm4", "gamm")) {
    res$gam$smooth
  } else res$smooth
  
  trmmap <- newtrmstrings
  names(trmmap) <- names(terms)
  labelmap <- as.list(trmmap)
  names(trmmap) <- names(terms)
  
  lbls <- sapply(res.smooth, function(x) x$label)
  if (length(where.par)) {
    for (w in where.par) labelmap[[w]] <- {
      where <- sapply(res.smooth, function(x) x$by) == names(labelmap)[w]
      sapply(res.smooth[where], function(x) x$label)
    }
    labelmap[-c(where.par)] <- lbls[pmatch(sapply(labelmap[-c(where.par)], function(x) {
      tmp <- eval(parse(text = x))
      if (is.list(tmp)) {
        return(tmp$label)
      } else {
        return(x)
      }
    }), lbls)]
  } else {
    labelmap[1:length(labelmap)] <- lbls[pmatch(sapply(labelmap, function(x) {
      tmp <- eval(parse(text = x))
      if (is.list(tmp)) {
        return(tmp$label)
      } else {
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
  } else {
    res$smooth <- res.smooth
  }
  
  termtype <- rep("par", length(terms))
  for (i in 1:length(specials))
    termtype[specials[[i]]-1] <- names(specials)[i]
  
  ret <- list(formula = formula, termmap = trmmap, labelmap = labelmap, 
              responsename = responsename, nobs = nobs,
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
