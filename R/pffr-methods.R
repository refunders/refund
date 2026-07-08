# methods for pffr-objects
#
#
# Author: fabians
# 16.08.2011, 13:01:24
###############################################################################

#' Prediction for penalized function-on-function regression
#'
#'  Takes a fitted \code{pffr}-object produced by \code{\link{pffr}()} and produces
#'  predictions given a new set of values for the model covariates or the original
#'  values used for the model fit. Predictions can be accompanied by standard errors,
#'  based on the posterior distribution of the model coefficients. This is a wrapper
#'  function for \code{\link[mgcv]{predict.gam}()}.
#'
#'  Index variables (i.e., evaluation points) for the functional covariates are reused
#'  from the fitted model object and cannot be supplied with \code{newdata}.
#'  Prediction is always for the entire index range of the responses as defined
#'  in the original fit. If the original fit was performed on sparse or irregular,
#'  non-gridded response data supplied via \code{pffr}'s \code{ydata}-argument
#'  and no \code{newdata} was supplied, this function will
#'  simply return fitted values for the original evaluation points of the response (in list form).
#'  If the original fit was performed on sparse or irregular data and \code{newdata} \emph{was}
#'  supplied, the function will return predictions on the grid of evaluation points given in
#'  \code{object$pffr$yind}.
#'
#' @param object a fitted \code{pffr}-object
#' @param newdata  A named list (or a \code{data.frame}) containing the values of the
#' model covariates at which predictions are required.
#' If no \code{newdata} is provided then predictions corresponding to the original data
#' are returned. If \code{newdata} is provided then it must contain all the variables needed
#' for prediction, in the format supplied to \code{pffr}, i.e., functional predictors must be
#'  supplied as matrices with each row corresponding to one observed function.
#'  See Details for more on index variables and prediction for models fit on
#'  irregular or sparse data.
#' @param reformat logical, defaults to TRUE. Should predictions be returned in matrix form (default) or
#' in the long vector shape returned by \code{predict.gam()}?
#' @param type see \code{\link[mgcv]{predict.gam}()} for details.
#'  Note that \code{type == "lpmatrix"} will force \code{reformat} to FALSE.
#' @param se.fit see \code{\link[mgcv]{predict.gam}()}
#' @param ...  additional arguments passed on to \code{\link[mgcv]{predict.gam}()}
#' @seealso \code{\link[mgcv]{predict.gam}()}
#' @return If \code{type == "lpmatrix"}, the design matrix for the supplied covariate values in long format.
#'  If \code{se == TRUE}, a list with entries \code{fit} and \code{se.fit} containing fits and standard errors, respectively.
#'  If \code{type == "terms"} or \code{"iterms"} each of these lists is a list of matrices of the same dimension as the response for \code{newdata}
#'  containing the linear predictor and its se for each term.
#' @export
#' @method predict pffr
#' @author Fabian Scheipl
#' @importFrom mgcv predict.gam predict.bam
predict.pffr <- function(
  object,
  newdata,
  reformat = TRUE,
  type = "link",
  se.fit = FALSE,
  ...
) {
  #browser()

  call <- match.call()
  nyindex <- object$pffr$nyindex

  ## warn if any entries in ... are not arguments for predict.gam
  dots <- list(...)
  if (length(dots)) {
    validDots <- c(names(formals(predict.gam)), "cluster")
    # should be
    # unique(c(names(formals(predict.gam)),
    #          names(formals(predict.bam))))
    # but predict.bam is not exported.
    notUsed <- names(dots)[!(names(dots) %in% validDots)]
    if (length(notUsed))
      warning(
        "Arguments <",
        paste(notUsed, collapse = ", "),
        "> supplied but not used."
      )
  }

  if (!missing(newdata)) {
    nobs <- nrow(as.matrix(newdata[[1]]))

    # check if the supplied data already has the shape expected by predict.gam
    # and dispatch immediately if so (need this so summary works as expected!)
    if (
      !(all(names(newdata) %in% names(object$model))) |
        !(paste0(object$pffr$yind_name, ".vec") %in% names(newdata))
    ) {
      # check lengths
      stopifnot(
        length(unique(sapply(
          newdata,
          function(x) ifelse(is.matrix(x), nrow(x), length(x))
        ))) ==
          1
      )
      #        #FIXME: better leave this check to predict.gam....
      #        covnames <- mapply(gsub,
      #                pattern=c(".[st]mat$"),
      #                replacement="", x=unique(unlist(sapply(object$smooth, function(x) x$term))))
      #        covnames <- unique(covnames[covnames != paste(object$pffr$yind_name, ".vec", sep="")])
      #        stopifnot(all(covnames %in% names(newdata)))

      #get newdata into the shape expected by predict gam:
      gamdata <- list()
      #y-index
      gamdata[[paste(object$pffr$yind_name, ".vec", sep = "")]] <- rep(
        object$pffr$yind,
        times = nobs
      )

      # which covariates occur in which terms?
      varmap <- sapply(
        names(object$pffr$label_map),
        function(x) all.vars(formula(paste("~", x)))
      )

      # don't include response
      covnames <- unique(names(newdata)[
        names(newdata) != deparse(object$formula[[2]])
      ])
      for (cov in covnames) {
        #find the term(s) <cov> is associated with
        trms <- which(sapply(
          varmap,
          function(x) any(grep(paste("^", cov, "$", sep = ""), x))
        ))
        if (!is.null(dots$terms)) trms <- trms[names(trms) %in% dots$terms]
        if (length(trms) != 0) {
          for (trm in trms) {
            is.ff <- trm %in% object$pffr$where$ff
            is.sff <- trm %in% object$pffr$where$sff
            is.ffpc <- trm %in% object$pffr$where$ffpc
            is.pcre <- trm %in% object$pffr$where$pcre
            #if ff(X) or sff(X), generate (X.mat), X.tmat, X.smat, L.X ...
            if (is.ff) {
              ff <- object$pffr$ff[[grep(
                paste(cov, "[,\\)]", sep = ""),
                names(object$pffr$ff)
              )]]
              #... but don't generate new data unless <cov> is the functional covariate.
              if (
                grepl(paste(cov, "\\.[st]mat", sep = ""), deparse(ff$call$x))
              ) {
                # make L-matrix for new obs:
                L <- ff$L
                if (any(apply(L, 2, function(x) length(unique(x))) != 1)) {
                  stop(
                    "Error for ",
                    names(varmap)[trm],
                    "-- Prediction for ff-terms with varying rows in integration operator L not implememented yet."
                  )
                }

                predL <- matrix(
                  L[1, ],
                  byrow = TRUE,
                  nrow = nrow(newdata[[cov]]),
                  ncol = ncol(L)
                )

                # Create s and t matrices for new predictions
                smat <- matrix(
                  ff$xind,
                  byrow = TRUE,
                  ncol = length(ff$xind),
                  nrow = nobs * nyindex
                )
                tmat <- matrix(
                  rep(object$pffr$yind, times = nobs),
                  ncol = length(ff$xind),
                  nrow = nobs * nyindex
                )
                LX_stacked <- (predL * newdata[[cov]])[
                  rep(1:nobs, each = nyindex),
                ]

                if (!is.null(ff$limits)) {
                  # Apply limits: set weights to 0 outside integration region
                  use <- ff$limits(smat, tmat)
                  LX_stacked <- LX_stacked * use

                  # Find windows and reduce matrix size if possible
                  windows <- compute_integration_windows(use)
                  max_width <- max(windows[, 3])
                  if (max_width < ncol(smat)) {
                    eff_windows <- expand_windows_to_maxwidth(
                      windows,
                      ncol(smat)
                    )
                    smat <- shift_and_shorten_matrix(smat, eff_windows)
                    tmat <- shift_and_shorten_matrix(tmat, eff_windows)
                    LX_stacked <- shift_and_shorten_matrix(
                      LX_stacked,
                      eff_windows
                    )
                  }
                }

                gamdata[[paste(cov, ".smat", sep = "")]] <- smat
                gamdata[[paste(cov, ".tmat", sep = "")]] <- tmat
                gamdata[[paste("L.", cov, sep = "")]] <- LX_stacked
              }
            }
            if (is.sff) {
              sff <- object$pffr$ff[[grep(
                paste(cov, "[,\\)]", sep = ""),
                names(object$pffr$ff)
              )]]
              #... but don't generate new data unless <cov> is the functional covariate.
              if (
                grepl(paste(cov, "\\.[st]mat", sep = ""), deparse(sff$call$x))
              ) {
                # make L-matrix for new obs:
                L <- sff$L
                if (any(apply(L, 2, function(x) length(unique(x))) != 1)) {
                  stop(
                    "Error for ",
                    names(varmap)[trm],
                    "-- Prediction for sff-terms with varying rows in integration operator L not implememented yet."
                  )
                }
                predL <- matrix(
                  L[1, ],
                  byrow = TRUE,
                  nrow = nrow(newdata[[cov]]),
                  ncol = ncol(L)
                )

                gamdata[[paste(cov, ".mat", sep = "")]] <- newdata[[cov]][
                  rep(1:nobs, e = nyindex),
                ]
                gamdata[[paste(cov, ".smat", sep = "")]] <-
                  matrix(
                    sff$xind,
                    byrow = TRUE,
                    ncol = length(sff$xind),
                    nrow = nobs * nyindex
                  )
                gamdata[[paste(cov, ".tmat", sep = "")]] <-
                  matrix(
                    rep(object$pffr$yind, times = nobs),
                    ncol = length(sff$xind),
                    nrow = nobs * nyindex
                  )
                gamdata[[paste("L.", cov, sep = "")]] <- predL[
                  rep(1:nobs, e = nyindex),
                ]
              }
            }
            if (is.pcre) {
              pcre <- object$pffr$pcre[[grep(cov, names(object$pffr$pcre))]]
              gamdata[[paste(cov, ".vec", sep = "")]] <- rep(
                newdata[[cov]],
                each = nyindex
              )
              for (nm in colnames(pcre$efunctions)) {
                tmp <- approx(
                  x = pcre$yind,
                  y = pcre$efunctions[, nm],
                  xout = object$pffr$yind,
                  method = "linear"
                )$y
                gamdata[[nm]] <- tmp[rep(1:nyindex, times = nobs)]
              }
            }
            if (is.ffpc) {
              ffpc <- object$pffr$ffpc[[grep(
                paste(cov, "[,\\)]", sep = ""),
                names(object$pffr$ffpc)
              )]]
              # Xc' = Phi xi' + error --> get loadings for new data:
              Xct <- t(newdata[[cov]]) - as.vector(ffpc$meanX)
              xiMat <- t(qr.coef(qr(ffpc$PCMat), Xct))
              colnames(xiMat) <- paste(
                make.names(cov),
                ".PC",
                1:ncol(xiMat),
                sep = ""
              )
              xiMat <- xiMat[rep(1:nobs, each = nyindex), , drop = FALSE]
              for (nm in colnames(xiMat)) {
                gamdata[[nm]] <- xiMat[, nm]
              }
            }
            if (!(is.ff | is.sff | is.ffpc | is.pcre)) {
              gamdata[[cov]] <- if (!is.matrix(drop(newdata[[cov]]))) {
                #just repeat each entry nyindex-times to correspond to vec(<Response>)
                drop(newdata[[cov]])[rep(1:nobs, each = nyindex)]
              } else {
                # stack the matrix of the functional covariate (row-wise!)
                as.vector(t(newdata[[cov]]))
              }
            }
          }
        }
      }
      gamdata <- list2df(gamdata)
      call[["newdata"]] <- gamdata
    }
  } else {
    call$newdata <- eval(call$newdata)
    nobs <- object$pffr$nobs
  }
  isIrregular <- missing(newdata) & object$pffr$is_sparse

  # Honor the fit-time sandwich choice for standard errors: resolve the robust
  # covariance through the accessor and inject it into the local `object` so
  # predict.gam's se.fit uses it (the caller's object is not modified).
  # sandwich = "none" fits are left untouched so predict.gam sees the fit's
  # own model-based matrices exactly as mgcv produced them.
  fit_sandwich_type <- normalize_sandwich_type(
    object$pffr$sandwich_info$type %||% object$pffr$sandwich
  )
  if (
    isTRUE(se.fit) &&
      !identical(type, "lpmatrix") &&
      !identical(fit_sandwich_type, "none")
  ) {
    vsw <- pffr_vcov(object, sandwich = NULL)
    object$Vp <- vsw
    object$Vc <- vsw
  }

  #call predict.gam
  call[[1]] <- if (inherits(object, "bam")) {
    mgcv::predict.bam
  } else mgcv::predict.gam
  call$object <- as.name("object")
  ret <- eval(call)

  if (type == "lpmatrix" && reformat) {
    reformat <- FALSE
    warning("Setting reformat to FALSE for type=\"lpmatrix\".")
  }

  #reformat into matrices with same shape as <Response>

  if (reformat) {
    if (!isIrregular) {
      if (missing(newdata) && !is.null(object$pffr$missing_indices)) {
        #pad with NAs at the appropriate locations so that fits are nobs x nyindex:
        insertNA <- function(x) {
          if (length(x) != nobs * object$pffr$nyindex) {
            tmp <- rep(NA, nobs * object$pffr$nyindex)
            tmp[-object$pffr$missing_indices] <- x
            return(tmp)
          } else {
            return(x)
          }
        }
      } else insertNA <- function(x) return(x)

      if (se.fit) {
        if (type %in% c("terms", "iterms")) {
          ret <- lapply(
            ret,
            function(x)
              do.call(
                list,
                sapply(1:ncol(x), function(i) {
                  #browser()
                  d <- list(I(matrix(
                    insertNA(x[, i]),
                    nrow = nobs,
                    ncol = object$pffr$nyindex,
                    byrow = TRUE
                  )))
                  names(d) <- colnames(x)[i]
                  return(d)
                })
              )
          )
        } else {
          ret <- lapply(
            ret,
            function(x)
              matrix(
                insertNA(x),
                nrow = nobs,
                ncol = object$pffr$nyindex,
                byrow = TRUE
              )
          )
        }
      } else {
        if (type %in% c("terms", "iterms")) {
          ret <- do.call(
            list,
            sapply(1:ncol(ret), function(i) {
              #browser()
              d <- list(I(matrix(
                insertNA(ret[, i]),
                nrow = nobs,
                ncol = object$pffr$nyindex,
                byrow = TRUE
              )))
              names(d) <- colnames(ret)[i]
              return(d)
            })
          )
        } else
          ret <- matrix(
            insertNA(ret),
            nrow = nobs,
            ncol = object$pffr$nyindex,
            byrow = TRUE
          )
      }
    } else {
      evalpoints <- object$pffr$ydata[, c(".obs", ".index")]
      if (se.fit) {
        if (type %in% c("terms", "iterms")) {
          ret <- lapply(
            ret,
            function(x)
              do.call(
                list,
                sapply(1:ncol(x), function(i) {
                  #browser()
                  d <- list(cbind(evalpoints, .value = x[, i]))
                  names(d) <- colnames(x)[i]
                  return(d)
                })
              )
          )
        } else {
          ret <- lapply(ret, function(x) cbind(evalpoints, .value = x))
        }
      } else {
        if (type %in% c("terms", "iterms")) {
          ret <- do.call(
            list,
            sapply(1:ncol(ret), function(i) {
              #browser()
              d <- list(cbind(evalpoints, .value = ret[, i]))
              names(d) <- colnames(ret)[i]
              return(d)
            })
          )
        } else ret <- cbind(evalpoints, .value = ret)
      }
    }
  }
  return(ret)
}

#' Obtain model matrix for a pffr fit
#'
#' @param object a fitted \code{pffr}-object
#' @param ... other arguments, passed to \code{\link[mgcv]{predict.gam}}.
#'
#' @return A model matrix
#' @method model.matrix pffr
#' @author Fabian Scheipl
model.matrix.pffr <- function(object, ...) {
  if (!inherits(object, "pffr")) stop("`object' is not of class \"pffr\"")
  predict(object, type = "lpmatrix", reformat = FALSE, ...)
}

#' Obtain residuals and fitted values for a pffr models
#'
#' See \code{\link{predict.pffr}} for alternative options to extract estimated
#' values from a \code{pffr} object.
#' "Fitted values" here refers to the estimated additive predictor values,
#' these will not be on the scale of the response for models with link functions.
#'
#' For \code{family = "gaulss"} (Gaussian location-scale models), the fitted
#' values matrix has two columns: means and log-standard deviations. Use the
#' \code{which} argument in \code{fitted.pffr} to control which values are
#' returned.
#'
#' @param object a fitted \code{pffr}-object
#' @param reformat logical, defaults to TRUE. Should residuals/fitted values be returned in
#'   \code{n x yindex} matrix form (regular grid data) or, respectively, in the
#'   shape of the originally supplied \code{ydata} argument (sparse/irregular
#'   data), or, if \code{FALSE}, simply as a long vector as returned by
#'   \code{resid.gam()} or \code{fitted.gam()}?
#' @param which For \code{fitted.pffr} with \code{family = "gaulss"} only:
#'   which fitted values to return. One of \code{"mean"} (default, returns
#'   predicted means), \code{"scale"} (returns predicted log-standard
#'   deviations), or \code{"both"} (returns list with both components).
#' @param ... other arguments, passed to \code{\link[mgcv]{residuals.gam}}.
#'
#' @return A matrix or \code{ydata}-like \code{data.frame} or a vector of
#'   residuals / fitted values (see \code{reformat}-argument). For
#'   \code{fitted.pffr} with \code{family = "gaulss"} and \code{which = "both"},
#'   returns a list with \code{mean} and \code{scale} components.
#' @export
#' @importFrom mgcv residuals.gam
#' @method residuals pffr
#' @aliases fitted.pffr
#' @author Fabian Scheipl
residuals.pffr <- function(object, reformat = TRUE, ...) {
  if (!inherits(object, "pffr")) stop("`object' is not of class \"pffr\"")
  ret <- mgcv::residuals.gam(object, ...)
  if (reformat) {
    if (!object$pffr$is_sparse) {
      if (!(length(ret) == object$pffr$nobs * object$pffr$nyindex)) {
        tmp <- rep(NA, object$pffr$nobs * object$pffr$nyindex)
        tmp[-object$pffr$missing_indices] <- ret
        ret <- tmp
      }
      ret <- matrix(
        ret,
        nrow = object$pffr$nobs,
        ncol = object$pffr$nyindex,
        byrow = TRUE
      )
    } else {
      tmp <- object$pffr$ydata
      tmp[, ".value"] <- ret
      ret <- tmp
    }
  }
  return(ret)
}

#' @method fitted pffr
#' @export
#' @rdname residuals.pffr
fitted.pffr <- function(
  object,
  reformat = TRUE,
  which = c("mean", "scale", "both"),
  ...
) {
  if (!inherits(object, "pffr")) {
    stop("`object' is not of class \"pffr\"")
  }
  which <- match.arg(which)

  ret <- object$fitted.values
  is_gaulss <- object$family$family == "gaulss"

  # Helper to reformat a single vector of fitted values
  reformat_fitted <- function(vals) {
    if (!object$pffr$is_sparse) {
      if (!(length(vals) == object$pffr$nobs * object$pffr$nyindex)) {
        tmp <- rep(NA, object$pffr$nobs * object$pffr$nyindex)
        tmp[-object$pffr$missing_indices] <- vals
        vals <- tmp
      }
      matrix(
        vals,
        nrow = object$pffr$nobs,
        ncol = object$pffr$nyindex,
        byrow = TRUE
      )
    } else {
      tmp <- object$pffr$ydata
      tmp[, ".value"] <- vals
      tmp
    }
  }

  if (is_gaulss && is.matrix(ret) && ncol(ret) >= 2) {
    # gaulss: column 1 = mean, column 2 = log(sd)
    mean_vals <- ret[, 1]
    scale_vals <- ret[, 2]

    if (reformat) {
      mean_mat <- reformat_fitted(mean_vals)
      scale_mat <- reformat_fitted(scale_vals)

      ret <- switch(
        which,
        mean = mean_mat,
        scale = scale_mat,
        both = list(mean = mean_mat, scale = scale_mat)
      )
    } else {
      ret <- switch(
        which,
        mean = mean_vals,
        scale = scale_vals,
        both = list(mean = mean_vals, scale = scale_vals)
      )
    }
  } else {
    # Non-gaulss or single-column case
    if (which != "mean" && !is_gaulss) {
      warning("'which' argument is ignored for non-gaulss families")
    }
    if (reformat) {
      ret <- reformat_fitted(ret)
    }
  }

  ret
}

#' Plot a pffr fit
#'
#' Plot a fitted pffr-object. Simply dispatches to \code{\link[mgcv]{plot.gam}}.
#'
#' @param x a fitted \code{pffr}-object
#' @param ... arguments handed over to \code{\link[mgcv]{plot.gam}}
#'
#' @return This function only generates plots (and invisibly returns
#'   \code{\link[mgcv]{plot.gam}}'s plot data).
#' @method plot pffr
#' @export
#' @importFrom mgcv plot.gam
#' @author Fabian Scheipl
plot.pffr <- function(x, ...) {
  call <- match.call()
  call[[1]] <- mgcv::plot.gam
  # Honor the fit-time sandwich choice: inject the resolved robust covariance
  # into the local `x` so plot.gam's standard-error bands use it. sandwich =
  # "none" fits keep their model-based matrices untouched.
  fit_sandwich_type <- normalize_sandwich_type(
    x$pffr$sandwich_info$type %||% x$pffr$sandwich
  )
  if (!identical(fit_sandwich_type, "none")) {
    vsw <- pffr_vcov(x, sandwich = NULL)
    x$Vp <- vsw
    x$Vc <- vsw
  }
  #drop "pffr" class and replace <object> with changed value s.t. method dispatch works without glitches
  class(x) <- class(x)[-1]
  # point the call at the modified local object (as in summary.pffr/
  # predict.pffr); otherwise eval() would resolve the caller's original
  # object and the class/covariance changes above would be ignored.
  call$x <- as.name("x")
  invisible(eval(call))
}


# -----------------------------------------------------------------------------
# Helper functions for coef.pffr (extracted for clarity and testability)
# -----------------------------------------------------------------------------

#' Safely compute range for coef.pffr
#'
#' Returns NA range for factors, otherwise numeric range.
#'
#' @param x A vector (possibly factor).
#' @returns Numeric vector of length 2 with range or c(NA, NA) for factors.
#' @keywords internal
coef_safe_range <- function(x) {
  if (is.factor(x)) return(c(NA, NA))
  range(x, na.rm = TRUE)
}

#' Generate evaluation grid for smooth term
#'
#' Creates a data frame grid over the range of the covariates for coefficient
#' evaluation.
#'
#' @param trm A smooth term object from object$smooth.
#' @param model_data The model data frame (object$model).
#' @param pffr_info List with pffr metadata: yind_name, pcre_terms.
#' @param grid_sizes Named list with n1, n2, n3 grid sizes.
#' @param is_pcre Logical, is this a pcre term?
#' @returns A data frame suitable for PredictMat, with xm/ym/zm attributes.
#' @keywords internal
#' @importFrom mgcv get.var
coef_make_data_grid <- function(
  trm,
  model_data,
  pffr_info,
  grid_sizes,
  is_pcre
) {
  x <- get.var(trm$term[1], model_data)

  # 1-dimensional smooth

  if (trm$dim == 1) {
    xg <- if (is.factor(x)) unique(x) else
      seq(min(x), max(x), length = grid_sizes$n1)
    d <- data.frame(xg)
    colnames(d) <- trm$term
    attr(d, "xm") <- xg
    return(finalize_grid_by_var(d, trm))
  }

  # PCRE term (special case)
  if (is_pcre) {
    ng <- grid_sizes$n2
    xg <- if (is.factor(x)) unique(x) else seq(min(x), max(x), length = ng)

    which_pcre <- which(
      sapply(pffr_info$pcre_terms, `[[`, "idname") == trm$term[1]
    )
    pcreterm <- pffr_info$pcre_terms[[which_pcre]]
    yg <- seq(min(pcreterm$yind), max(pcreterm$yind), l = ng)

    # Interpolate eigenfunctions to grid values
    efcts_grid <- sapply(colnames(pcreterm$efunctions), function(nm) {
      approx(
        x = pcreterm$yind,
        y = pcreterm$efunctions[, nm],
        xout = yg,
        method = "linear"
      )$y
    })
    efcts_grid <- data.frame(efcts_grid[rep(1:ng, each = length(xg)), ])
    colnames(efcts_grid) <- colnames(pcreterm$efunctions)

    d <- cbind(expand.grid(xg, yg), efcts_grid)
    colnames(d)[1:2] <- c(trm$term[1], paste0(pffr_info$yind_name, ".vec"))
    attr(d, "xm") <- xg
    attr(d, "ym") <- yg
    return(finalize_grid_by_var(d, trm))
  }

  # Multi-dimensional smooth (dim > 1)
  ng <- if (trm$dim == 2) grid_sizes$n2 else grid_sizes$n3

  xg <- if (is.factor(x)) unique(x) else seq(min(x), max(x), length = ng)
  y <- get.var(trm$term[2], model_data)
  yg <- if (is.factor(y)) unique(y) else seq(min(y), max(y), length = ng)

  if (length(trm$term) == 2) {
    d <- expand.grid(xg, yg)
    attr(d, "xm") <- xg
    attr(d, "ym") <- yg
  } else {
    z <- get.var(trm$term[3], model_data)
    zg <- if (is.factor(z)) unique(z) else seq(min(z), max(z), length = ng)
    d <- expand.grid(xg, yg, zg)
    attr(d, "xm") <- xg
    attr(d, "ym") <- yg
    attr(d, "zm") <- zg
  }
  colnames(d) <- trm$term
  finalize_grid_by_var(d, trm)
}

#' Add by-variable column to grid if needed
#'
#' @param d Data frame grid.
#' @param trm Smooth term object.
#' @returns Modified data frame with by column set to 1 if applicable.
#' @keywords internal
finalize_grid_by_var <- function(d, trm) {
  if (trm$by != "NA") {
    d$by <- 1
    colnames(d) <- c(head(colnames(d), -1), trm$by)
  }
  d
}

#' Resolve optional fixed evaluation grid for one smooth term
#'
#' @param eval_grid Optional list of per-term evaluation grids.
#' @param smooth_names Character vector of smooth names from object$smooth.
#' @param i Smooth term index.
#' @returns Data frame grid or NULL.
#' @keywords internal
resolve_eval_grid_for_term <- function(eval_grid, smooth_names, i) {
  if (is.null(eval_grid)) return(NULL)
  if (!is.list(eval_grid)) {
    stop("'eval_grid' must be a list when supplied.")
  }

  d <- NULL
  if (!is.null(names(eval_grid)) && smooth_names[i] %in% names(eval_grid)) {
    d <- eval_grid[[smooth_names[i]]]
  }
  if (is.null(d) && length(eval_grid) >= i) {
    d <- eval_grid[[i]]
  }
  d
}

#' Ensure evaluation grid has axis attributes used by coef extraction
#'
#' @param d Evaluation grid data frame.
#' @param trm Smooth term object.
#' @param is_pcre Logical, is this a pcre term?
#' @param pffr_info List with yind_name.
#' @returns Data frame with xm/ym/zm attributes set.
#' @keywords internal
ensure_grid_axis_attributes <- function(d, trm, is_pcre, pffr_info) {
  if (!is.data.frame(d)) {
    stop("Each 'eval_grid' entry must be a data.frame.")
  }
  if (is.null(attr(d, "xm")) && trm$term[1] %in% names(d)) {
    attr(d, "xm") <- unique(d[[trm$term[1]]])
  }

  if (trm$dim >= 2 || is_pcre) {
    y_term <- if (is_pcre) paste0(pffr_info$yind_name, ".vec") else trm$term[2]
    if (is.null(attr(d, "ym")) && y_term %in% names(d)) {
      attr(d, "ym") <- unique(d[[y_term]])
    }
  }

  if (trm$dim >= 3 && length(trm$term) >= 3 && trm$term[3] %in% names(d)) {
    if (is.null(attr(d, "zm"))) {
      attr(d, "zm") <- unique(d[[trm$term[3]]])
    }
  }

  d
}

#' Compute predictions for coefficient extraction
#'
#' Evaluates smooth term on grid and computes coefficients and standard errors.
#'
#' @param trm Smooth term object.
#' @param data_grid Data frame from coef_make_data_grid.
#' @param object_info List with: coefficients, cmX, Vp.
#' @param pffr_info List with: yind_name.
#' @param covmat Covariance matrix for SE computation.
#' @param se Logical, compute standard errors?
#' @param seWithMean Logical, include mean uncertainty?
#' @param is_pcre Logical, is this a pcre term?
#' @param ci One of "none", "pointwise", "simultaneous".
#' @param level Confidence level for intervals.
#' @param coef_draws Simulated coefficient perturbations (for simultaneous CIs).
#' @param t_scale Optional multiplier-t scaling, one positive value per
#'   simulation draw.
#' @returns List with x, y, z coordinates, value, se, coef data frame, dim.
#' @keywords internal
#' @importFrom mgcv PredictMat
coef_get_predictions <- function(
  trm,
  data_grid,
  object_info,
  pffr_info,
  covmat,
  se,
  seWithMean,
  is_pcre,
  ci = "none",
  level = 0.95,
  coef_draws = NULL,
  t_scale = NULL,
  crit_mode = "z",
  crit_df_const = NA_real_,
  df_ctx = NULL
) {
  X <- PredictMat(trm, data_grid)

  # For pcre terms, temporarily adjust term for axis setup
  if (is_pcre) {
    trm$dim <- 2
    trm$term[2] <- paste0(pffr_info$yind_name, ".vec")
  }

  # Build result structure based on dimensionality
  P <- build_coef_axes(trm, data_grid)

  # Compute predicted values
  trmind <- trm$first.para:trm$last.para
  P$value <- X %*% object_info$coefficients[trmind]
  P$coef <- cbind(data_grid, value = P$value)

  # Compute standard errors and intervals if requested
  if (se) {
    linear_map <- build_coef_linear_map(
      X = X,
      trmind = trmind,
      trm = trm,
      object_info = object_info,
      seWithMean = seWithMean
    )
    P$se <- compute_coef_se(linear_map = linear_map, covmat = covmat)
    P$coef <- cbind(P$coef, se = P$se)

    if (ci == "simultaneous") {
      crit <- compute_ci_critical(
        ci = ci,
        level = level,
        se_vec = P$se,
        linear_map = linear_map,
        coef_draws = coef_draws,
        t_scale = t_scale
      )
      P$crit <- crit
      ci_half <- crit * P$se
      P$coef <- cbind(
        P$coef,
        lower = P$value - ci_half,
        upper = P$value + ci_half
      )
    } else if (ci == "pointwise") {
      # Per-point pointwise reference (S3): crit may be a vector when
      # crit_mode = "satterthwaite" (a per-point Bell-McCaffrey df).
      pw <- compute_pointwise_ci(
        crit_mode = crit_mode,
        level = level,
        linear_map = linear_map,
        df_ctx = df_ctx,
        crit_df_const = crit_df_const
      )
      P$crit <- pw$crit
      ci_half <- pw$crit * P$se
      P$coef <- cbind(
        P$coef,
        lower = P$value - ci_half,
        upper = P$value + ci_half,
        df = pw$df
      )
    }
  }

  P$dim <- trm$dim
  P
}

#' Build coordinate axes for coefficient results
#'
#' @param trm Smooth term with dim attribute.
#' @param data_grid Grid data frame with xm/ym/zm attributes.
#' @returns List with x, y, z, xlab, ylab, zlab, xlim, ylim, zlim as appropriate.
#' @keywords internal
build_coef_axes <- function(trm, data_grid) {
  if (trm$dim == 1) {
    return(list(
      x = attr(data_grid, "xm"),
      xlab = trm$term,
      xlim = coef_safe_range(attr(data_grid, "xm"))
    ))
  }

  if (trm$dim == 2) {
    return(list(
      x = attr(data_grid, "xm"),
      y = attr(data_grid, "ym"),
      xlab = trm$term[1],
      ylab = trm$term[2],
      xlim = coef_safe_range(attr(data_grid, "xm")),
      ylim = coef_safe_range(attr(data_grid, "ym"))
    ))
  }

  if (trm$dim == 3) {
    return(list(
      x = attr(data_grid, "xm"),
      y = attr(data_grid, "ym"),
      z = attr(data_grid, "zm"),
      xlab = trm$term[1],
      ylab = trm$term[2],
      zlab = trm$term[3],
      xlim = coef_safe_range(attr(data_grid, "xm")),
      ylim = coef_safe_range(attr(data_grid, "ym")),
      zlim = coef_safe_range(attr(data_grid, "zm"))
    ))
  }

  NULL
}

#' Build linear map used for smooth-term uncertainty calculations
#'
#' @param X Prediction matrix for the smooth term.
#' @param trmind Index vector for term parameters.
#' @param trm Smooth term object (for nCons and meanL1).
#' @param object_info List with cmX, Vp.
#' @param seWithMean Logical, use seWithMean approach?
#' @returns List with linear map matrix and indexing metadata.
#' @keywords internal
build_coef_linear_map <- function(X, trmind, trm, object_info, seWithMean) {
  if (seWithMean && attr(trm, "nCons") > 0) {
    message("using seWithMean for ", trm$label, ".")
    X1 <- matrix(object_info$cmX, nrow(X), ncol(object_info$Vp), byrow = TRUE)
    meanL1 <- trm$meanL1
    if (!is.null(meanL1)) X1 <- X1 / meanL1
    X1[, trmind] <- X
    return(list(X = X1, use_full = TRUE, trmind = trmind))
  }
  list(X = X, use_full = FALSE, trmind = trmind)
}

#' Compute standard errors for coefficient extraction
#'
#' @param linear_map List returned by build_coef_linear_map().
#' @param covmat Covariance matrix.
#' @returns Numeric vector of standard errors.
#' @keywords internal
compute_coef_se <- function(linear_map, covmat) {
  if (linear_map$use_full) {
    return(sqrt(rowSums((linear_map$X %*% covmat) * linear_map$X)))
  }
  trmind <- linear_map$trmind
  sqrt(rowSums((linear_map$X %*% covmat[trmind, trmind]) * linear_map$X))
}

#' Draw coefficient perturbations for simultaneous intervals
#'
#' @param covmat Covariance matrix.
#' @param n_sim Number of simulation draws.
#' @param sim_seed Optional integer seed.
#' @param df Optional degrees of freedom for multiplier-t scaling.
#' @returns Matrix with one simulated perturbation vector per column, or a list
#'   with elements `draws` and `t_scale` when `df` is supplied.
#' @keywords internal
draw_coef_perturbations <- function(covmat, n_sim, sim_seed = NULL, df = NULL) {
  if (!is.null(sim_seed)) {
    has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (has_seed)
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(
      {
        if (has_seed) {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        } else if (
          exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        ) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      },
      add = TRUE
    )
    set.seed(sim_seed)
  }

  cov_sym <- 0.5 * (covmat + t(covmat))
  eig <- eigen(cov_sym, symmetric = TRUE)
  eig$values <- pmax(eig$values, 0)
  root_cov <- eig$vectors %*% diag(sqrt(eig$values), nrow = length(eig$values))

  draws <- root_cov %*%
    matrix(
      stats::rnorm(ncol(covmat) * n_sim),
      nrow = ncol(covmat),
      ncol = n_sim
    )
  if (is.null(df)) return(draws)
  if (!is.numeric(df) || length(df) != 1 || !is.finite(df) || df <= 0) {
    stop("`df` must be a single positive finite number.", call. = FALSE)
  }
  list(
    draws = draws,
    t_scale = sqrt(df / stats::rchisq(n_sim, df = df))
  )
}

#' Compute critical value for pointwise/simultaneous intervals
#'
#' @param ci One of "none", "pointwise", "simultaneous".
#' @param level Confidence level.
#' @param se_vec Standard errors for one term.
#' @param linear_map List returned by build_coef_linear_map().
#' @param coef_draws Simulated coefficient perturbations for simultaneous CIs.
#' @param t_scale Optional multiplier-t scaling, one positive value per
#'   simulation draw.
#' @returns Scalar critical value.
#' @keywords internal
compute_ci_critical <- function(
  ci,
  level,
  se_vec,
  linear_map,
  coef_draws = NULL,
  t_scale = NULL
) {
  if (ci == "none") return(NA_real_)
  if (ci == "pointwise") return(stats::qnorm((1 + level) / 2))

  eps <- sqrt(.Machine$double.eps)
  valid <- is.finite(se_vec) & (se_vec > eps)
  if (!any(valid)) return(0)

  if (is.null(coef_draws)) {
    stop("coef_draws must be supplied for simultaneous intervals.")
  }
  if (!is.null(t_scale)) {
    if (
      !is.numeric(t_scale) ||
        length(t_scale) != ncol(coef_draws) ||
        any(!is.finite(t_scale)) ||
        any(t_scale <= 0)
    ) {
      stop(
        "`t_scale` must be a positive finite vector with one value per simulation draw.",
        call. = FALSE
      )
    }
  }

  term_draws <- if (linear_map$use_full) {
    linear_map$X %*% coef_draws
  } else {
    linear_map$X %*% coef_draws[linear_map$trmind, , drop = FALSE]
  }

  max_stat <- apply(
    abs(term_draws[valid, , drop = FALSE] / se_vec[valid]),
    2,
    max
  )
  if (!is.null(t_scale)) max_stat <- max_stat * t_scale
  as.numeric(stats::quantile(
    max_stat,
    probs = level,
    names = FALSE,
    type = 8,
    na.rm = TRUE
  ))
}

#' Expand a term's linear map to full-coefficient-space contrasts
#'
#' Returns the `n_points x p` matrix whose rows are the contrast vectors `a`
#' (one per evaluation point) in the full model coefficient space, as consumed
#' by the Satterthwaite df kernel. For the `seWithMean` path the linear map is
#' already full-width; otherwise the term columns are scattered into a zero
#' matrix at the term's coefficient indices.
#'
#' @param linear_map List returned by [build_coef_linear_map()].
#' @param p Number of model coefficients (full covariance dimension).
#' @returns An `n_points x p` contrast matrix.
#' @keywords internal
pointwise_full_contrasts <- function(linear_map, p) {
  if (isTRUE(linear_map$use_full)) {
    return(linear_map$X)
  }
  Xp <- matrix(0, nrow = nrow(linear_map$X), ncol = p)
  Xp[, linear_map$trmind] <- linear_map$X
  Xp
}

#' Per-point pointwise critical value and reference df
#'
#' Computes the pointwise interval half-width multiplier for one term under the
#' chosen reference (S3): `"z"` (Gaussian, reported df `Inf`), `"tG1"`
#' (\eqn{t_{G-1}}, constant df), or `"satterthwaite"` (per-point Bell-McCaffrey
#' df from `df_ctx`, with a Gaussian fallback at any zero-variance / undefined
#' point).
#'
#' @param crit_mode One of `"z"`, `"tG1"`, `"satterthwaite"`.
#' @param level Confidence level.
#' @param linear_map List returned by [build_coef_linear_map()] (only needed for
#'   `"satterthwaite"`).
#' @param df_ctx A [pffr_df_context()] result (only needed for
#'   `"satterthwaite"`).
#' @param crit_df_const Constant df for `"tG1"` (\eqn{G-1}).
#' @returns A list with `crit` (scalar or per-point vector) and `df` (per-point
#'   vector; `Inf` for `"z"`).
#' @keywords internal
compute_pointwise_ci <- function(
  crit_mode,
  level,
  linear_map,
  df_ctx = NULL,
  crit_df_const = NA_real_
) {
  prob <- (1 + level) / 2
  n <- nrow(linear_map$X)
  if (crit_mode == "z" || is.null(df_ctx) && crit_mode == "satterthwaite") {
    return(list(crit = stats::qnorm(prob), df = rep(Inf, n)))
  }
  if (crit_mode == "tG1") {
    return(list(
      crit = stats::qt(prob, crit_df_const),
      df = rep(crit_df_const, n)
    ))
  }
  # satterthwaite: per-point df from the full-space contrasts and df context.
  Xp <- pointwise_full_contrasts(linear_map, p = ncol(df_ctx$Vp))
  df <- pffr_df_from_context(df_ctx, Xp)
  crit <- ifelse(
    is.finite(df),
    stats::qt(prob, pmax(df, 1)),
    stats::qnorm(prob)
  )
  list(crit = crit, df = df)
}


#' Get estimated coefficients from a pffr fit
#'
#' Returns estimated coefficient functions/surfaces \eqn{\beta(t), \beta(s,t)}
#' and estimated smooth effects \eqn{f(z), f(x,z)} or \eqn{f(x, z, t)} and their point-wise estimated standard errors.
#' Not implemented for smooths in more than 3 dimensions.
#'
#' The \code{seWithMean}-option corresponds to the \code{"iterms"}-option in \code{\link[mgcv]{predict.gam}}.
#' The \code{sandwich}-option computes robust standard errors. With
#' \code{sandwich="cluster"}, a cluster-robust sandwich (clustering by curve)
#' is used, which handles both heteroskedasticity and within-curve correlation.
#' With \code{sandwich="cl2"}, a leverage-adjusted cluster-robust sandwich
#' (Bell-McCaffrey style CL2) is used.
#' With \code{sandwich="hc"}, mgcv's observation-level HC sandwich is used.
#' If the model was fitted with a matching sandwich option in
#' \code{\link{pffr}}, the pre-computed covariance matrices are used directly.
#'
#'
#' @param object a fitted \code{pffr}-object
#' @param raw logical, defaults to FALSE. If TRUE, the function simply returns \code{object$coefficients}
#' @param se logical, defaults to TRUE. Return estimated standard error of the estimates?
#' @param freq logical, defaults to FALSE. If FALSE, use Bayesian posterior covariance for
#'   variability estimates: \code{object$Vc} if available (includes correction for smoothing
#'   parameter uncertainty), otherwise \code{object$Vp}. If TRUE, use frequentist
#'   covariance \code{object$Ve}. See \code{\link[mgcv]{gamObject}}.
#' @param sandwich Type of sandwich-corrected covariance for standard errors.
#'   \code{"cluster"} (default): cluster-robust sandwich (clustering by
#'   curve).
#'   \code{"cl2"}: leverage-adjusted cluster-robust sandwich (clustering by
#'   curve).
#'   \code{"hc"}: observation-level HC sandwich via \code{\link[mgcv]{vcov.gam}}.
#'   \code{"none"}: use model's default covariance.
#'   If the model was fitted with a matching sandwich type, the pre-computed
#'   covariance matrices are used directly.
#' @param cluster optional grouping for the cluster-robust sandwich
#'   (\code{sandwich = "cluster"} or \code{"cl2"}): a vector with one entry per
#'   curve (functional observation) mapping each curve to its independent unit.
#'   Defaults to \code{NULL}, i.e. each curve is its own cluster. Supply this for
#'   nested / repeated-measures designs where several curves share a higher-level
#'   unit (e.g. a subject id with multiple visits), so the sandwich clusters at
#'   the correct level. Only supported for densely-observed responses. When
#'   supplied, the pre-computed-covariance shortcut is bypassed.
#' @param dof_correction Optional CR1 small-sample dof correction for
#'   \code{sandwich = "cluster"}: \code{"none"} or \code{"edf"} (see
#'   \code{\link{pffr}}). Defaults to \code{NULL}, i.e. inherit whatever the
#'   model was fitted with, so the stored covariance is reused. Supplying a value
#'   that differs from the fit forces recomputation with the requested option.
#'   Ignored (with a warning) for \code{sandwich} other than \code{"cluster"}.
#' @param edf_type Which EDF the \code{"edf"} correction uses
#'   (\code{"trace"}/\code{"edf2"}/\code{"basis"}; see \code{\link{pffr}}).
#'   Defaults to \code{NULL} (inherit from the fit).
#' @param seWithMean logical, defaults to TRUE. Include uncertainty about the intercept/overall mean in  standard errors returned for smooth components?
#' @param n1 see below
#' @param n2 see below
#' @param n3 \code{n1, n2, n3} give the number of gridpoints for 1-/2-/3-dimensional smooth terms
#' used in the marginal equidistant grids over the range of the covariates at which the estimated effects are evaluated.
#' @param ci Type of confidence intervals to return in addition to standard
#'   errors. One of \code{"none"} (default), \code{"pointwise"}, or
#'   \code{"simultaneous"}.
#' @param ci_ref Reference distribution for \code{ci = "simultaneous"}.
#'   \code{"t"} (default) uses a finite-sample
#'   \eqn{t_{G-1}}{t_(G-1)} multiplier reference, where \eqn{G} is the number of
#'   independent curves or user-supplied clusters. This widens simultaneous
#'   bands at small \eqn{G} and converges to the Gaussian multiplier reference
#'   as \eqn{G} grows. \code{"normal"} restores the previous Gaussian
#'   multiplier reference. Pointwise intervals are governed instead by
#'   \code{crit} (below).
#' @param crit Reference distribution for the \emph{pointwise} critical value
#'   (\code{ci = "pointwise"}); the pointwise counterpart of \code{ci_ref}.
#'   \code{"auto"} (default) uses the per-point Satterthwaite reference when the
#'   standard errors come from a cluster-robust sandwich
#'   (\code{sandwich = "cluster"} or \code{"cl2"}) and the number of independent
#'   curves/clusters is moderate (\eqn{G < 150}), and the Gaussian reference
#'   otherwise. \code{"z"} always uses the Gaussian quantile (the historical
#'   behaviour). \code{"tG1"} uses a \eqn{t_{G-1}}{t_(G-1)} reference (constant
#'   df, the pointwise analogue of \code{ci_ref = "t"}). \code{"satterthwaite"}
#'   uses the per-point Bell-McCaffrey (Satterthwaite) df, \eqn{\nu(a) =
#'   (\sum_g \lVert q_g\rVert^2)^2 / \sum_g \lVert q_g\rVert^4} with
#'   \eqn{q_g = A_g \tilde X_g V_p a}; requested on a non-cluster covariance it
#'   degrades to \code{"z"} with a warning. This is the missing (df) half of the
#'   CL2 leverage adjustment. \strong{Honesty note:} the df uses a working-iid
#'   Satterthwaite shortcut that drops the same cross-cluster residual terms as
#'   the shipped \eqn{(I-H_{gg})^{-1/2}} CL2 covariance (paper Appendix C); it
#'   therefore returns \eqn{\approx G} for a perfectly balanced design where the
#'   exact Bell-McCaffrey df is \eqn{G-1} (the exact-BM df is future work).
#'   Simultaneous bands are unaffected.
#' @param level Confidence level for confidence intervals, defaults to
#'   \code{0.95}.
#' @param n_sim Number of simulations for simultaneous intervals, defaults to
#'   \code{2000}. Ignored unless \code{ci = "simultaneous"}.
#' @param sim_seed Optional integer seed for simultaneous interval simulation.
#' @param ... other arguments, not used.
#'
#' @return If \code{raw==FALSE}, a list containing \itemize{
#'  \item \code{pterms} a matrix containing the parametric / non-functional coefficients (and, optionally, their se's)
#'  \item \code{smterms} a named list with one entry for each smooth term in the model. Each entry contains
#'     \itemize{
#'          \item \code{coef} a matrix giving the grid values over the covariates, the estimated effect (and, optionally, the se's).
#'                          The first covariate varies the fastest.
#'          \item \code{x, y, z} the unique gridpoints used to evaluate the smooth/coefficient function/coefficient surface
#'          \item \code{xlim, ylim, zlim} the extent of the x/y/z-axes
#'          \item \code{xlab, ylab, zlab} the names of the covariates for the x/y/z-axes
#'          \item \code{dim} the dimensionality of the effect
#'          \item \code{main} the label of the smooth term (a short label, same as the one used in \code{summary.pffr})
#' }}
#' If \code{ci != "none"}, the returned matrices include columns \code{lower}
#' and \code{upper}. For \code{ci = "pointwise"} they also include a \code{df}
#' column giving the per-point reference degrees of freedom used for the
#' critical value (\code{Inf} for \code{crit = "z"}, \eqn{G-1} for
#' \code{crit = "tG1"}, and the per-point Satterthwaite df for
#' \code{crit = "satterthwaite"}/\code{"auto"}). The returned list also includes
#' \code{ci_meta} with CI settings (including \code{crit} and the resolved
#' \code{crit_used}).
#' @method coef pffr
#' @export
#' @importFrom mgcv PredictMat get.var
#' @seealso \code{\link[mgcv]{plot.gam}}, \code{\link[mgcv]{predict.gam}} which this routine is
#'   based on.
#' @author Fabian Scheipl
coef.pffr <- function(
  object,
  raw = FALSE,
  se = TRUE,
  freq = FALSE,
  sandwich = c("cluster", "cl2", "hc", "none"),
  cluster = NULL,
  dof_correction = NULL,
  edf_type = NULL,
  seWithMean = TRUE,
  n1 = 100,
  n2 = 40,
  n3 = 20,
  ci = c("none", "pointwise", "simultaneous"),
  ci_ref = c("t", "normal"),
  crit = c("auto", "z", "tG1", "satterthwaite"),
  level = 0.95,
  n_sim = 2000,
  sim_seed = NULL,
  ...
) {
  sandwich_missing <- missing(sandwich)
  # Backward compat: TRUE -> "cluster", FALSE -> "none"
  if (is.logical(sandwich)) sandwich <- if (sandwich) "cluster" else "none"
  sandwich <- match.arg(sandwich)
  ci <- match.arg(ci)
  ci_ref <- match.arg(ci_ref)
  crit <- match.arg(crit)

  # dof_correction / edf_type default to inheriting whatever the model was
  # fitted with (so coef() with no override returns the stored covariance);
  # an explicit value triggers recomputation (see cache logic below).
  dof_explicitly_set <- !is.null(dof_correction)
  model_dof_correction <- object$pffr$dof_correction %||% "none"
  model_edf_type <- object$pffr$edf_type %||% "trace"
  if (is.null(dof_correction)) dof_correction <- model_dof_correction
  if (is.null(edf_type)) edf_type <- model_edf_type
  dof_correction <- match.arg(dof_correction, c("none", "edf"))
  edf_type <- match.arg(edf_type, c("trace", "edf2", "basis"))
  # Only warn when the user *explicitly* asked for a dof correction on a
  # non-cluster sandwich; an inherited "edf" (from the fit) stays silent.
  if (dof_explicitly_set && dof_correction != "none" && sandwich != "cluster") {
    warning(
      "dof_correction = \"",
      dof_correction,
      "\" only applies to sandwich = \"cluster\" and is ignored for ",
      "sandwich = \"",
      sandwich,
      "\".",
      call. = FALSE
    )
  }

  is_gls_fit <- !is.null(object$pffr$hatSigma)
  if (is_gls_fit) {
    if (sandwich_missing) {
      sandwich <- "none"
    } else if (sandwich != "none") {
      warning(
        "sandwich = \"",
        sandwich,
        "\" is not supported for legacy pffr_gls fits. ",
        "Use pffr() with sandwich CIs instead. ",
        "Using sandwich = \"none\".",
        call. = FALSE
      )
      sandwich <- "none"
    }
  }

  if (
    !is.numeric(level) ||
      length(level) != 1 ||
      !is.finite(level) ||
      level <= 0 ||
      level >= 1
  ) {
    stop("'level' must be a single number in (0, 1).")
  }
  if (!is.numeric(n_sim) || length(n_sim) != 1 || n_sim < 2) {
    stop("'n_sim' must be a single integer >= 2.")
  }
  n_sim <- as.integer(n_sim)

  if (ci != "none" && !se) {
    warning("Setting se = TRUE because confidence intervals were requested.")
    se <- TRUE
  }
  if (!is.null(sim_seed) && (!is.numeric(sim_seed) || length(sim_seed) != 1)) {
    stop("'sim_seed' must be NULL or a single integer value.")
  }
  if (!is.null(sim_seed)) sim_seed <- as.integer(sim_seed)

  dots <- list(...)
  eval_grid <- dots$eval_grid %||% NULL

  # Internal, non-user-facing sandwich ablation switches (X5/X6), reachable
  # through `...` so the public coef() signature is unchanged. `b2 = FALSE`
  # drops the additive B2 term; `center_scores = TRUE` centers the per-cluster
  # score sums before the meat. Both default to the shipped behavior.
  b2 <- dots$b2 %||% TRUE
  center_scores <- dots$center_scores %||% FALSE

  # Warn if deprecated Ktt argument is passed
  if ("Ktt" %in% names(dots)) {
    warning(
      "The 'Ktt' argument is deprecated and ignored. ",
      "Use sandwich=\"cluster\" for robust standard errors.",
      call. = FALSE
    )
  }
  if (raw) {
    return(object$coefficients)
  } else {
    # Prepare info structures for helper functions
    pffr_info <- list(
      yind_name = object$pffr$yind_name,
      pcre_terms = object$pffr$pcre_terms
    )
    grid_sizes <- list(n1 = n1, n2 = n2, n3 = n3)
    object_info <- list(
      coefficients = object$coefficients,
      cmX = object$cmX,
      Vp = object$Vp
    )

    getCoefs <- function(i) {
      ## Constructs a grid over the range of the covariates
      ## and returns estimated values on this grid, with
      ## by-variables set to 1.
      ## Uses extracted helper functions for modularity.
      trm <- object$smooth[[i]]
      is_pcre <- "pcre.random.effect" %in% class(trm)

      # Check for unsupported dimensions
      if (trm$dim > 3 && !is_pcre) {
        warning(
          "can't deal with smooths with more than 3 dimensions, returning NULL for ",
          shrtlbls[names(object$smooth)[i] == unlist(object$pffr$label_map)]
        )
        return(NULL)
      }

      # Generate evaluation grid and compute predictions
      d <- resolve_eval_grid_for_term(
        eval_grid = eval_grid,
        smooth_names = names(object$smooth),
        i = i
      )
      if (is.null(d)) {
        d <- coef_make_data_grid(
          trm,
          object$model,
          pffr_info,
          grid_sizes,
          is_pcre
        )
      } else {
        d <- ensure_grid_axis_attributes(
          d = d,
          trm = trm,
          is_pcre = is_pcre,
          pffr_info = pffr_info
        )
      }
      P <- coef_get_predictions(
        trm,
        d,
        object_info,
        pffr_info,
        covmat,
        se,
        seWithMean,
        is_pcre,
        ci = ci,
        level = level,
        coef_draws = coef_draws,
        t_scale = t_scale,
        crit_mode = crit_mode,
        crit_df_const = crit_df_const,
        df_ctx = df_ctx
      )

      # Add proper labeling
      P$main <- shrtlbls[
        names(object$smooth)[i] == unlist(object$pffr$label_map)
      ]

      # Fix axis labels for ff and sff terms
      which <- match(names(object$smooth)[i], object$pffr$label_map)
      if (which %in% object$pffr$where$ff) {
        which_ff <- which(object$pffr$where$ff == which)
        P$ylab <- object$pffr$yind_name
        xlab <- deparse(
          as.call(formula(paste("~", names(object$pffr$ff)[which_ff]))[[
            2
          ]])$xind
        )
        P$xlab <- if (xlab == "NULL") "xindex" else xlab
      }
      if (which %in% object$pffr$where$sff) {
        which_sff <- which(object$pffr$where$sff == which)
        P$ylab <- object$pffr$yind_name
        xlab <- deparse(
          as.call(formula(paste("~", names(object$pffr$ff)[which_sff]))[[
            2
          ]])$xind
        )
        P$xlab <- if (xlab == "NULL") "xindex" else xlab
        P$zlab <- gsub(".mat$", "", object$pffr$ff[[which_sff]]$xname)
      }

      P
    }

    # Resolve the covariance through the single pffr accessor: $Vp/$Vc/$Ve are
    # model-based, so recomputing a sandwich here never double-applies it.
    # sandwich = "none" returns the genuinely model-based covariance; any other
    # value (or a custom `cluster`) recomputes from the model-based bread.
    covmat <- pffr_vcov(
      object,
      sandwich = sandwich,
      freq = freq,
      cluster = cluster,
      dof_correction = dof_correction,
      edf_type = edf_type,
      b2 = b2,
      center_scores = center_scores
    )

    # Pointwise critical-value reference (S3). `crit` selects the pointwise
    # reference distribution (independent of the simultaneous-band `ci_ref`):
    # "z" (Gaussian), "tG1" (t_{G-1}, the pointwise counterpart of the
    # simultaneous ci_ref = "t"), "satterthwaite" (per-point Bell-McCaffrey df),
    # or "auto" (satterthwaite for cluster/cl2 SEs at G < 150, else z). df is a
    # pointwise concept, so this only engages for ci = "pointwise"; simultaneous
    # bands keep their existing multiplier machinery.
    crit_mode <- "z"
    crit_df_const <- NA_real_
    df_ctx <- NULL
    if (ci == "pointwise") {
      df_G <- length(unique(build_cluster_id(object$pffr, cluster = cluster)))
      crit_mode <- resolve_crit_reference(crit, sandwich, df_G)
      if (crit_mode == "tG1") {
        crit_df_const <- df_G - 1
        if (!is.finite(crit_df_const) || crit_df_const < 1) {
          warning(
            "crit = \"tG1\" requires at least two independent curves or ",
            "clusters; using crit = \"z\" instead.",
            call. = FALSE
          )
          crit_mode <- "z"
        }
      } else if (crit_mode == "satterthwaite") {
        df_ctx <- pffr_df_context(object, sandwich, cluster = cluster)
        if (!isTRUE(df_ctx$ok)) {
          # No whitened score path for this family; degrade to the Gaussian
          # reference (still an honest pointwise interval from the robust SE).
          crit_mode <- "z"
          df_ctx <- NULL
        }
      }
    }

    coef_draws <- NULL
    t_scale <- NULL
    ci_ref_n_clusters <- NA_integer_
    ci_ref_df <- NA_real_
    ci_ref_used <- NA_character_
    if (ci == "simultaneous") {
      ci_cluster_id <- build_cluster_id(object$pffr, cluster = cluster)
      ci_ref_n_clusters <- length(unique(ci_cluster_id))
      ci_ref_df <- ci_ref_n_clusters - 1
      draw_df <- NULL
      ci_ref_used <- "normal"
      if (ci_ref == "t") {
        if (ci_ref_df >= 1) {
          draw_df <- ci_ref_df
          ci_ref_used <- "t"
        } else {
          warning(
            "ci_ref = \"t\" requires at least two independent curves or clusters; ",
            "using ci_ref = \"normal\" for this simultaneous band.",
            call. = FALSE
          )
        }
      }
      coef_draws <- draw_coef_perturbations(
        covmat = covmat,
        n_sim = n_sim,
        sim_seed = sim_seed,
        df = draw_df
      )
      if (is.list(coef_draws)) {
        t_scale <- coef_draws$t_scale
        coef_draws <- coef_draws$draws
      }
    }

    ret <- list()
    smind <- unlist(sapply(object$smooth, function(x) {
      seq(x$first.para, x$last.para)
    }))
    ret$pterms <- cbind(value = object$coefficients[-smind])
    if (se) ret$pterms <- cbind(ret$pterms, se = sqrt(diag(covmat)[-smind]))

    if (se && ci != "none") {
      p_se <- ret$pterms[, "se"]
      p_df <- NULL
      if (ci == "pointwise") {
        prob <- (1 + level) / 2
        if (crit_mode == "tG1") {
          p_crit <- stats::qt(prob, crit_df_const)
          p_df <- rep(crit_df_const, length(p_se))
        } else if (crit_mode == "satterthwaite") {
          # Per-parametric-coefficient df: contrasts are unit vectors e_j at the
          # non-smooth coefficient indices (same order as ret$pterms rows).
          pind <- seq_along(object$coefficients)[-smind]
          if (length(pind) > 0) {
            Xp_p <- matrix(0, length(pind), length(object$coefficients))
            Xp_p[cbind(seq_along(pind), pind)] <- 1
            p_df <- pffr_df_from_context(df_ctx, Xp_p)
          } else {
            p_df <- numeric(0)
          }
          p_crit <- ifelse(
            is.finite(p_df),
            stats::qt(prob, pmax(p_df, 1)),
            stats::qnorm(prob)
          )
        } else {
          p_crit <- stats::qnorm(prob)
          p_df <- rep(Inf, length(p_se))
        }
      } else if (
        length(p_se) == 0 ||
          (length(p_se) <= 1 && is.null(t_scale))
      ) {
        p_crit <- stats::qnorm((1 + level) / 2)
      } else {
        eps <- sqrt(.Machine$double.eps)
        valid <- is.finite(p_se) & (p_se > eps)
        if (!any(valid)) {
          p_crit <- 0
        } else {
          p_draws <- coef_draws[-smind, , drop = FALSE]
          max_stat <- apply(
            abs(p_draws[valid, , drop = FALSE] / p_se[valid]),
            2,
            max
          )
          if (!is.null(t_scale)) max_stat <- max_stat * t_scale
          p_crit <- as.numeric(stats::quantile(
            max_stat,
            probs = level,
            names = FALSE,
            type = 8,
            na.rm = TRUE
          ))
        }
      }
      p_half <- p_crit * p_se
      ret$pterms <- cbind(
        ret$pterms,
        lower = ret$pterms[, "value"] - p_half,
        upper = ret$pterms[, "value"] + p_half
      )
      if (!is.null(p_df)) {
        ret$pterms <- cbind(ret$pterms, df = p_df)
      }
    }

    shrtlbls <- object$pffr$short_labels

    ret$smterms <- lapply(1:length(object$smooth), getCoefs)
    names(ret$smterms) <- sapply(seq_along(ret$smterms), function(i) {
      ret$smterms[[i]]$main
    })
    ret$ci_meta <- list(
      type = ci,
      level = level,
      n_sim = if (ci == "simultaneous") n_sim else NA_integer_,
      sim_seed = sim_seed,
      ci_ref = if (ci == "simultaneous") ci_ref else NA_character_,
      ci_ref_used = ci_ref_used,
      ci_ref_n_clusters = ci_ref_n_clusters,
      ci_ref_df = ci_ref_df,
      crit = crit,
      crit_used = if (ci == "pointwise") crit_mode else NA_character_
    )
    return(ret)
  }
}

#' Covariance matrix for a pffr fit
#'
#' Dispatches to \code{\link[mgcv]{vcov.gam}()} on the fit's \emph{model-based}
#' covariance. Since refund now keeps \code{$Vp}/\code{$Vc}/\code{$Ve}
#' model-based on every fit (the robust covariance lives in
#' \code{object$pffr$Vsandwich}), \code{vcov(object)} returns the Bayesian
#' posterior (model-based) covariance even when the fit was created with a
#' sandwich option. With \code{sandwich = TRUE}, mgcv recomputes an
#' observation-level HC sandwich from that model-based bread. For the
#' cluster-robust estimators (or to obtain the fit-time robust covariance) use
#' \code{\link{coef.pffr}} with the \code{sandwich} argument.
#'
#' @param object a fitted \code{pffr}-object
#' @param sandwich compute an observation-level HC sandwich covariance? See
#'   \code{\link[mgcv]{vcov.gam}()}.
#' @param ... see \code{\link[mgcv]{vcov.gam}()} for options.
#'
#' @return A covariance matrix, see \code{\link[mgcv]{vcov.gam}()}.
#' @export
#' @method vcov pffr
vcov.pffr <- function(object, sandwich = FALSE, ...) {
  object <- pffr_model_based_gam(object)
  stats::vcov(object, sandwich = sandwich, ...)
}

#' Summarize the Satterthwaite pointwise-CI df of a pffr fit
#'
#' For a fit whose sandwich type is cluster-robust (`"cluster"`/`"cl2"`) and
#' whose curve/cluster count is moderate (\eqn{G < 150}, so the `crit = "auto"`
#' pointwise reference would be Satterthwaite), computes the per-point
#' Bell-McCaffrey df over a coarse coefficient grid and returns its median and
#' minimum. Returns `NULL` (so [print.summary.pffr()] prints nothing) for
#' non-cluster fits, large \eqn{G}, families without a whitened score path, or
#' any error. Uses a coarse grid to stay cheap.
#'
#' @param object A fitted pffr model.
#' @returns `NULL`, or a list with `type`, `G`, `median`, `min`, `n_points`.
#' @keywords internal
pffr_summary_df <- function(object) {
  type <- normalize_sandwich_type(
    object$pffr$sandwich_info$type %||% object$pffr$sandwich
  )
  if (!type %in% c("cluster", "cl2")) {
    return(NULL)
  }
  G <- tryCatch(
    length(unique(build_cluster_id(object$pffr))),
    error = function(e) NA_integer_
  )
  if (!is.finite(G) || G >= 150) {
    return(NULL)
  }
  # The suppressWarnings() below would otherwise silently consume the
  # once-per-session approx-score disclosure (pffr_warn_approx_score, fired via
  # pffr_df_context) before the user ever saw it: the warn-once key is set even
  # when the warning is muffled. Snapshot the approx_score_* keys and drop any
  # that appear during the suppressed call, so a later user-facing computation
  # still discloses.
  keys_before <- grep(
    "^approx_score_",
    ls(envir = .pffr_state),
    value = TRUE
  )
  df_vals <- tryCatch(
    {
      # Call coef.pffr() directly: summary.pffr() strips the "pffr" class from
      # `object` before this runs, so the generic would dispatch elsewhere.
      co <- suppressMessages(suppressWarnings(coef.pffr(
        object,
        se = TRUE,
        ci = "pointwise",
        crit = "satterthwaite",
        sandwich = type,
        n1 = 20,
        n2 = 12,
        n3 = 8
      )))
      vals <- unlist(lapply(co$smterms, function(tm) tm$coef$df))
      if (!is.null(co$pterms) && "df" %in% colnames(co$pterms)) {
        vals <- c(vals, co$pterms[, "df"])
      }
      vals[is.finite(vals)]
    },
    error = function(e) numeric(0)
  )
  keys_new <- setdiff(
    grep("^approx_score_", ls(envir = .pffr_state), value = TRUE),
    keys_before
  )
  if (length(keys_new) > 0) {
    rm(list = keys_new, envir = .pffr_state)
  }
  if (length(df_vals) == 0) {
    return(NULL)
  }
  list(
    type = type,
    G = G,
    median = stats::median(df_vals),
    min = min(df_vals),
    n_points = length(df_vals)
  )
}

#' Summary for a pffr fit
#'
#' Take a fitted \code{pffr}-object and produce summaries from it.
#' See \code{\link[mgcv]{summary.gam}()} for details.
#'
#' @param object a fitted \code{pffr}-object
#' @param ... see \code{\link[mgcv]{summary.gam}()} for options.
#'
#' @return A list with summary information, see \code{\link[mgcv]{summary.gam}()}
#' @export
#' @method summary pffr
#' @importFrom mgcv summary.gam
#' @author Fabian Scheipl, adapted from \code{\link[mgcv]{summary.gam}()} by Simon Wood, Henric Nilsson
summary.pffr <- function(object, ...) {
  call <- match.call()
  call[[1]] <- mgcv::summary.gam
  ar1rho <- object$AR1.rho
  ## drop "pffr" class and replace <object> with changed value s.t. method dispatch works without glitches
  ## if we don't do this, summary.gam will call predict on the object if n>3000 & freq==TRUE
  ## and this predict-call gets dispatched to predict.pffr which dispatches back
  ## to predict.gam. Somewhere along the way an index variable get's lost and
  ## shit breaks down.
  class(object) <- class(object)[!(class(object) %in% "pffr")]
  call$object <- as.name("object")
  ret <- eval(call)

  ret$formula <- object$pffr$formula

  # Use pre-computed short labels
  shrtlbls <- object$pffr$short_labels

  if (!is.null(ret$s.table)) {
    rownames(ret$s.table) <- vapply(
      rownames(ret$s.table),
      \(x) {
        # Direct lookup in shortlabels
        if (x %in% names(shrtlbls)) {
          shrtlbls[[x]]
        } else {
          # Fallback: try partial match against labelmap for backwards compat
          idx <- pmatch(x, unlist(object$pffr$label_map))
          if (!is.na(idx)) shrtlbls[[idx]] else x
        }
      },
      character(1)
    )
  }

  # Handle parametric effects for multi-linear-predictor families (e.g., gaulss)
  # These have names like "(Intercept).1", "grpB.1" in p.table
  # Only apply this transformation for families with multiple linear predictors
  is_multi_lp <- !is.null(object$family$nlp) && object$family$nlp > 1
  if (is_multi_lp && !is.null(ret$p.table)) {
    # Use log(SD) label only for gaulss, generic lpN label for other families
    is_gaulss <- identical(object$family$family, "gaulss")
    rownames(ret$p.table) <- vapply(
      rownames(ret$p.table),
      \(x) {
        # Check for .N suffix indicating additional linear predictor (N > 0)
        if (grepl("\\.([0-9]+)$", x)) {
          # Extract base name and suffix
          base_name <- sub("\\.([0-9]+)$", "", x)
          lp_num <- sub(".*\\.([0-9]+)$", "\\1", x)
          if (is_gaulss) {
            paste0("log(SD): ", base_name)
          } else {
            paste0("lp", as.integer(lp_num) + 1, ": ", base_name)
          }
        } else {
          x
        }
      },
      character(1)
    )
  }
  class(ret) <- c("summary.pffr", class(ret))
  if (!object$pffr$is_sparse) {
    ret$n <- paste(
      ret$n,
      " (",
      object$pffr$nobs,
      " x ",
      object$pffr$nyindex,
      ")",
      sep = ""
    )
  } else {
    ret$n <- paste(ret$n, " (in ", object$pffr$nobs, " curves)", sep = "")
  }
  ret$sandwich <- object$pffr$sandwich
  ret$satterthwaite_df <- pffr_summary_df(object)
  # Descriptive within-curve dependence flag (S5); never fatal to summary().
  ret$dependence <- tryCatch(
    pffr_dependence_check(object),
    error = function(e) NULL
  )
  if (!is.null(ar1rho)) {
    ret$AR1.rho <- ar1rho
  }
  return(ret)
}

#' Print method for summary of a pffr fit
#'
#' Pretty printing for a \code{summary.pffr}-object.
#' See \code{\link[mgcv]{print.summary.gam}()} for details.
#'
#' @param x a fitted \code{pffr}-object
#' @param digits controls number of digits printed in output.
#' @param signif.stars Should significance stars be printed alongside output?
#' @param ... not used
#'
#' @return A \code{\link{summary.pffr}} object
#' @method print summary.pffr
#' @importFrom stats printCoefmat
#' @export
#' @author Fabian Scheipl, adapted from \code{\link[mgcv]{print.summary.gam}()} by Simon Wood, Henric Nilsson
print.summary.pffr <- function(
  x,
  digits = max(3, getOption("digits") - 3),
  signif.stars = getOption("show.signif.stars"),
  ...
) {
  # mostly identical to print.summary.gam
  print(x$family)
  cat("Formula:\n")
  print(x$formula)
  if (length(x$p.coeff) > 0) {
    cat("\nConstant coefficients:\n")
    printCoefmat(
      x$p.table,
      digits = digits,
      signif.stars = signif.stars,
      na.print = "NA",
      ...
    )
  }
  cat("\n")
  if (!is.null(x$AR1.rho) && is.finite(x$AR1.rho) && abs(x$AR1.rho) > 0) {
    cat(
      "AR(1) residual correlation (rho):",
      formatC(x$AR1.rho, digits = digits, format = "fg"),
      "\n\n"
    )
  }
  if (x$m > 0) {
    cat("Smooth terms & functional coefficients:\n")
    printCoefmat(
      x$s.table,
      digits = digits,
      signif.stars = signif.stars,
      has.Pvalue = TRUE,
      na.print = "NA",
      cs.ind = 1,
      ...
    )
  }
  cat("\nR-sq.(adj) = ", formatC(x$r.sq, digits = 3, width = 5))
  if (length(x$dev.expl) > 0)
    cat(
      "   Deviance explained = ",
      formatC(x$dev.expl * 100, digits = 3, width = 4),
      "%\n",
      sep = ""
    )

  if (!is.null(x$method) && !(x$method %in% c("PQL", "lme.ML", "lme.REML")))
    cat(x$method, " score = ", formatC(x$sp.criterion, digits = 5), sep = "")

  cat(
    "  Scale est. = ",
    formatC(x$scale, digits = 5, width = 8, flag = "-"),
    "  n = ",
    x$n,
    "\n",
    sep = ""
  )
  sw <- normalize_sandwich_type(x$sandwich)
  if (sw != "none") {
    cat(
      "Model fitted with sandwich = \"",
      sw,
      "\"; this summary reports model-based uncertainty. ",
      "Use coef(., sandwich = \"",
      sw,
      "\") for robust intervals.\n",
      sep = ""
    )
  }
  if (!is.null(x$satterthwaite_df)) {
    st <- x$satterthwaite_df
    cat(sprintf(
      paste0(
        "Satterthwaite df for %s pointwise CIs (crit = \"auto\"): ",
        "median %s, min %s (G = %d).\n"
      ),
      st$type,
      formatC(st$median, digits = digits, format = "fg"),
      formatC(st$min, digits = digits, format = "fg"),
      st$G
    ))
  }
  if (!is.null(x$dependence)) {
    print_dependence_line(x$dependence, digits = digits)
  }
  invisible(x)
}

#' Within-curve dependence diagnostic for a pffr fit
#'
#' @description
#' A quick, DESCRIPTIVE flag for how strongly the working residuals are
#' correlated \emph{within} each functional response curve. It is \strong{not a
#' test and not an estimator}: it exists only to tell users which inference
#' regime they are in, since the model-based vs. robust interval trade-off
#' hinges on within-curve dependence. When dependence is weak, model-based
#' intervals are broadly valid and the robust (cluster/CL2) sandwich costs a few
#' points of coverage; when it is strong, model-based intervals are
#' anti-conservative and the robust path is needed.
#'
#' @details
#' For each curve the working residuals are ordered by the functional index
#' \eqn{t} (irregular grids are sorted; the sandwich's own curve alignment is
#' reused) and their lag-1..\code{max(lags)} sample autocorrelations are
#' computed. Reported summaries:
#' \itemize{
#'   \item \code{rho1_mean}, \code{rho1_iqr}: across-curve mean and IQR of the
#'     per-curve lag-1 autocorrelation.
#'   \item \code{Dbar}: mean number of residual points per curve.
#'   \item \code{DE}: a crude plug-in design effect
#'     \eqn{DE = 1 + (\bar D - 1)\,\max(\bar{\bar\rho}, 0)}, where
#'     \eqn{\bar{\bar\rho}} is the across-curve mean of each curve's mean SIGNED
#'     autocorrelation over \code{lags}. This deviates from the brief's literal
#'     "mean ABSOLUTE autocorrelation": \eqn{|\rho|} has a positive sampling-noise
#'     floor \eqn{\sim\sqrt{2/(\pi D)}} per lag that grows with the grid, so the
#'     absolute version never approaches 1 under independence; averaging the
#'     signed autocorrelations across curves cancels that mean-zero noise and a
#'     single clip at 0 handles alternating dependence. It is a deliberately
#'     crude flag (a working-independence effective-sample-size heuristic), NOT a
#'     variance estimate.
#'   \item \code{N_eff}: implied effective sample size \eqn{N / DE}, shown next
#'     to the number of curves \eqn{G}.
#' }
#' The advisory follows a two-regime rule: \code{DE < 1.5} ("dependence looks
#' weak") vs. \code{DE >= 1.5} ("dependence detected"). Caveat: because the
#' across-curve average uses SIGNED autocorrelations, opposite-sign per-curve
#' autocorrelations can cancel each other in \code{DE} --- consistent with its
#' role as a crude descriptive flag, not an estimator. For a single-curve fit
#' (\eqn{G < 2}) the design effect is meaningless; \code{DE}/\code{N_eff} are
#' \code{NA} and the advisory says so.
#'
#' @param fit A fitted \code{\link{pffr}} model.
#' @param lags Integer lags whose per-curve SIGNED autocorrelations are averaged
#'   (then averaged across curves and clipped at 0) in \code{DE} (default
#'   \code{1:3}). The lag-1 summaries always use lag 1.
#' @returns An object of class \code{"pffr_dependence_check"}: a list with
#'   \code{rho1_mean}, \code{rho1_iqr}, per-curve lag-1 \code{rho1}, per-curve
#'   over-lag mean \code{rho_avg}, the across-curve mean \code{rhobar_avg},
#'   \code{Dbar}, \code{DE}, \code{G}, \code{N}, \code{N_eff}, \code{lags},
#'   \code{regime} (\code{"weak"}/\code{"detected"}/\code{"undetermined"}),
#'   \code{advisory}, and \code{n_curves_used}.
#' @seealso \code{\link{pffr}}, \code{\link{coef.pffr}}
#' @export
#' @author Fabian Scheipl
pffr_dependence_check <- function(fit, lags = 1:3) {
  if (is.null(fit$pffr)) {
    stop("`fit` must be a fitted pffr model.", call. = FALSE)
  }
  lags <- sort(unique(as.integer(lags)))
  if (length(lags) < 1L || any(!is.finite(lags)) || any(lags < 1L)) {
    stop("`lags` must be positive integers.", call. = FALSE)
  }
  meta <- fit$pffr

  # Working residuals as a plain vector (bypass residuals.pffr dispatch so this
  # also works on the class-stripped object summary.pffr passes internally).
  gamobj <- fit
  class(gamobj) <- setdiff(class(gamobj), "pffr")
  wr <- as.numeric(stats::residuals(gamobj, type = "working"))

  # Map each residual to its curve and functional index, reusing the sandwich's
  # curve alignment (build_cluster_id) so the ordering is guaranteed consistent.
  cid <- build_cluster_id(meta)
  if (isTRUE(meta$is_sparse)) {
    tval <- meta$ydata$.index
    # Sparse fits record no missing_indices, but mgcv silently drops rows with
    # NA response (na.omit) while ydata keeps them; filter those rows so
    # cid/tval align with the fitted residuals (verified empirically:
    # fit$y == ydata$.value[!is.na(.value)] on a sparse fit with NA .value).
    na_y <- is.na(meta$ydata$.value)
    if (any(na_y)) {
      cid <- cid[!na_y]
      tval <- tval[!na_y]
    }
  } else {
    tval <- rep(meta$yind, times = meta$nobs)
    if (!is.null(meta$missing_indices)) {
      tval <- tval[-meta$missing_indices]
    }
  }
  if (length(wr) != length(cid) || length(tval) != length(cid)) {
    stop(
      "could not align working residuals to curves (length mismatch); ",
      "the dependence diagnostic is unavailable for this fit.",
      call. = FALSE
    )
  }

  # per-curve residual series, sorted by t within curve
  ord <- order(cid, tval)
  series <- split(wr[ord], cid[ord])

  Lmax <- max(lags)
  per_curve <- lapply(series, function(x) {
    x <- x[is.finite(x)]
    D <- length(x)
    ac <- rep(NA_real_, Lmax)
    if (D >= 2L && stats::sd(x) > 0) {
      k <- min(Lmax, D - 1L)
      acf_vals <- tryCatch(
        stats::acf(x, lag.max = k, plot = FALSE, demean = TRUE)$acf[-1L],
        error = function(e) rep(NA_real_, k)
      )
      ac[seq_len(k)] <- acf_vals
    }
    list(D = D, acf = ac)
  })

  D_vec <- vapply(per_curve, function(z) z$D, numeric(1))
  rho1 <- vapply(per_curve, function(z) z$acf[1L], numeric(1)) # lag-1 (signed)
  # Per-curve mean SIGNED autocorrelation over `lags`. Deviation from the
  # brief's literal "mean ABSOLUTE autocorrelation": |rho| has expectation
  # ~sqrt(2/(pi D)) per lag under independence, a positive noise floor that
  # grows with the grid so the absolute-value DE never approaches 1 for iid
  # data (and its pmax(., 0) is vacuous on non-negative values). Averaging the
  # SIGNED autocorrelations ACROSS curves cancels that mean-zero noise, and a
  # single pmax(., 0) at the end clips net negative (alternating) dependence to
  # DE = 1. It is a deliberately crude regime flag, not a variance estimate.
  rho_avg <- vapply(
    per_curve,
    function(z) {
      a <- z$acf[lags]
      if (all(is.na(a))) NA_real_ else mean(a, na.rm = TRUE)
    },
    numeric(1)
  )

  Dbar <- mean(D_vec)
  N <- sum(D_vec)
  G <- length(series)
  rhobar_avg <- mean(rho_avg, na.rm = TRUE)
  if (!is.finite(rhobar_avg)) rhobar_avg <- 0

  if (G < 2) {
    # Single curve: an across-curve design effect and cluster-robust intervals
    # are both meaningless here; NA the derived quantities and say so instead
    # of printing a misleading "dependence looks weak" line.
    DE <- NA_real_
    N_eff <- NA_real_
    regime <- "undetermined"
    advisory <- paste0(
      "only one curve: the dependence diagnostic and cluster-robust ",
      "intervals are not meaningful"
    )
  } else {
    DE <- 1 + (Dbar - 1) * max(rhobar_avg, 0)
    N_eff <- if (is.finite(DE) && DE > 0) N / DE else NA_real_
    regime <- if (is.finite(DE) && DE >= 1.5) "detected" else "weak"
    advisory <- if (regime == "detected") {
      paste0(
        "within-curve dependence detected; model-based intervals would be ",
        # TODO(S7): point at ?pffr_inference once that help topic/vignette
        # ships with the release task; it does not exist yet.
        "anti-conservative (see ?pffr_dependence_check)"
      )
    } else {
      paste0(
        "within-curve dependence looks weak; model-based and robust ",
        "intervals should broadly agree (robust costs a few points of ",
        "coverage here)"
      )
    }
  }

  structure(
    list(
      rho1_mean = mean(rho1, na.rm = TRUE),
      rho1_iqr = stats::IQR(rho1, na.rm = TRUE),
      rho1 = rho1,
      rho_avg = rho_avg,
      rhobar_avg = rhobar_avg,
      Dbar = Dbar,
      DE = DE,
      G = G,
      N = N,
      N_eff = N_eff,
      lags = lags,
      regime = regime,
      advisory = advisory,
      n_curves_used = sum(!is.na(rho1))
    ),
    class = "pffr_dependence_check"
  )
}

# Shared one-line + advisory printer, used by both print.pffr_dependence_check
# and print.summary.pffr.
print_dependence_line <- function(x, digits = 3) {
  if (!is.finite(x$DE)) {
    # Single-curve (or otherwise undetermined) case: only the advisory.
    cat(
      "Within-curve dependence (descriptive flag): ",
      x$advisory,
      ".\n",
      sep = ""
    )
    return(invisible(x))
  }
  fmt <- function(v) formatC(v, digits = digits, format = "fg")
  cat(sprintf(
    paste0(
      "Within-curve dependence (descriptive flag): mean rho1 = %s ",
      "(IQR %s); design effect DE = %s, implied N_eff = %s vs G = %d.\n"
    ),
    fmt(x$rho1_mean),
    fmt(x$rho1_iqr),
    fmt(x$DE),
    formatC(x$N_eff, digits = 0, format = "f"),
    x$G
  ))
  cat("  ", x$advisory, "\n", sep = "")
  invisible(x)
}

#' Print method for a pffr within-curve dependence diagnostic
#'
#' @param x A \code{"pffr_dependence_check"} object from
#'   \code{\link{pffr_dependence_check}}.
#' @param digits Number of significant digits for the printed summaries.
#' @param ... Not used.
#' @returns \code{x}, invisibly.
#' @method print pffr_dependence_check
#' @export
print.pffr_dependence_check <- function(
  x,
  digits = max(3, getOption("digits") - 3),
  ...
) {
  print_dependence_line(x, digits = digits)
}

#' QQ plots for pffr model residuals
#'
#' This is simply a wrapper for \code{\link[mgcv]{qq.gam}()}.
#'
#' @param object a fitted \code{\link{pffr}}-object
#' @inheritParams mgcv::qq.gam
#' @return None, called for its side effect of producing QQ plots.
#' @importFrom mgcv qq.gam
#' @export
pffr_qq <- function(
  object,
  rep = 0,
  level = 0.9,
  s.rep = 10,
  type = c("deviance", "pearson", "response"),
  pch = ".",
  rl.col = 2,
  rep.col = "gray80",
  ...
) {
  if (!inherits(object, "pffr")) stop("`object' is not of class \"pffr\"")
  call <- match.call()
  # drop pffr-class so only gam-methods are used on object
  class(object) <- class(object)[-1]
  call$object <- object
  call[[1]] <- mgcv::qq.gam
  eval(call)
}


#' QQ plots for pffr model residuals (deprecated)
#'
#' @description
#' **Deprecated**
#'
#' `qq.pffr()` was renamed to [pffr_qq()] for consistency with the
#' package naming conventions.
#'
#' @inheritParams pffr_qq
#' @return None, called for its side effect of producing QQ plots.
#' @export
#' @keywords internal
qq.pffr <- function(
  object,
  rep = 0,
  level = 0.9,
  s.rep = 10,
  type = c("deviance", "pearson", "response"),
  pch = ".",
  rl.col = 2,
  rep.col = "gray80",
  ...
) {
  .Deprecated("pffr_qq")
  pffr_qq(
    object = object,
    rep = rep,
    level = level,
    s.rep = s.rep,
    type = type,
    pch = pch,
    rl.col = rl.col,
    rep.col = rep.col,
    ...
  )
}

#' Some diagnostics for a fitted pffr model
#'
#' This is simply a wrapper for \code{\link[mgcv]{gam.check}()}.
#'
#' @param b a fitted \code{\link{pffr}}-object
#' @inheritParams mgcv::gam.check
#' @param rep passed to \code{\link[mgcv]{qq.gam}} when \code{old.style} is \code{FALSE}.
#' @param level passed to \code{\link[mgcv]{qq.gam}} when \code{old.style} is \code{FALSE}.
#' @param rl.col passed to \code{\link[mgcv]{qq.gam}} when \code{old.style} is \code{FALSE}.
#' @param rep.col passed to \code{\link[mgcv]{qq.gam}} when \code{old.style} is \code{FALSE}.
#' @return None, called for its side effect of producing diagnostic plots and
#'   printing basis dimension checks.
#' @export
pffr_check <- function(
  b,
  old.style = FALSE,
  type = c("deviance", "pearson", "response"),
  k.sample = 5000,
  k.rep = 200,
  rep = 0,
  level = 0.9,
  rl.col = 2,
  rep.col = "gray80",
  ...
) {
  if (!inherits(b, "pffr")) stop("`object' is not of class \"pffr\"")
  call <- match.call()
  # drop pffr-class so only gam-methods are used on b
  class(b) <- class(b)[-1]
  call$b <- b
  call[[1]] <- mgcv::gam.check
  eval(call)
}


#' Some diagnostics for a fitted pffr model (deprecated)
#'
#' @description
#' **Deprecated**
#'
#' `pffr.check()` was renamed to [pffr_check()] for consistency with the
#' package naming conventions.
#'
#' @inheritParams pffr_check
#' @return None, called for its side effect of producing diagnostic plots.
#' @export
#' @keywords internal
pffr.check <- function(
  b,
  old.style = FALSE,
  type = c("deviance", "pearson", "response"),
  k.sample = 5000,
  k.rep = 200,
  rep = 0,
  level = 0.9,
  rl.col = 2,
  rep.col = "gray80",
  ...
) {
  .Deprecated("pffr_check")
  pffr_check(
    b = b,
    old.style = old.style,
    type = type,
    k.sample = k.sample,
    k.rep = k.rep,
    rep = rep,
    level = level,
    rl.col = rl.col,
    rep.col = rep.col,
    ...
  )
}
