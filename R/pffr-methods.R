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
        !(paste0(object$pffr$yindname, ".vec") %in% names(newdata))
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
      #        covnames <- unique(covnames[covnames != paste(object$pffr$yindname, ".vec", sep="")])
      #        stopifnot(all(covnames %in% names(newdata)))

      #get newdata into the shape expected by predict gam:
      gamdata <- list()
      #y-index
      gamdata[[paste(object$pffr$yindname, ".vec", sep = "")]] <- rep(
        object$pffr$yind,
        times = nobs
      )

      # which covariates occur in which terms?
      varmap <- sapply(
        names(object$pffr$labelmap),
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
#' @return This function only generates plots.
#' @method plot pffr
#' @importFrom mgcv plot.gam
#' @author Fabian Scheipl
plot.pffr <- function(x, ...) {
  call <- match.call()
  call[[1]] <- mgcv::plot.gam
  #drop "pffr" class and replace <object> with changed value s.t. method dispatch works without glitches
  class(x) <- class(x)[-1]
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
#' @param pffr_info List with pffr metadata: yindname, pcreterms.
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
      sapply(pffr_info$pcreterms, `[[`, "idname") == trm$term[1]
    )
    pcreterm <- pffr_info$pcreterms[[which_pcre]]
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
    colnames(d)[1:2] <- c(trm$term[1], paste0(pffr_info$yindname, ".vec"))
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

#' Compute predictions for coefficient extraction
#'
#' Evaluates smooth term on grid and computes coefficients and standard errors.
#'
#' @param trm Smooth term object.
#' @param data_grid Data frame from coef_make_data_grid.
#' @param object_info List with: coefficients, cmX, Vp.
#' @param pffr_info List with: yindname.
#' @param covmat Covariance matrix for SE computation.
#' @param se Logical, compute standard errors?
#' @param seWithMean Logical, include mean uncertainty?
#' @param is_pcre Logical, is this a pcre term?
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
  is_pcre
) {
  X <- PredictMat(trm, data_grid)

  # For pcre terms, temporarily adjust term for axis setup
  if (is_pcre) {
    trm$dim <- 2
    trm$term[2] <- paste0(pffr_info$yindname, ".vec")
  }

  # Build result structure based on dimensionality
  P <- build_coef_axes(trm, data_grid)

  # Compute predicted values
  trmind <- trm$first.para:trm$last.para
  P$value <- X %*% object_info$coefficients[trmind]
  P$coef <- cbind(data_grid, value = P$value)

  # Compute standard errors if requested
  if (se) {
    P$se <- compute_coef_se(
      X,
      trmind,
      trm,
      object_info,
      covmat,
      seWithMean
    )
    P$coef <- cbind(P$coef, se = P$se)
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

#' Compute standard errors for coefficient extraction
#'
#' @param X Prediction matrix.
#' @param trmind Index vector for term parameters.
#' @param trm Smooth term object (for nCons and meanL1).
#' @param object_info List with cmX, Vp.
#' @param covmat Covariance matrix.
#' @param seWithMean Logical, use seWithMean approach?
#' @returns Numeric vector of standard errors.
#' @keywords internal
compute_coef_se <- function(X, trmind, trm, object_info, covmat, seWithMean) {
  if (seWithMean && attr(trm, "nCons") > 0) {
    cat("using seWithMean for ", trm$label, ".\n")
    X1 <- matrix(object_info$cmX, nrow(X), ncol(object_info$Vp), byrow = TRUE)
    meanL1 <- trm$meanL1
    if (!is.null(meanL1)) X1 <- X1 / meanL1
    X1[, trmind] <- X
    sqrt(rowSums((X1 %*% covmat) * X1))
  } else {
    sqrt(rowSums((X %*% covmat[trmind, trmind]) * X))
  }
}


#' Get estimated coefficients from a pffr fit
#'
#' Returns estimated coefficient functions/surfaces \eqn{\beta(t), \beta(s,t)}
#' and estimated smooth effects \eqn{f(z), f(x,z)} or \eqn{f(x, z, t)} and their point-wise estimated standard errors.
#' Not implemented for smooths in more than 3 dimensions.
#'
#' The \code{seWithMean}-option corresponds to the \code{"iterms"}-option in \code{\link[mgcv]{predict.gam}}.
#' The \code{sandwich}-option uses \code{\link[mgcv]{vcov.gam}} with \code{sandwich=TRUE} to compute
#' robust standard errors. If the model was fitted with \code{sandwich=TRUE} in \code{\link{pffr}},
#' the pre-computed sandwich covariance matrices are used directly.
#'
#'
#' @param object a fitted \code{pffr}-object
#' @param raw logical, defaults to FALSE. If TRUE, the function simply returns \code{object$coefficients}
#' @param se logical, defaults to TRUE. Return estimated standard error of the estimates?
#' @param freq logical, defaults to FALSE. If FALSE, use posterior variance \code{object$Vp} for variability estimates,
#'  else use \code{object$Ve}. See \code{\link[mgcv]{gamObject}}
#' @param sandwich logical, defaults to FALSE. Use sandwich-corrected covariance for standard errors.
#'   Uses \code{\link[mgcv]{vcov.gam}} with \code{sandwich=TRUE}. If the model was fitted with
#'   \code{sandwich=TRUE}, the pre-computed covariance matrices are used.
#' @param seWithMean logical, defaults to TRUE. Include uncertainty about the intercept/overall mean in  standard errors returned for smooth components?
#' @param n1 see below
#' @param n2 see below
#' @param n3 \code{n1, n2, n3} give the number of gridpoints for 1-/2-/3-dimensional smooth terms
#' used in the marginal equidistant grids over the range of the covariates at which the estimated effects are evaluated.
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
  sandwich = FALSE,
  seWithMean = TRUE,
  n1 = 100,
  n2 = 40,
  n3 = 20,
  ...
) {
  # Warn if deprecated Ktt argument is passed
  if ("Ktt" %in% names(list(...))) {
    warning(
      "The 'Ktt' argument is deprecated and ignored. ",
      "Use sandwich=TRUE for robust standard errors via mgcv::vcov.gam().",
      call. = FALSE
    )
  }
  if (raw) {
    return(object$coefficients)
  } else {
    # Prepare info structures for helper functions
    pffr_info <- list(
      yindname = object$pffr$yindname,
      pcreterms = object$pffr$pcreterms
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
          shrtlbls[names(object$smooth)[i] == unlist(object$pffr$labelmap)]
        )
        return(NULL)
      }

      # Generate evaluation grid and compute predictions
      d <- coef_make_data_grid(
        trm,
        object$model,
        pffr_info,
        grid_sizes,
        is_pcre
      )
      P <- coef_get_predictions(
        trm,
        d,
        object_info,
        pffr_info,
        covmat,
        se,
        seWithMean,
        is_pcre
      )

      # Add proper labeling
      P$main <- shrtlbls[
        names(object$smooth)[i] == unlist(object$pffr$labelmap)
      ]

      # Fix axis labels for ff and sff terms
      which <- match(names(object$smooth)[i], object$pffr$labelmap)
      if (which %in% object$pffr$where$ff) {
        which_ff <- which(object$pffr$where$ff == which)
        P$ylab <- object$pffr$yindname
        xlab <- deparse(
          as.call(formula(paste("~", names(object$pffr$ff)[which_ff]))[[
            2
          ]])$xind
        )
        P$xlab <- if (xlab == "NULL") "xindex" else xlab
      }
      if (which %in% object$pffr$where$sff) {
        which_sff <- which(object$pffr$where$sff == which)
        P$ylab <- object$pffr$yindname
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

    if (sandwich) {
      # Use sandwich-corrected covariance matrix
      # If model was fitted with sandwich=TRUE, use stored matrices
      # Otherwise, compute sandwich correction now
      model_has_sandwich <- isTRUE(object$pffr$sandwich)
      if (model_has_sandwich) {
        covmat <- if (freq) object$Ve else object$Vp
      } else {
        # Strip pffr class to avoid predict.pffr warnings during vcov
        object_stripped <- object
        class(object_stripped) <- setdiff(class(object_stripped), "pffr")
        covmat <- stats::vcov(object_stripped, sandwich = TRUE, freq = freq)
      }
    } else {
      covmat <- if (freq) object$Ve else object$Vp
    }

    ret <- list()
    smind <- unlist(sapply(object$smooth, function(x) {
      seq(x$first.para, x$last.para)
    }))
    ret$pterms <- cbind(value = object$coefficients[-smind])
    if (se) ret$pterms <- cbind(ret$pterms, se = sqrt(diag(covmat)[-smind]))

    shrtlbls <- object$pffr$shortlabels

    ret$smterms <- lapply(1:length(object$smooth), getCoefs)
    names(ret$smterms) <- sapply(seq_along(ret$smterms), function(i) {
      ret$smterms[[i]]$main
    })
    return(ret)
  }
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
  shrtlbls <- object$pffr$shortlabels

  if (!is.null(ret$s.table)) {
    rownames(ret$s.table) <- vapply(
      rownames(ret$s.table),
      \(x) {
        # Direct lookup in shortlabels
        if (x %in% names(shrtlbls)) {
          shrtlbls[[x]]
        } else {
          # Fallback: try partial match against labelmap for backwards compat
          idx <- pmatch(x, unlist(object$pffr$labelmap))
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
  ret$sandwich <- isTRUE(object$pffr$sandwich)
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
  if (isTRUE(x$sandwich)) {
    cat("Sandwich correction applied to covariance matrices.\n")
  }
  invisible(x)
}

#' QQ plots for pffr model residuals
#'
#' This is simply a wrapper for \code{\link[mgcv]{qq.gam}()}.
#'
#' @param object a fitted \code{\link{pffr}}-object
#' @inheritParams mgcv::qq.gam
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
