#' Extract coefficient function estimate from a fitted fgam model
#' 
#' This function may be used to extract a functional coefficient from
#' a fitted \code{fgam} model, at user-specified evaluation points. It works
#' for \code{lf}, \code{af}, and \code{lf.vd} terms. This function
#' will be depricated once we write a proper \code{coef.fgam} method.
#' Can also perform a linear combination of multiple smooth terms,
#' which is useful for \code{lf.vd} terms with parametric interactions.
#' 
#' @param mod A fitted \code{fgam} model.
#' @param data A data frame containing the locations at which the term will
#'    be evaluated. See Details.
#' @param sm.id The index of the smooth term in the model. Only smooth
#'    terms are given an index number. May be supplied as a vector; see
#'    Details.
#' @param var Should the pointwise variance also be returned?
#' @param lincom Matrix of coefficients to use for a linear combination
#'    across multiple terms. Should have the same number of columns as
#'    \code{length(sm.id)}.This feature assumes all bases for the smooth
#'    terms defined by sm.id are identical. See Details.
#' 
#' @details The index indicated by \code{sm.id} corresponds to the index
#' within the \code{mod$smooth} object. There will be one entry in this
#' list for each smooth term. Use \code{names(mod$smooth)} to see a label
#' for each smooth term in the model, which should allow the user to
#' determine the appropriate \code{sm.id}.
#' 
#' The labels produced by \code{names(mod$smooth)} also contains the name
#' of all variables that are required to be included in the \code{data}
#' argument. If a "by" matrix (which will have a variable name beginning with
#' the letter "L") is not supplied, it will be assumed to be 1 (no weights),
#' which is usually desired. The other variables that are supplied will include
#' the indices of the appropriate evaluation points. For univariate
#' smooths (e.g., \code{lf} terms and \code{lf.vd} terms with parametric
#' interactions), the evaluation points should be supplied as vector columns
#' in the data frame. For bivariate smooths (e.g., \code{af} terms and fully
#' nonparametric \code{lf.vd} terms), the evaluation points should be supplied
#' as matrix-columns in the data frame, with the same dimension for each one.
#' Note that matrices can be included as columns in the data frame by using
#' the \code{I()} ("AsIs") function, or by simply assigning them to an existing
#' dataframe (\code{df$matcolumn <- myMatrix}).
#' 
#' If \code{sm.id} is supplied as a vector with \code{lincom} equal to
#' \code{NULL}, then multiple smooth terms will be returned. If
#' \code{lincom} is supplied, its rows will contain the weights for a linear
#' combination of multiple estimates. e.g., to perform the linear combination
#' \eqn{\beta_1(t) + T_i \beta_2(t) + T_i^2 \beta_3(t)}, \code{lincom} will be
#' a matrix of 3 columns: a column of 1's, a column of \eqn{T_i} values, and
#' a column of \eqn{T_i^2} values. \code{sm.id} will be a vector corresponding
#' to the indecies of \eqn{\beta_1(t)}, \eqn{\beta_2(t)}, and \eqn{\beta_3(t)}
#' in the model, and \code{data} will be a data frame with vector columns for
#' each variable. See examples.
#' 
#' @return A list, matrix, or vector containing the resulting coefficient
#' function estimate(s) and/or pointwise variance estimate(s).
#' 
#' @export
#' @author Jonathan E. Gellar <jgellar1@@jhu.edu>
#' @seealso \code{\link{fgam}}, \code{\link{lf.vd}} for examples.
#' 

getEst <- function(mod, data, sm.id=1, var=FALSE, lincom=NULL) {
  if (!is.null(lincom)) {
    if (ncol(lincom)!=length(sm.id)) {
      stop("length(sm.id) and ncol(lincom) must match")
    } else if (length(sm.id)==1) {
      stop("Need at least two smooths for lincom")
    }
  }
  
  m <- length(sm.id)
  est.list=var.list <- vector(mode = "list", length = m)
  names(est.list) = names(var.list) <- names(mod$smooth)[sm.id]
  for (j in 1:m) {
    if (!is.null(lincom) | j==1) {
      # Calculate Xp matrix and get estimate
      i <- sm.id[j]
      sm.i <- mod$smooth[[i]]
      rng.i <- sm.i$first.para:sm.i$last.para
      trms.i <- sm.i$term
      
      # Make sure we have all the data
      if (all(trms.i %in% names(data))) {
        data.i <- data[,trms.i, drop=FALSE]
        term1 <- data.i[[1]]
        if (sm.i$by!="NA") {
          # Add in by variable, or create it if not supplied
          data.i[[sm.i$by]] <- if (sm.i$by %in% names(data)) {
            data[[sm.i$by]]
          } else if (is.matrix(term1)) {
            matrix(1, nrow=nrow(term1), ncol=ncol(term1))
          } else {
            1
          }
        }
        if (is.matrix(term1)) {
          # Vectorize matrices, and remember locations of NA's
          is.mat <- TRUE
          data.i <- as.data.frame(sapply(data.i, as.vector))
          idx <- apply(data.i,1,function(x) !any(is.na(x)))
          data.i <- data.i[idx,]
        } else is.mat <- FALSE
        
        Xp <- PredictMat(sm.i, data=data.i)
        bhat <- as.vector(Xp%*%(mod$coef)[rng.i])
        est.list[[j]] <- if (is.mat & is.null(lincom)) {
          # Convert back to matrix
          tmp <- matrix(nrow=nrow(term1), ncol=ncol(term1))
          tmp[idx] <- bhat
          tmp
        } else bhat
        if (var & is.null(lincom)) {
          myXprod <- function(A,B) {
            apply(A, 1, function(a) {
              a %*% B %*% a
            })
          }
          vhat <- myXprod(Xp, mod$Vp[rng.i,rng.i])
          var.list[[j]] <- if (is.mat) {
            # Convert back to matrix
            tmp <- matrix(nrow=nrow(term1), ncol=ncol(term1))
            tmp[idx] <- vhat
            tmp
          } else vhat
        }
          
      } else {
        stop("ERROR: Not all data is supplied!")
      }
    } else {
      # Re-use existing Xp matrix to get estimate (assumes exact same basis)
      sm.i <- mod$smooth[[sm.id[j]]]
      rng.i <- sm.i$first.para:sm.i$last.para
      est.list[[j]] <- as.vector(Xp%*%(mod$coef)[rng.i])
      if (var)
        var.list[[j]] <- myXprod(Xp, mod$Vp[rng.i,rng.i])
    }
  }
  
  if (!is.null(lincom)) {
    # Linear combination step
    ret <- t(apply(lincom, 1, function(x) {
      rowSums(sapply(1:length(x), function(i) {
        x[i] * est.list[[i]]
      }))
    }))
    if (var) {
      rng.mat <- sapply(mod$smooth[sm.id], function(x) {x$first.para:x$last.para})
      tmp  <- matrix(nrow=nrow(rng.mat), ncol=nrow(rng.mat))
      vHat <- t(apply(lincom, 1, function(x) {
        tmp <- matrix(0,nrow=nrow(rng.mat), ncol=nrow(rng.mat))
        for (i in 1:m)
          for (j in 1:m)
            tmp <- tmp + x[i]*x[j]*mod$Vp[rng.mat[,i],rng.mat[,j]]
        diag(Xp %*% tmp %*% t(Xp))
      }))
      ret <- list(est=ret, var=vHat)
    }
  } else {
    # Return estimates and variances
    if (length(est.list)==1) {
      # List not needed; collapse it
      ret <- if (var) {
        list(est=est.list[[1]], var=var.list[[1]])
      } else {
        est.list[[1]]
      }
    } else {
      ret <- if (var) {
        list(est=est.list, var=var.list)
      } else est.list
    }
  }
  
  ret
}