##' Functional principal components analysis by smoothed covariance
##'
##' Decomposes functional observations using functional principal components
##' analysis. A mixed model framework is used to estimate scores and obtain
##' variance estimates.
##'
##' This function computes a FPC decomposition for a set of observed curves,
##' which may be sparsely observed and/or measured with error. A mixed model
##' framework is used to estimate curve-specific scores and variances.
##'
##' FPCA via kernel smoothing of the covariance function, with the diagonal
##' treated separately, was proposed in Staniswalis and Lee (1998) and much
##' extended by Yao et al. (2005), who introduced the "PACE" method.
##' \code{fpca.sc} uses penalized splines to smooth the covariance function, as
##' developed by Di et al. (2009) and Goldsmith et al. (2013).
##'
##' The functional data must be supplied as either \itemize{ \item an \eqn{n
##' \times d} matrix \code{Y}, each row of which is one functional observation,
##' with missing values allowed; or \item a data frame \code{ydata}, with
##' columns \code{'.id'} (which curve the point belongs to, say \eqn{i}),
##' \code{'.index'} (function argument such as time point \eqn{t}), and
##' \code{'.value'} (observed function value \eqn{Y_i(t)}).}
##'
##' @param Y,ydata the user must supply either \code{Y}, a matrix of functions
##' observed on a regular grid, or a data frame \code{ydata} representing
##' irregularly observed functions. See Details.
##' @param Y.pred if desired, a matrix of functions to be approximated using
##' the FPC decomposition.
##' @param argvals function argument.
##' @param random.int If \code{TRUE}, the mean is estimated by
##' \code{\link[gamm4]{gamm4}} with random intercepts. If \code{FALSE} (the
##' default), the mean is estimated by \code{\link[mgcv]{gam}} treating all the
##' data as independent.
##' @param nbasis number of B-spline basis functions used for estimation of the
##' mean function and bivariate smoothing of the covariance surface.
##' @param pve proportion of variance explained: used to choose the number of
##' principal components.
##' @param npc prespecified value for the number of principal components (if
##' given, this overrides \code{pve}).
##' @param var \code{TRUE} or \code{FALSE} indicating whether model-based
##' estimates for the variance of FPCA expansions should be computed.
##' @param simul logical: should critical values be estimated for simultaneous
##' confidence intervals?
##' @param sim.alpha 1 - coverage probability of the simultaneous intervals.
##' @param useSymm logical, indicating whether to smooth only the upper
##' triangular part of the naive covariance (when \code{cov.est.method==2}).
##' This can save computation time for large data sets, and allows for
##' covariance surfaces that are very peaked on the diagonal.
##' @param makePD logical: should positive definiteness be enforced for the
##' covariance surface estimate?
##' @param center logical: should an estimated mean function be subtracted from
##' \code{Y}? Set to \code{FALSE} if you have already demeaned the data using
##' your favorite mean function estimate.
##' @param cov.est.method covariance estimation method. If set to \code{1}, a
##' one-step method that applies a bivariate smooth to the \eqn{y(s_1)y(s_2)}
##' values. This can be very slow. If set to \code{2} (the default), a two-step
##' method that obtains a naive covariance estimate which is then smoothed.
##' @param integration quadrature method for numerical integration; only
##' \code{"trapezoidal"} is currently supported.
##' @return \item{Yhat}{FPC approximation (projection onto leading components)
##' of \code{Y.pred} if specified, or else of \code{Y}.} \item{scores}{\eqn{n
##' \times npc} matrix of estimated FPC scores.} \item{mu}{estimated mean
##' function (or a vector of zeroes if \code{center==FALSE}).} \item{efunctions
##' }{\eqn{d \times npc} matrix of estimated eigenfunctions of the functional
##' covariance, i.e., the FPC basis functions.} \item{evalues}{estimated
##' eigenvalues of the covariance operator, i.e., variances of FPC scores.}
##' \item{npc }{number of FPCs: either the supplied \code{npc}, or the minimum
##' number of basis functions needed to explain proportion \code{pve} of the
##' variance in the observed curves.} \item{sigma2}{estimated measurement error
##' variance.} \item{diag.var}{diagonal elements of the covariance matrices for
##' each estimated curve.} \item{VarMats}{a list containing the estimated
##' covariance matrices for each curve in \code{Y}.} \item{crit.val}{estimated
##' critical values for constructing simultaneous confidence intervals.}
##' @author Jeff Goldsmith \email{jeff.goldsmith@@columbia.edu}, Sonja Greven
##' \email{sonja.greven@@stat.uni-muenchen.de}, Lan Huo
##' \email{Lan.Huo@@nyumc.org}, Lei Huang \email{huangracer@@gmail.com}, and
##' Philip Reiss \email{phil.reiss@@nyumc.org}
##' @references Di, C., Crainiceanu, C., Caffo, B., and Punjabi, N. (2009).
##' Multilevel functional principal component analysis. \emph{Annals of Applied
##' Statistics}, 3, 458--488.
##'
##' Goldsmith, J., Greven, S., and Crainiceanu, C. (2013). Corrected confidence
##' bands for functional data using principal components. \emph{Biometrics},
##' 69(1), 41--51.
##'
##' Staniswalis, J. G., and Lee, J. J. (1998). Nonparametric regression
##' analysis of longitudinal data. \emph{Journal of the American Statistical
##' Association}, 93, 1403--1418.
##'
##' Yao, F., Mueller, H.-G., and Wang, J.-L. (2005). Functional data analysis
##' for sparse longitudinal data. \emph{Journal of the American Statistical
##' Association}, 100, 577--590.
##' @examples
##' \dontrun{
##' library(ggplot2)
##' library(reshape2)
##' data(cd4)
##'
##' Fit.MM = fpca.sc(cd4, var = TRUE, simul = TRUE)
##'
##' Fit.mu = data.frame(mu = Fit.MM$mu,
##'                     d = as.numeric(colnames(cd4)))
##' Fit.basis = data.frame(phi = Fit.MM$efunctions,
##'                        d = as.numeric(colnames(cd4)))
##'
##' ## for one subject, examine curve estimate, pointwise and simultaneous itervals
##' EX = 1
##' EX.MM = data.frame(fitted = Fit.MM$Yhat[EX,],
##'            ptwise.UB = Fit.MM$Yhat[EX,] + 1.96 * sqrt(Fit.MM$diag.var[EX,]),
##'            ptwise.LB = Fit.MM$Yhat[EX,] - 1.96 * sqrt(Fit.MM$diag.var[EX,]),
##'            simul.UB = Fit.MM$Yhat[EX,] + Fit.MM$crit.val[EX] * sqrt(Fit.MM$diag.var[EX,]),
##'            simul.LB = Fit.MM$Yhat[EX,] - Fit.MM$crit.val[EX] * sqrt(Fit.MM$diag.var[EX,]),
##'            d = as.numeric(colnames(cd4)))
##'
##' ## plot data for one subject, with curve and interval estimates
##' EX.MM.m = melt(EX.MM, id = "d")
##' ggplot(EX.MM.m, aes(x = d, y = value, group = variable, color = variable, linetype = variable)) +
##'   geom_path() +
##'   scale_linetype_manual(values = c(fitted = 1, ptwise.UB = 2,
##'                         ptwise.LB = 2, simul.UB = 3, simul.LB = 3)) +
##'   scale_color_manual(values = c(fitted = 1, ptwise.UB = 2,
##'                      ptwise.LB = 2, simul.UB = 3, simul.LB = 3)) +
##'   labs(x = "Months since seroconversion", y = "Total CD4 Cell Count")
##'
##' ## plot estimated mean function
##' ggplot(Fit.mu, aes(x = d, y = mu)) + geom_path() +
##'   labs(x = "Months since seroconversion", y = "Total CD4 Cell Count")
##'
##' ## plot the first two estimated basis functions
##' Fit.basis.m = melt(Fit.basis, id = "d")
##' ggplot(subset(Fit.basis.m, variable %in% c("phi.1", "phi.2")), aes(x = d,
##' y = value, group = variable, color = variable)) + geom_path()
##'
##' ## input a dataframe instead of a matrix
##' nid <- 20
##' nobs <- sample(10:20, nid, rep=TRUE)
##' ydata <- data.frame(
##'     .id = rep(1:nid, nobs),
##'     .index = round(runif(sum(nobs), 0, 1), 3))
##' ydata$.value <- unlist(tapply(ydata$.index,
##'                               ydata$.id,
##'                               function(x)
##'                                   runif(1, -.5, .5) +
##'                                   dbeta(x, runif(1, 6, 8), runif(1, 3, 5))
##'                               )
##'                        )
##'
##' Fit.MM = fpca.sc(ydata=ydata, var = TRUE, simul = FALSE)
##'
##' }
##' @export
##' @importFrom Matrix nearPD Matrix t as.matrix
##' @importFrom mgcv gam predict.gam
##' @importFrom gamm4 gamm4
## npc=1 seems to give error
fpca.sc <- function(Y=NULL, ydata = NULL, Y.pred=NULL, argvals = NULL, random.int = FALSE,
         nbasis = 10, pve = .99, npc = NULL, var = FALSE, simul = FALSE, sim.alpha = .95,
         useSymm = FALSE, makePD = FALSE, center=TRUE, cov.est.method = 2,
         integration="trapezoidal") {

    stopifnot((!is.null(Y) && is.null(ydata))||(is.null(Y) && !is.null(ydata)))

    # if data.frame version of ydata is provided
    sparseOrNongrid <- !is.null(ydata)
    if (sparseOrNongrid) {
        stopifnot(ncol(ydata) == 3)
        stopifnot(c(".id", ".index", ".value") == colnames(ydata))
        stopifnot(is.null(argvals))
        Y = irreg2mat(ydata)
        argvals = sort(unique(ydata$.index))
    }

    if (is.null(Y.pred)) Y.pred = Y
    D = NCOL(Y)
    I = NROW(Y)
    I.pred = NROW(Y.pred)

    if (is.null(argvals)) argvals = seq(0,1,,D)

    d.vec = rep(argvals, each = I)
    id = rep(1:I, rep(D, I))

    if (center) {
        if (random.int) gam0 = gamm4(as.vector(Y) ~ s(d.vec, k = nbasis), random=~(1|id))$gam
        else gam0 = gam(as.vector(Y) ~ s(d.vec, k = nbasis))
        mu = predict(gam0, newdata = data.frame(d.vec = argvals))
        Y.tilde = Y - matrix(mu, I, D, byrow = TRUE)
    }
    else {
        Y.tilde = Y
        mu = rep(0, D)
    }

    if (cov.est.method==2) {   # smooth raw covariance estimate
        cov.sum = cov.count = cov.mean = matrix(0, D, D)
        for (i in 1:I) {
            obs.points = which(!is.na(Y[i, ]))
            cov.count[obs.points, obs.points] = cov.count[obs.points, obs.points] + 1
            cov.sum[obs.points, obs.points] = cov.sum[obs.points, obs.points] + tcrossprod(Y.tilde[i, obs.points])
        }
        G.0 = ifelse(cov.count == 0, NA, cov.sum/cov.count)
        diag.G0 = diag(G.0)
        diag(G.0) = NA
        if (!useSymm) {
            row.vec = rep(argvals, each = D)
            col.vec = rep(argvals, D)
            npc.0 = matrix(predict(gam(as.vector(G.0) ~ te(row.vec, col.vec, k = nbasis), weights =as.vector(cov.count)), newdata = data.frame(row.vec = row.vec,col.vec = col.vec)), D, D)
            npc.0 = (npc.0 + t(npc.0))/2
        }
        else {
            use <- upper.tri(G.0, diag = TRUE)
            use[2, 1] <- use[ncol(G.0), ncol(G.0) - 1] <- TRUE
            usecov.count <- cov.count
            usecov.count[2, 1] <- usecov.count[ncol(G.0), ncol(G.0) - 1] <- 0
            usecov.count <- as.vector(usecov.count)[use]
            use <- as.vector(use)
            vG.0 <- as.vector(G.0)[use]
            row.vec <- rep(argvals, each = D)[use]
            col.vec <- rep(argvals, times = D)[use]
            mCov <- gam(vG.0 ~ te(row.vec, col.vec, k = nbasis), weights = usecov.count)
            npc.0 <- matrix(NA, D, D)
            spred <- rep(argvals, each = D)[upper.tri(npc.0, diag = TRUE)]
            tpred <- rep(argvals, times = D)[upper.tri(npc.0, diag = TRUE)]
            smVCov <- predict(mCov, newdata = data.frame(row.vec = spred, col.vec = tpred))
            npc.0[upper.tri(npc.0, diag = TRUE)] <- smVCov
            npc.0[lower.tri(npc.0)] <- t(npc.0)[lower.tri(npc.0)]
        }
    }

    else if (cov.est.method==1) {  # smooth y(s1)y(s2) values to obtain covariance estimate
        row.vec = col.vec = G.0.vec = c()
        cov.sum = cov.count = cov.mean = matrix(0, D, D)
        for (i in 1:I) {
            obs.points = which(!is.na(Y[i, ]))
            temp = tcrossprod(Y.tilde[i, obs.points])
            diag(temp) = NA
            row.vec = c(row.vec, rep(argvals[obs.points], each = length(obs.points)))
            col.vec = c(col.vec, rep(argvals[obs.points], length(obs.points)))
            G.0.vec = c(G.0.vec, as.vector(temp))
            # still need G.O raw to calculate to get the raw to get the diagonal
            cov.count[obs.points, obs.points] = cov.count[obs.points, obs.points] + 1
            cov.sum[obs.points, obs.points] = cov.sum[obs.points, obs.points] + tcrossprod(Y.tilde[i, obs.points])
        }
        row.vec.pred = rep(argvals, each = D)
        col.vec.pred = rep(argvals, D)
        npc.0 = matrix(predict(gam(G.0.vec ~ te(row.vec, col.vec, k = nbasis)), newdata = data.frame(row.vec = row.vec.pred, col.vec = col.vec.pred)), D, D)
        npc.0 = (npc.0 + t(npc.0))/2
        G.0 = ifelse(cov.count == 0, NA, cov.sum/cov.count)
        diag.G0 = diag(G.0)
    }

    if (makePD) {
        npc.0 <- {
            tmp <- Matrix::nearPD(npc.0, corr = FALSE, keepDiag = FALSE, do2eigen = TRUE, trace = TRUE)
            as.matrix(tmp$mat)
        }
    }
    ### numerical integration for calculation of eigenvalues (see Ramsay & Silverman, Chapter 8)
    w <- quadWeights(argvals, method=integration)
    Wsqrt <- diag(sqrt(w))
    Winvsqrt <- diag(1/(sqrt(w)))
    V <- Wsqrt %*% npc.0 %*% Wsqrt
    evalues = eigen(V, symmetric = TRUE, only.values = TRUE)$values
    ###
    evalues = replace(evalues, which(evalues <= 0), 0)
    npc = ifelse(is.null(npc), min(which(cumsum(evalues)/sum(evalues) > pve)), npc)
    efunctions = matrix(Winvsqrt%*%eigen(V, symmetric = TRUE)$vectors[, seq(len = npc)], nrow = D, ncol = npc)
    evalues = eigen(V, symmetric = TRUE, only.values = TRUE)$values[1:npc]  # use correct matrix for eigenvalue problem
    cov.hat = efunctions %*% diag(evalues) %*% t(efunctions)
    ### numerical integration for estimation of sigma2
    T.len <- argvals[D] - argvals[1] # total interval length
    T1.min <- min(which(argvals >= argvals[1] + 0.25*T.len)) # left bound of narrower interval T1
    T1.max <- max(which(argvals <= argvals[D] - 0.25*T.len)) # right bound of narrower interval T1
    DIAG = (diag.G0 - diag(cov.hat))[T1.min :T1.max] # function values
    w2 <- quadWeights(argvals[T1.min:T1.max], method = integration)
    sigma2 <- max(weighted.mean(DIAG, w=w2, na.rm = TRUE), 0)

    ####
    D.inv = diag(1/evalues)
    Z = efunctions
    Y.tilde = Y.pred - matrix(mu, I.pred, D, byrow = TRUE)
    Yhat = matrix(0, nrow = I.pred, ncol = D)
    scores = matrix(NA, nrow = I.pred, ncol = npc)
    VarMats = vector("list", I.pred)
    for (i in 1:I.pred) VarMats[[i]] = matrix(NA, nrow = D, ncol = D)
    diag.var = matrix(NA, nrow = I.pred, ncol = D)
    crit.val = rep(0, I.pred)
    for (i.subj in 1:I.pred) {
        obs.points = which(!is.na(Y.pred[i.subj, ]))
        if (sigma2 == 0 & length(obs.points) < npc)
            stop("Measurement error estimated to be zero and there are fewer observed points than PCs; scores cannot be estimated.")
        Zcur = matrix(Z[obs.points, ], nrow = length(obs.points), ncol = dim(Z)[2])
        scores[i.subj, ] = solve(crossprod(Zcur) + sigma2 * D.inv,
          t(Zcur) %*% (Y.tilde[i.subj, obs.points]))
        Yhat[i.subj, ] = t(as.matrix(mu)) + scores[i.subj, ] %*% t(efunctions)
        if (var) {
            VarMats[[i.subj]] = sigma2 * Z %*% ZtZ_sD.inv %*% t(Z)
            diag.var[i.subj, ] = diag(VarMats[[i.subj]])
            if (simul & sigma2 != 0) {
                norm.samp = mvrnorm(2500, mu = rep(0, D), Sigma = VarMats[[i.subj]])/matrix(sqrt(diag(VarMats[[i.subj]])), nrow = 2500, ncol = D, byrow = TRUE)
                crit.val[i.subj] = quantile(apply(abs(norm.samp), 1, max), sim.alpha)
            }
        }
    }

    ret.objects = c("Yhat", "scores", "mu", "efunctions", "evalues", "npc")
    if (var) {
        ret.objects = c(ret.objects, "sigma2", "diag.var", "VarMats")
        if (simul) ret.objects = c(ret.objects, "crit.val")
    }
    ret = lapply(1:length(ret.objects), function(u) get(ret.objects[u]))
    names(ret) = ret.objects
    return(ret)
}
