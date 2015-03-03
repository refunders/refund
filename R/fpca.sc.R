# npc=1 seems to give error
fpca.sc <- function(Y=NULL, ydata = NULL, Y.pred=NULL, argvals = NULL, random.int = FALSE, 
         nbasis = 10, pve = .99, npc = NULL, var = FALSE, simul = FALSE, sim.alpha = .95,
         useSymm = FALSE, makePD = FALSE, center=TRUE, cov.est.method = 2) {
  
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
  
    if (is.null(argvals)) argvals = 1:D
  
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
            npc.0 = matrix(predict(gam(as.vector(G.0) ~ te(row.vec, col.vec, k = nbasis), weights = as.vector(cov.count)),     
                                   newdata = data.frame(row.vec = row.vec,col.vec = col.vec)), 
                                   D, D)
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
    evalues = eigen(npc.0, symmetric = TRUE, only.values = TRUE)$values
    evalues = replace(evalues, which(evalues <= 0), 0)
    npc = ifelse(is.null(npc), min(which(cumsum(evalues)/sum(evalues) > pve)), npc)
    efunctions = matrix(eigen(npc.0, symmetric = TRUE)$vectors[, seq(len = npc)], nrow = D, ncol = npc)
    evalues = eigen(npc.0, symmetric = TRUE, only.values = TRUE)$values[1:npc]
    cov.hat = efunctions %*% tcrossprod(diag(evalues, nrow = npc, ncol = npc), efunctions)
    DIAG = (diag.G0 - diag(cov.hat))[floor(D * 0.2):ceiling(D * 0.8)]
    sigma2 = max(mean(DIAG, na.rm = TRUE), 0)
    D.inv = diag(1/evalues, nrow = npc, ncol = npc)
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
        ZtZ_sD.inv = solve(crossprod(Zcur) + sigma2 * D.inv)
        scores[i.subj, ] = ZtZ_sD.inv %*% t(Zcur) %*% (Y.tilde[i.subj, obs.points])
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

