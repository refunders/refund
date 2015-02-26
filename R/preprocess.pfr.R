preprocess.pfr <- function (subj=NULL, covariates = NULL, funcs, kz = NULL, kb = NULL, nbasis=10, funcs.new=NULL, smooth.option="fpca.sc",pve=0.99){
## TBD:  make kz,kb take the value from pfr.obj in predict.pfr() calls.
 
  
  N_subj = length(unique(subj))
  p = ifelse(is.null(covariates), 0, dim(covariates)[2])
  ## handle if funcs is a list (multiple predictors) or matrix (single predictor)
  if (is.matrix(funcs)) {
    Funcs = list(length = 1)
    Funcs[[1]] = funcs
  }else {
    Funcs = funcs
  }
  ## handle if funcs.new is a list (multiple predictors) or matrix (single predictor)
  if (is.matrix(funcs.new)) {
    Funcs.new = list(length = 1)
    Funcs.new[[1]] = funcs.new
  }else {
    Funcs.new = funcs.new
  }
  N.Pred = length(Funcs)
  ### change: July 11, 2013
  if(!is.null(kz)){
  if(length(kz)==1) kz = rep(kz,N.Pred)
  if(length(kz)!=N.Pred) stop("Length of kz is not the number of predictors\n")
  }
  kz.adj = rep(NA, N.Pred)
  ### end of change
  
  ## outcome length "o.len" is number of rows of original data if no predictive data provided
  ## if predictive data provided, then must be the number of rows for predictive data
  if(is.null(funcs.new)){o.len <- nrow(as.matrix(Funcs[[1]]))
                       }else{o.len <- nrow(as.matrix(Funcs.new[[1]]))}
  t <- phi <- FPCA <- psi <- C <- J <- CJ <- D <- list()
  
  ## obtain FPCA decomposition of each predictor
  if (smooth.option=="fpca.sc"){
      # using fpca.sc()
      for(i in 1:N.Pred){
        t[[i]] = seq(0, 1, length = dim(Funcs[[i]])[2])
        ## for fpca() (not fpca.sc()):  when I have K=kz, I have to lower kz a bit
        FPCA[[i]] = fpca.sc(Y = Funcs[[i]], Y.pred = Funcs.new[[i]], pve=pve, nbasis=nbasis, npc=kz[i]) 
        psi[[i]] = FPCA[[i]]$efunctions
        C[[i]]=FPCA[[i]]$scores
        kz.adj[i]=FPCA[[i]]$npc
        ## what psi and C are if using fpca()
        ##    psi[[i]] = FPCA[[i]]$phi.hat
        ##    C[[i]]=FPCA[[i]]$xi.hat
      }
  }
  
  if (smooth.option=="fpca.face"){
    # using face
    for(i in 1:N.Pred){
  
        Funcs[[i]] = apply(Funcs[[i]],2,function(x){x-0*mean(x,na.rm=TRUE)})
        if(!is.null(Funcs.new[[i]])) Funcs.new[[i]] = apply(Funcs.new[[i]],2,function(x){x-0*mean(x,na.rm=TRUE)})
        
        t[[i]] = seq(0, 1, length = dim(Funcs[[i]])[2])
        #if (length(t[[i]])>70) nbasis = max(nbasis,35)  
        FPCA[[i]] = fpca.face(Y = Funcs[[i]], Y.pred = Funcs.new[[i]], knots=nbasis,pve = pve)
        if (is.null(kz[i]) || kz[i]>dim(FPCA[[i]]$eigenvectors)[2]){
            psi[[i]] = FPCA[[i]]$eigenvectors
            C[[i]]=FPCA[[i]]$scores*sqrt(dim(Funcs[[i]])[2])
            kz.adj[i] =  dim(FPCA[[i]]$eigenvectors)[2]
            cat("For the ", i, "-th functional predictor, the number of PCs changes to", kz.adj[i],"\n");
            cat("For details, see the manual\n");
        }
        else {
            #cat(kz[i],"\n",dim(FPCA[[i]]$eigenvectors),"\n")
            psi[[i]] = FPCA[[i]]$eigenvectors[,1:kz[i]]
            C[[i]]=FPCA[[i]]$scores[,1:kz[i]]*sqrt(dim(Funcs[[i]])[2])
            kz.adj[i] =  kz[i]
        }
        
    }
        
  }
  
  
  ## construct phi for b-splines; J and CJ.
  for(i in 1:N.Pred){
    phi[[i]] = cbind(1, bs(t[[i]], df=kb-1, intercept=FALSE, degree=3))
    J[[i]] = t(psi[[i]]) %*% phi[[i]]
    CJ[[i]] = C[[i]] %*% J[[i]]
  }
  ## setup gam() penalty matrices
  if(!is.null(subj)){
    ## setup matrix for random intercepts
    Z1 = matrix(0, nrow = o.len, ncol = N_subj)
    for (i in 1:length(unique(subj))) {
      Z1[which(subj == unique(subj)[i]), i] = 1
    }
    colnames(Z1)=c(paste("i",1:dim(Z1)[2], sep=""))
    ## D[[1]] is the penalty on random intercepts
    D[[1]] = diag(c(rep(0, 1 + p),
       rep(1, N_subj),
       rep(0, length = N.Pred * (kb))))
    totD <- N.Pred+1
    startD <- 2
  }else{
    ## in the case that is.null(subj) set Z1 NULL
    Z1 <- NULL
    totD <- N.Pred
    startD <- 1
  }
  ## set up penalty sub-matrix for b-splines
  ## to be embedded in D[[startD]],...,D[[totD]]
  temp=matrix(0, nrow=kb-1, ncol=kb-1)
  for(ii in 1:(kb-1)){
    for(jj in 1:(kb-1)){
      temp[ii,jj]=min(ii,jj)-1
    }
  }
  spl.pen = matrix(1, nrow=kb-1, ncol=kb-1)+temp
  Dinv=solve(spl.pen)
  for (i in startD:totD) {
    D[[i]] = magic::adiag( diag(c(rep(0, 1 + p + N_subj*!is.null(subj)), rep(0, kb * (i - startD)) , rep(0, 1))),
       Dinv,
       diag(rep(0, kb * (totD - i))))
  }
  ## set up gam() design matrix X
  X = cbind(rep(1,o.len), covariates, Z1)
  for (i in 1:N.Pred) {
    X = cbind(X, CJ[[i]])
  }
  ## set up fixed and random effects matrix for (R)LRTSim(), lme()
  fixed.mat = X[,1:(1+p)]
  rand.mat = Z1
  for (i in 1:N.Pred) {
    fixed.mat = cbind(fixed.mat, CJ[[i]][,1])
    rand.mat  = cbind(rand.mat , CJ[[i]][,2:kb])
  }
  ## only need (Z1, C, psi) for predict.pfr
  ## having fixed.mat, rand.mat useful for rlrt.pfr
  ## pfr() requires 
  ret <- list( X, D,  phi, psi, C, J, CJ,
              Z1, subj,
              fixed.mat, rand.mat, N_subj, p, N.Pred,
              kz, kz.adj, kb, nbasis,
              totD,
              funcs, covariates, smooth.option)
  names(ret) <- c("X", "D", "phi", "psi", "C", "J", "CJ",
                  "Z1", "subj",
                  "fixed.mat", "rand.mat", "N_subj", "p", "N.Pred",
                  "kz", "kz.adj", "kb", "nbasis",
                  "totD",
                  "funcs", "covariates", "smooth.option")
  ret
}
