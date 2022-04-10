#' Multilevel functional principal components analysis with fast covariance estimation
#' 
#' Decompose dense or sparse multilevel functional observations using multilevel 
#' functional principal component analysis with the fast covariance estimation 
#' approach. 
#' 
#' The fast MFPCA approach (Cui et al., 2022+) uses FACE (Xiao et al., 2016) to estimate 
#' covariance functions and mixed model equations (MME) to predict 
#' scores for each level. As a result, it has lower computational complexity than 
#' MFPCA (Di et al., 2009) implemented in the \code{mfpca.sc} function, and
#' can be applied to decompose data sets with over 10000 subjects and over 10000 
#' dimensions.
#' 
#' 
#' @param Y A multilevel functional dataset on a regular grid stored in a matrix.
#' Each row of the data is the functional observations at one visit for one subject.
#' Missingness is allowed and need to be labeled as NA. The data must be specified.
#' @param id A vector containing the id information to identify the subjects. The
#' data must be specified.
#' @param visit A vector containing information used to identify the visits. 
#' If not provided, assume the visit id are 1,2,... for each subject.
#' @param twoway Logical, indicating whether to carry out twoway ANOVA and 
#' calculate visit-specific means. Defaults to \code{TRUE}.
#' @param weight The way of calculating covariance. \code{weight = "obs"} indicates
#' that the sample covariance is weighted by observations. \code{weight = "subj"} 
#' indicates that the sample covariance is weighted equally by subjects. Defaults to \code{"obs"}.
#' @param argvals A vector containing observed locations on the functional domain.
#' @param pve Proportion of variance explained. This value is used to choose the 
#' number of principal components for both levels.
#' @param npc Pre-specified value for the number of principal components. 
#' If given, this overrides \code{pve}.
#' @param p The degree of B-splines functions to use. Defaults to 3.
#' @param m The order of difference penalty to use. Defaults to 2.
#' @param knots Number of knots to use or the vectors of knots. Defaults to 35.
#' @param silent Logical, indicating whether to not display the name of each step.
#' Defaults to \code{TRUE}.
#' 
#' @return An object of class \code{mfpca} containing:
#' \item{Xhat}{FPC approximation (projection onto leading components)
#' of \code{Y}, estimated curves for all subjects and visits}
#' \item{Xhat.subject}{Estimated subject specific curves for all subjects}
#' \item{Y}{The observed data}
#' \item{scores}{A matrix of estimated FPC scores for level1 and level2.} 
#' \item{mu}{estimated mean function (or a vector of zeroes if \code{center==FALSE}).} 
#' \item{efunctions}{A matrix of estimated eigenfunctions of the functional
#' covariance, i.e., the FPC basis functions for levels 1 and 2.} 
#' \item{evalues}{Estimated eigenvalues of the covariance operator, i.e., variances 
#' of FPC scores for levels 1 and 2.}
#' \item{npc}{Number of FPCs: either the supplied \code{npc}, or the minimum
#' number of basis functions needed to explain proportion \code{pve} of the
#' variance in the observed curves for levels 1 and 2.} 
#' \item{sigma2}{Estimated measurement error variance.} 
#' \item{eta}{The estimated visit specific shifts from overall mean.}
#' 
#' @author Ruonan Li \email{rli20@@ncsu.edu}, Erjia Cui \email{ecui1@@jhmi.edu}
#' 
#' @references Cui, E., Li, R., Crainiceanu, C., and Xiao, L. (2022+). Fast multilevel
#' functional principal component analysis.
#' 
#' Di, C., Crainiceanu, C., Caffo, B., and Punjabi, N. (2009).
#' Multilevel functional principal component analysis. \emph{Annals of Applied
#' Statistics}, 3, 458--488.
#' 
#' Xiao, L., Ruppert, D., Zipunnikov, V., and Crainiceanu, C. (2016).
#' Fast covariance estimation for high-dimensional functional data.  
#' \emph{Statistics and Computing}, 26, 409-421.
#' 
#' @export
#' @importFrom splines spline.des
#' @importFrom mgcv s smooth.construct gam predict.gam
#' @importFrom MASS ginv
#' @importFrom simex diag.block
#' @importFrom Matrix crossprod
#' @importFrom stats smooth.spline optim
#' 
#' 
#' @examples 
#' data(DTI)
#' mfpca.DTI <- mfpca.face(Y = DTI$cca, id = DTI$ID, twoway = TRUE)

mfpca.face <- function(Y, id, visit = NULL, twoway = TRUE, weight = "obs", argvals = NULL,
                        pve = 0.99, npc = NULL, p = 3, m = 2, knots = 35, silent = TRUE){
  
  pspline.setting.mfpca <- function(x,knots=35,p=3,m=2,weight=NULL,type="full",
                                    knots.option="equally-spaced"){
    # x: the marginal data points
    # knots: the list of interior knots or the numbers of interior knots
    # p: degrees for B-splines, with defaults values 3
    # m: orders of difference penalty, with default values 2
    # knots.option: type of knots placement, with default values "equally-spaced"
    
    #require(splines)
    #require(Matrix)
    
    ### design matrix 
    K = length(knots)-2*p-1
    B = spline.des(knots=knots, x=x, ord=p+1, outer.ok=TRUE,sparse=TRUE)$design
    
    bs = "ps"
    
    if(knots.option == "quantile"){
      bs = "bs"
    }
    
    s.object = s(x=x, bs=bs, k=K+p, m=c(p-1,2), sp=NULL)
    object  = smooth.construct(s.object,data = data.frame(x=x),knots=list(x=knots))
    P =  object$S[[1]]
    if(knots.option == "quantile") P = P / max(abs(P))*10 # rescaling
    
    if(is.null(weight)) weight <- rep(1,length(x))
    
    if(type=="full"){
      
      Sig = crossprod(matrix.multiply.mfpca(B,weight,option=2),B)
      eSig = eigen(Sig)
      V = eSig$vectors
      E = eSig$values
      if(min(E)<=0.0000001) {#cat("Warning! t(B)%*%B is singular!\n");
        #cat("A small identity matrix is added!\n");
        E <- E + 0.000001;
      }
      Sigi_sqrt = matrix.multiply.mfpca(V,1/sqrt(E))%*%t(V)
      
      tUPU = Sigi_sqrt%*%(P%*%Sigi_sqrt)
      Esig = eigen(tUPU,symmetric=TRUE)
      U = Esig$vectors
      s = Esig$values
      s[(K+p-m+1):(K+p)]=0
      A = B%*%(Sigi_sqrt%*%U)
    }
    
    if(type=="simple"){
      A = NULL
      s = NULL
      Sigi_sqrt = NULL
      U = NULL
    }
    
    List = list("A" = A, "B" = B, "s" = s, "Sigi.sqrt" = Sigi_sqrt, "U" = U, "P" = P)
    
    return(List)
  }
  
  ## a function supposed to be in the refund package
  quadWeights.mfpca <- function(argvals, method = "trapezoidal")
  {
    ret <- switch(method,
                  trapezoidal = {D <- length(argvals)
                  1/2*c(argvals[2] - argvals[1], argvals[3:D] -argvals[1:(D-2)], argvals[D] - argvals[D-1])},
                  midpoint = c(0,diff(argvals)),
                  stop("function quadWeights: choose either trapezoidal or midpoint quadrature rule"))
    
    return(ret)  
  }
  
  ##################################################################################
  ## Organize the input
  ##################################################################################
  if(silent == FALSE) print("Organize the input")
  
  stopifnot((!is.null(Y) & !is.null(id)))
  stopifnot(is.matrix(Y))
  
  ## specify visit variable if not provided
  if (!is.null(visit)){ 
    visit <- as.factor(visit)
  }else{ ## if visit is not provided, assume the visit id are 1,2,... for each subject
    visit <- as.factor(ave(id, id, FUN=seq_along))
  }
  id <- as.factor(id) ## convert id into a factor
  
  ## organize data into one data frame
  df <- data.frame(id = id, visit = visit, Y = I(Y))
  rm(id, visit, Y)
  
  ## derive several variables that will be used later
  J <- length(levels(df$visit)) ## number of visits
  L <- ncol(df$Y) ## number of observations along the domain
  nVisits <- data.frame(table(df$id))  ## calculate number of visits for each subject
  colnames(nVisits) = c("id", "numVisits")
  ID = sort(unique(df$id)) ## id of each subject
  I <- length(ID) ## number of subjects
  ## assume observations are equally-spaced on [0,1] if not specified
  if (is.null(argvals))  argvals <- seq(0, 1, length.out=L) 
  
  
  ##################################################################################
  ## Estimate population mean function (mu) and visit-specific mean function (eta)
  ##################################################################################
  if(silent == FALSE) print("Estimate population and visit-specific mean functions")
  
  meanY <- colMeans(df$Y, na.rm = TRUE)
  fit_mu <- gam(meanY ~ s(argvals))
  mu <- as.vector(predict(fit_mu, newdata = data.frame(argvals = argvals)))
  rm(meanY, fit_mu)
  
  mueta = matrix(0, L, J) 
  eta = matrix(0, L, J) ## matrix to store visit-specific means
  colnames(mueta) <- colnames(eta) <- levels(df$visit)
  Ytilde <- matrix(NA, nrow = nrow(df$Y), ncol = ncol(df$Y))
  if(twoway==TRUE) {
    for(j in 1:J) {
      ind_j <- which(df$visit == levels(df$visit)[j])
      if(length(ind_j) > 1){
        meanYj <- colMeans(df$Y[ind_j,], na.rm=TRUE)
      }else{
        meanYj <- df$Y[ind_j,]
      }
      fit_mueta <- gam(meanYj ~ s(argvals))
      mueta[,j] <- predict(fit_mueta, newdata = data.frame(argvals = argvals))
      eta[,j] <- mueta[,j] - mu
      Ytilde[ind_j,] <- df$Y[ind_j,] - matrix(mueta[,j], nrow = length(ind_j), ncol = L, byrow = TRUE)
    }
    rm(meanYj, fit_mueta, ind_j, j)
  } else{
    Ytilde <- df$Y - matrix(mu, nrow = nrow(df$Y), ncol = L, byrow = TRUE)
  }
  df$Ytilde <- I(Ytilde) ## Ytilde is the centered multilevel functional data
  rm(Ytilde)
  
  
  ##################################################################################
  ## FACE preparation: see Xiao et al. (2016) for details
  ##################################################################################
  if(silent == FALSE) print("Prepare ingredients for FACE")
  
  ## Specify the knots of B-spline basis
  if(length(knots)==1){
    if(knots+p>=L) cat("Too many knots!\n")
    stopifnot(knots+p<L)
    
    K.p <- knots
    knots <- seq(-p, K.p+p, length=K.p+1+2*p)/K.p
    knots <- knots*(max(argvals)-min(argvals)) + min(argvals)
  }
  if(length(knots)>1) K.p <- length(knots)-2*p-1
  if(K.p>=L) cat("Too many knots!\n")
  stopifnot(K.p < L)
  c.p <- K.p + p #the number of B-spline basis functions
  
  ## Precalculation for smoothing
  List <- pspline.setting.mfpca(argvals, knots, p, m)
  B <- List$B #B is the J by c design matrix
  Sigi.sqrt <- List$Sigi.sqrt #(t(B)B)^(-1/2)
  s <- List$s #eigenvalues of Sigi_sqrt%*%(P%*%Sigi_sqrt)
  U <- List$U #eigenvectors of Sigi_sqrt%*%(P%*%Sigi_sqrt)
  A0 <- Sigi.sqrt %*% U
  G <- crossprod(B) / nrow(B)
  eig_G <- eigen(G, symmetric = T)
  G_half <- eig_G$vectors %*% diag(sqrt(eig_G$values)) %*% t(eig_G$vectors)
  G_invhalf <- eig_G$vectors %*% diag(1/sqrt(eig_G$values)) %*% t(eig_G$vectors)
  Bnew <- as.matrix(B %*% G_invhalf) #the transformed basis functions
  Anew <- G_half %*% A0
  rm(List, Sigi.sqrt, U, G, eig_G, G_half)
  
  ##################################################################################
  ## Impute missing data of Y using FACE and estimate the total covariance (Kt)
  ##################################################################################
  if(silent == FALSE) print("Estimate the total covariance (Kt)")
  
  Ji <- as.numeric(table(df$id))
  diagD <- rep(Ji, Ji)
  smooth.Gt = face.Cov.mfpca(Y=unclass(df$Ytilde), argvals, A0, B, Anew, Bnew, G_invhalf, s)
  ## impute missing data of Y using FACE approach
  if(sum(is.na(df$Ytilde))>0){
    df$Ytilde[which(is.na(df$Ytilde))] <- smooth.Gt$Yhat[which(is.na(df$Ytilde))]
  }
  if(weight=="subj"){
    YH <- unclass(df$Ytilde)*sqrt(nrow(df$Ytilde)/(I*diagD))
    smooth.Gt <- face.Cov.mfpca(Y=YH, argvals, A0, B, Anew, Bnew, G_invhalf, s)
    rm(YH)
  }
  diag_Gt <- colMeans(df$Ytilde^2)
  
  
  ##################################################################################
  ## Estimate principal components of the within covariance (Kw)
  ##################################################################################
  if(silent == FALSE) print("Estimate principal components of the within covariance (Kw)")
  
  inx_row_ls <- split(1:nrow(df$Ytilde), f=factor(df$id, levels=unique(df$id)))
  Ysubm <- t(vapply(inx_row_ls, function(x) colSums(df$Ytilde[x,,drop=FALSE],na.rm=TRUE), numeric(L)))
  if(weight=="obs"){
    weights <- sqrt(nrow(df$Ytilde)/(sum(diagD) - nrow(df$Ytilde)))
    YR <-  do.call("rbind",lapply(1:I, function(x) {
      weights*sqrt(Ji[x]) * t(t(df$Ytilde[inx_row_ls[[x]],,drop=FALSE]) - Ysubm[x,]/Ji[x])
    }))
  }
  if(weight=="subj"){
    weights <- sqrt(nrow(df$Ytilde)/sum(Ji>1))
    YR <-  do.call("rbind",lapply(1:I, function(x) {
      if(Ji[x]>1) return((weights/sqrt(Ji[x]-1)) * t(t(df$Ytilde[inx_row_ls[[x]],,drop=FALSE]) - Ysubm[x,]/Ji[x]))
    }))
  }
  smooth.Gw <- face.Cov.mfpca(Y=YR, argvals, A0, B, Anew, Bnew, G_invhalf, s)
  sigma.Gw <- smooth.Gw$evalues # raw eigenvalues of Gw
  per <- cumsum(sigma.Gw)/sum(sigma.Gw)
  N.Gw <- ifelse (is.null(npc), min(which(per>pve)), min(npc, length(sigma.Gw)))
  smooth.Gw$efunctions <- smooth.Gw$efunctions[,1:N.Gw]
  smooth.Gw$evalues <- smooth.Gw$evalues[1:N.Gw]
  rm(Ji, diagD, inx_row_ls, weights, weight, Ysubm, YR, B, Anew, G_invhalf, s, per, N.Gw, sigma.Gw)
  
  
  ##################################################################################
  ## Estimate principal components of the between covariance (Kb)
  ##################################################################################
  if(silent == FALSE) print("Estimate principal components of the between covariance (Kb)")
  
  temp = smooth.Gt$decom - smooth.Gw$decom
  Eigen <- eigen(temp,symmetric=TRUE)
  Sigma <- Eigen$values
  d <- Sigma[1:c.p]
  d <- d[d>0]
  per <- cumsum(d)/sum(d)
  N.Gb <- ifelse (is.null(npc), min(which(per>pve)), min(npc, length(d)))
  smooth.Gb <- list(evalues=Sigma[1:N.Gb], efunctions=Bnew%*%Eigen$vectors[,1:N.Gb])
  rm(smooth.Gt, temp, Eigen, Sigma, d, per, N.Gb, Bnew)
  
  
  ###########################################################################################
  ## Estimate eigenvalues and eigenfunctions at two levels by calling the 
  ## eigenfunction (in R "base" package) on discretized covariance matrices.
  ###########################################################################################
  if(silent == FALSE) print("Estimate eigenvalues and eigenfunctions at two levels")
  
  efunctions <- list(level1=as.matrix(smooth.Gb$efunctions), level2=as.matrix(smooth.Gw$efunctions))
  evalues <- list(level1=smooth.Gb$evalues,level2=smooth.Gw$evalues)
  npc <- list(level1=length(evalues[[1]]), level2=length(evalues[[2]]))
  names(efunctions) <- names(evalues) <- names(npc) <- c("level1", "level2")
  rm(smooth.Gb, smooth.Gw)
  
  
  ###################################################################
  # Estimate the measurement error variance (sigma^2)
  ###################################################################
  if(silent == FALSE) print("Estimate the measurement error variance (sigma^2)")
  
  cov.hat <- lapply(c("level1", "level2"), function(x) colSums(t(efunctions[[x]]^2)*evalues[[x]]))
  T.len <- argvals[L] - argvals[1]
  T1.min <- min(which(argvals >= argvals[1] + 0.25 * T.len))
  T1.max <- max(which(argvals <= argvals[L] - 0.25 * T.len))
  DIAG <- (diag_Gt - cov.hat[[1]] - cov.hat[[2]])[T1.min:T1.max]
  w2 <- quadWeights.mfpca(argvals[T1.min:T1.max], method = "trapezoidal")
  sigma2 <- max(weighted.mean(DIAG, w = w2, na.rm = TRUE), 0) ## estimated measurement error variance
  rm(cov.hat, T.len, T1.min, T1.max, DIAG, w2)
  
  
  ###################################################################
  # Estimate the principal component scores
  ###################################################################
  if(silent == FALSE) print("Estimate principal component scores")
  
  ## estimated subject-visit and subject-specific underlying smooth curves
  Xhat <- Xhat.subject <- matrix(0, nrow(df$Y), L) 
  phi1 <- efunctions[[1]] ## extract eigenfunctions for simpler notations
  phi2 <- efunctions[[2]]
  score1 <- matrix(0, I, npc[[1]]) ## matrices storing scores of two levels
  score2 <- matrix(0, nrow(df$Y), npc[[2]])
  
  unVisits <- unique(nVisits$numVisits) ## unique number of visits
  if(length(unVisits) < I){
    for(j in 1:length(unVisits)){
      Jm <- unVisits[j]
      ## calculate block matrices
      if(sigma2 < 1e-4){
        A <- Jm * (t(phi1) %*% phi1)
        B <- matrix(rep(t(phi1) %*% phi2, Jm), nrow = npc[[1]])
        temp <- ginv(t(phi2) %*% phi2)
      }else{
        if(length(evalues[[1]])==1){
          A <- Jm * (t(phi1) %*% phi1) / sigma2 + 1 / evalues[[1]]
        }else{
          A <- Jm * (t(phi1) %*% phi1) / sigma2 + diag(1 / evalues[[1]])
        }
        B = matrix(rep(t(phi1) %*% phi2 / sigma2, Jm), nrow = npc[[1]])
        if(length(evalues[[2]])==1){
          temp = ginv(t(phi2) %*% phi2 / sigma2 + 1 / evalues[[2]])
        }else{
          temp = ginv(t(phi2) %*% phi2 / sigma2 + diag(1 / evalues[[2]]))
        }
      }
      C <- t(B)
      invD = diag.block(temp, Jm)
      
      ## calculate inverse of each block components separately
      MatE <- ginv(A-B%*%invD%*%C)
      MatF <- -invD%*%C%*%MatE
      MatG <- -MatE%*%B%*%invD
      MatH <- invD - invD%*%C%*%MatG
      Mat1 <- cbind(MatE,MatG)
      Mat2 <- cbind(MatF,MatH)
      
      ## estimate the principal component scores using MME
      ind.Jm <- nVisits$id[which(nVisits$numVisits == Jm)]
      YJm <- matrix(df$Ytilde[which(df$id %in% ind.Jm),], ncol = L)
      int1 <- rowsum(df$Ytilde[which(df$id %in% ind.Jm),] %*% phi1, rep(1:length(ind.Jm), each = Jm))
      int2 <- t(matrix(t(df$Ytilde[which(df$id %in% ind.Jm),] %*% phi2), nrow = npc[[2]]*Jm))
      int <- cbind(int1, int2)
      if(sigma2 >= 1e-4){
        int <- int / sigma2
      }
      score1[which(nVisits$id %in% ind.Jm),] <- int %*% t(Mat1)
      score2[which(df$id %in% ind.Jm),] <- t(matrix(Mat2 %*% t(int), nrow = npc[[2]]))
      
      temp <- score1[which(nVisits$id %in% ind.Jm),] %*% t(phi1)
      Xhat.subject[which(df$id %in% ind.Jm),] <- temp[rep(1:length(ind.Jm), each = Jm),]
      Xhat[which(df$id %in% ind.Jm),] <- Xhat.subject[which(df$id %in% ind.Jm),] + 
        score2[which(df$id %in% ind.Jm),] %*% t(phi2)
    }
    for(g in 1:length(levels(df$visit))){
      ind.visit <- which(df$visit == levels(df$visit)[g])
      Xhat.subject[ind.visit,] <- t(t(Xhat.subject[ind.visit,]) + mu + eta[,levels(df$visit)[g]])
      Xhat[ind.visit,] <- t(t(Xhat[ind.visit,]) + mu + eta[,levels(df$visit)[g]])
    }
    rm(YJm, g, ind.visit, ind.Jm)
    
  }else{
    for(m in 1:I){
      Jm <- nVisits[m, 2]  ## number of visits for mth subject
      ## calculate block matrices
      if(sigma2 < 1e-4){
        A <- Jm * (t(phi1) %*% phi1)
        B <- matrix(rep(t(phi1) %*% phi2, Jm), nrow = npc[[1]])
        temp <- ginv(t(phi2) %*% phi2)
      }else{
        if(length(evalues[[1]])==1){
          A <- Jm * (t(phi1) %*% phi1) / sigma2 + 1 / evalues[[1]]
        }else{
          A <- Jm * (t(phi1) %*% phi1) / sigma2 + diag(1 / evalues[[1]])
        }
        B = matrix(rep(t(phi1) %*% phi2 / sigma2, Jm), nrow = npc[[1]])
        if(length(evalues[[2]])==1){
          temp = ginv(t(phi2) %*% phi2 / sigma2 + 1 / evalues[[2]])
        }else{
          temp = ginv(t(phi2) %*% phi2 / sigma2 + diag(1 / evalues[[2]]))
        }
      }
      C <- t(B)
      invD = diag.block(temp, Jm)
      
      ## calculate inverse of each block components separately
      MatE <- ginv(A-B%*%invD%*%C)
      MatF <- -invD%*%C%*%MatE
      MatG <- -MatE%*%B%*%invD
      MatH <- invD - invD%*%C%*%MatG
      Mat1 <- cbind(MatE,MatG)
      Mat2 <- cbind(MatF,MatH)
      
      ## estimate the principal component scores
      int1 <- colSums(matrix(df$Ytilde[df$id==ID[m],], ncol = L) %*% phi1)
      int2 <- matrix(df$Ytilde[df$id==ID[m],], ncol = L) %*% phi2
      if(sigma2 < 1e-4){
        int <- c(int1, as.vector(t(int2)))
      }else{
        int <- c(int1, as.vector(t(int2))) / sigma2
      }
      score1[m,] <- Mat1 %*% int
      score2[which(df$id==ID[m]),] <- matrix(Mat2 %*% int, ncol = npc[[2]], byrow=TRUE)
      for (j in which(df$id==ID[m])) {
        Xhat.subject[j,] <- as.matrix(mu) + eta[,df$visit[j]] + as.vector(phi1 %*% score1[m,])
        Xhat[j,] <- Xhat.subject[j,] + as.vector(phi2 %*% score2[j,])
      }
    }
  }
  scores <- list(level1 = score1, level2 = score2)
  
  rm(A, B, C, int, int1, int2, invD, Mat1, Mat2, MatE, MatF, MatG, MatH, temp, j, Jm, unVisits, phi1, phi2, score1, score2)
  
  
  ###################################################################
  # Organize the results
  ###################################################################
  if(silent == FALSE) print("Organize the results")
  
  res <- list(Xhat = Xhat, Xhat.subject = Xhat.subject, mu = mu, eta = eta, scores = scores, 
              efunctions = efunctions, evalues = evalues, npc = npc, sigma2 = sigma2, Y = df$Y)
  
  rm(df, efunctions, eta, evalues, mueta, nVisits, npc, scores, Xhat, Xhat.subject, 
     argvals, diag_Gt, I, ID, J, mu, pve, L, sigma2)
  
  return(res)
}
