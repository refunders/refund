#' Penalized flexible functional regression
#'
#' Implements additive regression for functional and scalar covariates and
#' functional responses. This function is a wrapper for \code{mgcv}'s
#' \code{\link[mgcv]{gam}} and its siblings to fit models of the general form
#' \cr \eqn{E(Y_i(t)) = g(\mu(t) + \int X_i(s)\beta(s,t)ds + f(z_{1i}, t) +
#' f(z_{2i}) + z_{3i} \beta_3(t) + \dots )}\cr with a functional (but not
#' necessarily continuous) response \eqn{Y(t)}, response function \eqn{g},
#' (optional) smooth intercept \eqn{\mu(t)}, (multiple) functional covariates
#' \eqn{X(t)} and scalar covariates \eqn{z_1}, \eqn{z_2}, etc.
#'
#' @section Details: The routine can estimate \enumerate{ \item linear
#'   functional effects of scalar (numeric or factor) covariates that vary
#'   smoothly over \eqn{t} (e.g. \eqn{z_{1i} \beta_1(t)}, specified as
#'   \code{~z1}), \item nonlinear, and possibly multivariate functional effects
#'   of (one or multiple) scalar covariates \eqn{z} that vary smoothly over the
#'   index \eqn{t} of \eqn{Y(t)} (e.g. \eqn{f(z_{2i}, t)}, specified in the
#'   \code{formula} simply as \code{~s(z2)}) \item (nonlinear) effects of scalar
#'   covariates that are constant over \eqn{t} (e.g. \eqn{f(z_{3i})}, specified
#'   as \code{~c(s(z3))}, or \eqn{\beta_3 z_{3i}}, specified as \code{~c(z3)}),
#'   \item function-on-function regression terms (e.g. \eqn{\int
#'   X_i(s)\beta(s,t)ds}, specified as \code{~ff(X, yindex=t, xindex=s)}, see
#'   \code{\link{ff}}). Terms given by \code{\link{sff}} and \code{\link{ffpc}}
#'   provide nonlinear and FPC-based effects of functional covariates,
#'   respectively. \item concurrent effects of functional covariates \code{X}
#'   measured on the same grid as the response  are specified as follows:
#'   \code{~s(x)} for a smooth, index-varying effect \eqn{f(X(t),t)}, \code{~x}
#'   for a linear index-varying effect \eqn{X(t)\beta(t)}, \code{~c(s(x))} for a
#'   constant nonlinear effect \eqn{f(X(t))}, \code{~c(x)} for a constant linear
#'   effect \eqn{X(t)\beta}. \item Smooth functional random intercepts
#'   \eqn{b_{0g(i)}(t)} for a grouping variable \code{g} with levels \eqn{g(i)}
#'   can be specified via \code{~s(g, bs="re")}), functional random slopes
#'   \eqn{u_i b_{1g(i)}(t)} in a numeric variable \code{u} via \code{~s(g, u,
#'   bs="re")}). Scheipl, Staicu, Greven (2013) contains code examples for
#'   modeling correlated functional random intercepts using
#'   \code{\link[mgcv]{mrf}}-terms. } Use the \code{c()}-notation to denote
#'   model terms that are constant over the index of the functional response.\cr
#'
#'   Internally, univariate smooth terms without a \code{c()}-wrapper are
#'   expanded into bivariate smooth terms in the original covariate and the
#'   index of the functional response. Bivariate smooth terms (\code{s(), te()}
#'   or \code{t2()}) without a \code{c()}-wrapper are expanded into trivariate
#'   smooth terms in the original covariates and the index of the functional
#'   response. Linear terms for scalar covariates or categorical covariates are
#'   expanded into varying coefficient terms, varying smoothly over the index of
#'   the functional response. For factor variables, a separate smooth function
#'   with its own smoothing parameter is estimated for each level of the
#'   factor.\cr \cr The marginal spline basis used for the index of the the
#'   functional response is specified via the \emph{global} argument
#'   \code{bs.yindex}. If necessary, this can be overriden for any specific term
#'   by supplying a \code{bs.yindex}-argument to that term in the formula, e.g.
#'   \code{~s(x, bs.yindex=list(bs="tp", k=7))} would yield a tensor product
#'   spline over \code{x} and the index of the response in which the marginal
#'   basis for the index of the response are 7 cubic thin-plate spline functions
#'   (overriding the global default for the basis and penalty on the index of
#'   the response given by the \emph{global} \code{bs.yindex}-argument).\cr Use
#'   \code{~-1 + c(1) + ...} to specify a model with only a constant and no
#'   functional intercept. \cr
#'
#'   The functional covariates have to be supplied as a \eqn{n} by <no. of
#'   evaluations> matrices, i.e. each row is one functional observation. For
#'   data on a regular grid, the functional response is supplied in the same
#'   format, i.e. as a matrix-valued entry in \code{data},  which can contain
#'   missing values.\cr
#'
#'   If the functional responses are \emph{sparse or irregular} (i.e., not
#'   evaluated on the same evaluation points across all observations), the
#'   \code{ydata}-argument can be used to specify the responses: \code{ydata}
#'   must be a \code{data.frame} with 3 columns called \code{'.obs', '.index',
#'   '.value'} which specify which curve the point belongs to
#'   (\code{'.obs'}=\eqn{i}), at which \eqn{t} it was observed
#'   (\code{'.index'}=\eqn{t}), and the observed value
#'   (\code{'.value'}=\eqn{Y_i(t)}). Note that the vector of unique sorted
#'   entries in \code{ydata$.obs} must be equal to \code{rownames(data)} to
#'   ensure the correct association of entries in \code{ydata} to the
#'   corresponding rows of \code{data}. For both regular and irregular
#'   functional responses, the model is then fitted with the data in long
#'   format, i.e., for data on a grid the rows of the matrix of the functional
#'   response evaluations \eqn{Y_i(t)} are stacked into one long vector and the
#'   covariates are expanded/repeated correspondingly. This means the models get
#'   quite big fairly fast, since the effective number of rows in the design
#'   matrix is number of observations times number of evaluations of \eqn{Y(t)}
#'   per observation.\cr
#'
#'   Note that \code{pffr} does not use \code{mgcv}'s default identifiability
#'   constraints (i.e., \eqn{\sum_{i,t} \hat f(z_i, x_i, t) = 0} or
#'   \eqn{\sum_{i,t} \hat f(x_i, t) = 0}) for tensor product terms whose
#'   marginals include the index \eqn{t} of the functional response.  Instead,
#'   \eqn{\sum_i \hat f(z_i, x_i, t) = 0} for all \eqn{t} is enforced, so that
#'   effects varying over \eqn{t} can be interpreted as local deviations from
#'   the global functional intercept. This is achieved by using
#'   \code{\link[mgcv]{ti}}-terms with a suitably modified \code{mc}-argument.
#'   Note that this is not possible if \code{algorithm='gamm4'} since only
#'   \code{t2}-type terms can then be used and these modified constraints are
#'   not available for \code{t2}. We recommend using centered scalar covariates
#'   for terms like \eqn{z \beta(t)} (\code{~z}) and centered functional
#'   covariates with \eqn{\sum_i X_i(t) = 0} for all \eqn{t} in \code{ff}-terms
#'   so that the global functional intercept can be interpreted as the global
#'   mean function.
#'
#'   The \code{family}-argument can be used to specify all of the response
#'   distributions and link functions described in
#'   \code{\link[mgcv]{family.mgcv}}. Note that  \code{family = "gaulss"} is
#'   treated in a special way: Users can supply the formula for the variance by
#'   supplying a special argument \code{varformula}, but this is not modified in
#'   the way that the \code{formula}-argument is but handed over to the fitter
#'   directly, so this is for expert use only. If \code{varformula} is not
#'   given, \code{pffr} will use the parameters from argument \code{bs.int} to
#'   define a spline basis along the index of the response, i.e., a smooth
#'   variance function over $t$ for responses $Y(t)$.
#'
#' @param formula a formula with special terms as for \code{\link[mgcv]{gam}},
#'   with additional special terms \code{\link{ff}(), \link{sff}(),
#'   \link{ffpc}(), \link{pcre}()} and \code{c()}.
#' @param yind a vector with length equal to the number of columns of the matrix
#'   of functional responses giving the vector of evaluation points \eqn{(t_1,
#'   \dots ,t_{G})}. If not supplied, \code{yind} is set to
#'   \code{1:ncol(<response>)}.
#' @param algorithm the name of the function used to estimate the model.
#'   Defaults to \code{\link[mgcv]{gam}} if the matrix of functional responses
#'   has less than \code{2e5} data points and to \code{\link[mgcv]{bam}} if not.
#'   \code{'\link[mgcv]{gamm}'}, \code{'\link[gamm4]{gamm4}'} and
#'   \code{'\link[mgcv]{jagam}'} are valid options as well. See Details for
#'   \code{'\link[gamm4]{gamm4}'} and \code{'\link[mgcv]{jagam}'}.
#' @param data an (optional) \code{data.frame} containing the data. Can also be
#'   a named list for regular data. Functional covariates have to be supplied as
#'   <no. of observations> by <no. of evaluations> matrices, i.e. each row is
#'   one functional observation.
#' @param ydata an (optional) \code{data.frame} supplying functional responses
#'   that are not observed on a regular grid. See Details.
#' @param method Defaults to \code{"REML"}-estimation, including of unknown
#'   scale. If \code{algorithm="bam"}, the default is switched to
#'   \code{"fREML"}. See \code{\link[mgcv]{gam}} and \code{\link[mgcv]{bam}} for
#'   details.
#' @param bs.yindex a named (!) list giving the parameters for spline bases on
#'   the index of the functional response. Defaults to \code{list(bs="ps", k=5,
#'   m=c(2, 1))}, i.e. 5 cubic B-splines bases with first order difference
#'   penalty.
#' @param bs.int a named (!) list giving the parameters for the spline basis for
#'   the global functional intercept. Defaults to \code{list(bs="ps", k=20,
#'   m=c(2, 1))}, i.e. 20 cubic B-splines bases with first order difference
#'   penalty.
#' @param tensortype which typ of tensor product splines to use. One of
#'   "\code{\link[mgcv]{ti}}" or "\code{\link[mgcv]{t2}}", defaults to
#'   \code{ti}. \code{t2}-type terms do not enforce the more suitable special
#'   constraints for functional regression, see Details.
#' @param ... additional arguments that are valid for \code{\link[mgcv]{gam}},
#'   \code{\link[mgcv]{bam}}, \code{'\link[gamm4]{gamm4}'} or
#'   \code{'\link[mgcv]{jagam}'}. \code{subset} is not implemented.
#' @return A fitted \code{pffr}-object, which is a
#'   \code{\link[mgcv]{gam}}-object with some additional information in an
#'   \code{pffr}-entry. If \code{algorithm} is \code{"gamm"} or \code{"gamm4"},
#'   only the \code{$gam} part of the returned list is modified in this way.\cr
#'   Available methods/functions to postprocess fitted models:
#'   \code{\link{summary.pffr}}, \code{\link{plot.pffr}},
#'   \code{\link{coef.pffr}}, \code{\link{fitted.pffr}},
#'   \code{\link{residuals.pffr}}, \code{\link{predict.pffr}},
#'   \code{\link{model.matrix.pffr}},  \code{\link{qq.pffr}},
#'   \code{\link{pffr.check}}.\cr If \code{algorithm} is \code{"jagam"}, only
#'   the location of the model file and the usual
#'   \code{\link[mgcv]{jagam}}-object are returned, you have to run the sampler
#'   yourself.\cr
#' @author Fabian Scheipl, Sonja Greven
#' @seealso \code{\link[mgcv]{smooth.terms}} for details of \code{mgcv} syntax
#'   and available spline bases and penalties.
#' @references Ivanescu, A., Staicu, A.-M., Scheipl, F. and Greven, S. (2015).
#'   Penalized function-on-function regression. Computational Statistics,
#'   30(2):539--568. \url{https://biostats.bepress.com/jhubiostat/paper254/}
#'
#'   Scheipl, F., Staicu, A.-M. and Greven, S. (2015). Functional Additive Mixed
#'   Models. Journal of Computational & Graphical Statistics, 24(2): 477--501.
#'   \url{ https://arxiv.org/abs/1207.5947}
#'
#'   F. Scheipl, J. Gertheiss, S. Greven (2016):  Generalized Functional Additive Mixed Models,
#'   Electronic Journal of Statistics, 10(1), 1455--1492.
#'   \url{https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-10/issue-1/Generalized-functional-additive-mixed-models/10.1214/16-EJS1145.full}
#' @export
#' @importFrom mgcv ti jagam gam gam.fit bam gamm
#' @importFrom gamm4 gamm4
#' @importFrom lme4 lmer
#' @examples
#' ###############################################################################
#' # univariate model:
#' # Y(t) = f(t)  + \int X1(s)\beta(s,t)ds + eps
#' set.seed(2121)
#' data1 <- pffrSim(scenario="ff", n=40)
#' t <- attr(data1, "yindex")
#' s <- attr(data1, "xindex")
#' m1 <- pffr(Y ~ ff(X1, xind=s), yind=t, data=data1)
#' summary(m1)
#' plot(m1, pages=1)
#'
#' \dontrun{
#' ###############################################################################
#' # multivariate model:
#' # E(Y(t)) = \beta_0(t)  + \int X1(s)\beta_1(s,t)ds + xlin \beta_3(t) +
#' #        f_1(xte1, xte2) + f_2(xsmoo, t) + \beta_4 xconst
#' data2 <- pffrSim(scenario="all", n=200)
#' t <- attr(data2, "yindex")
#' s <- attr(data2, "xindex")
#' m2 <- pffr(Y ~  ff(X1, xind=s) + #linear function-on-function
#'                 xlin  +  #varying coefficient term
#'                 c(te(xte1, xte2)) + #bivariate smooth term in xte1 & xte2, const. over Y-index
#'                 s(xsmoo) + #smooth effect of xsmoo varying over Y-index
#'                 c(xconst), # linear effect of xconst constant over Y-index
#'         yind=t,
#'         data=data2)
#' summary(m2)
#' plot(m2)
#' str(coef(m2))
#' # convenience functions:
#' preddata <- pffrSim(scenario="all", n=20)
#' str(predict(m2, newdata=preddata))
#' str(predict(m2, type="terms"))
#' cm2 <- coef(m2)
#' cm2$pterms
#' str(cm2$smterms, 2)
#' str(cm2$smterms[["s(xsmoo)"]]$coef)
#'
#' #############################################################################
#' # sparse data (80% missing on a regular grid):
#' set.seed(88182004)
#' data3 <- pffrSim(scenario=c("int", "smoo"), n=100, propmissing=0.8)
#' t <- attr(data3, "yindex")
#' m3.sparse <- pffr(Y ~ s(xsmoo), data=data3$data, ydata=data3$ydata, yind=t)
#' summary(m3.sparse)
#' plot(m3.sparse,pages=1)
#' }
pffr <- function(
  formula,
  yind,
  data=NULL,
  ydata=NULL,
  algorithm = NA,
  method="REML",
  tensortype = c("ti", "t2"),
  bs.yindex = list(bs="ps", k=5, m=c(2, 1)), # only bs, k, m are propagated...
  bs.int = list(bs="ps", k=20, m=c(2, 1)), # only bs, k, m are propagated...
  ...
){
  # TODO: subset args!

  call <- match.call()
  tensortype <- as.symbol(match.arg(tensortype))
  # make sure we use values for the args that were defined as close to the
  #  actual function call as possible:
  lapply(names(head(call, -1))[-1], function(nm)
    try(assign(nm, eval(nm, parent.frame()))))
  ## TODO: does this make sense? useful if pffr is called from a function that
  ## supplies args as variables that are also defined differently in GlobalEnv:
  ## this should then ensure that the args as defined in the calling function,
  ## and not in GlobalEnv get used....

  ## warn if entries in ... aren't arguments for gam/gam.fit/jagam or gamm4/lmer
  ## check for special case of gaulss family
  dots <- list(...)
  gaulss <- FALSE
  if(length(dots)){
    validDots <- if(!is.na(algorithm) && algorithm=="gamm4"){
      c(names(formals(gamm4)), names(formals(lmer)))
    } else {
      c(names(formals(gam)), names(formals(bam)), names(formals(gam.fit)),
        names(formals(jagam)))
    }
    if(!is.null(dots$family) && dots$family == "gaulss") {
      validDots <- c(validDots, "varformula")
      gaulss <- TRUE
    }
    notUsed <- names(dots)[!(names(dots) %in% validDots)]
    if(length(notUsed))
      warning("Arguments <", paste(notUsed, collapse=", "), "> supplied but not used." )
  }




  sparseOrNongrid <- !is.null(ydata)
  if(sparseOrNongrid){
    stopifnot(ncol(ydata)==3)
    stopifnot(c(".obs", ".index", ".value") == colnames(ydata))
  }

  pffrspecials <- c("s", "te", "ti", "t2", "ff", "c", "sff", "ffpc", "pcre")
  tf <- terms.formula(formula, specials=pffrspecials)
  trmstrings <- attr(tf, "term.labels")
  terms <- sapply(trmstrings, function(trm) as.call(parse(text=trm))[[1]], simplify=FALSE)
  #ugly, but getTerms(formula)[-1] does not work for terms like I(x1:x2)
  frmlenv <- environment(formula)

  #which terms are which type:
  where.specials <- sapply(pffrspecials, function(sp) attr(tf, "specials")[[sp]]-1)
  if(length(trmstrings)) {
    where.specials$par <- which(!(1:length(trmstrings) %in%
                                    (unlist(attr(tf, "specials")) - 1)))
    # indices of linear/factor terms with functional coefficients over yind
  } else where.specials$par <- numeric(0)

  responsename <- attr(tf,"variables")[2][[1]]

  #start new formula
  newfrml <- paste(responsename, "~", sep="")
  newfrmlenv <- new.env()
  evalenv <- if("data" %in% names(call)) eval.parent(call$data) else NULL

  if(sparseOrNongrid){
    nobs <- length(unique(ydata$.obs))
    stopifnot(all(ydata$.obs %in% rownames(data)))
    # FIXME: allow for non-1:nobs .obs-formats!
    stopifnot(all(ydata$.obs %in% 1:nobs))


    #works for data-lists or matrix-valued covariates as well:
    nobs.data <- nrow(as.matrix(data[[1]]))
    stopifnot(nobs == nobs.data)
    ntotal <- nrow(ydata)

    #generate yind for estimates/predictions etc
    yind <- if(length(unique(ydata$.index))>100){
      seq(min(ydata$.index), max(ydata$.index), l=100)
    } else {
      sort(unique(ydata$.index))
    }
    nyindex <- length(yind)
  } else {
    nobs <- nrow(eval(responsename,  envir=evalenv, enclos=frmlenv))
    nyindex <- ncol(eval(responsename,  envir=evalenv, enclos=frmlenv))
    ntotal <- nobs*nyindex
  }

  if(missing(algorithm)||is.na(algorithm)){
    algorithm <- ifelse(ntotal > 1e5, "bam", "gam")
  }
  if(algorithm == "bam" & missing(method)){
      call$method <- "fREML"
  }
  algorithm <- as.symbol(algorithm)
  if(as.character(algorithm)=="bam" && !("chunk.size" %in% names(call))){
    call$chunk.size <- 10000
    #same default as in bam
  }
  ## no te-terms possible in gamm4:
  if(as.character(algorithm)=="gamm4"){
    stopifnot(length(unlist(where.specials[c("te","ti")]))<1)
  }


  if(!sparseOrNongrid){
    #if missing, define y-index or get it from first ff/sff-term, then assign expanded versions to newfrmlenv
    if(missing(yind)){
      if(length(c(where.specials$ff, where.specials$sff))){
        if(length(where.specials$ff)){
          ffcall <- expand.call(ff, as.call(terms[where.specials$ff][1])[[1]])
        }  else ffcall <- expand.call(sff, as.call(terms[where.specials$sff][1])[[1]])
        if(!is.null(ffcall$yind)){
          yind <- eval(ffcall$yind, envir=evalenv, enclos=frmlenv)
          yindname <- deparse(ffcall$yind)
        } else {
          yind <- 1:nyindex
          yindname <- "yindex"
        }
      } else {
        yind <- 1:nyindex
        yindname <- "yindex"
      }
    } else {
      if (is.symbol(substitute(yind)) | is.character(yind)) {
        yindname <- deparse(substitute(yind))
        if(!is.null(data) && !is.null(data[[yindname]])){
          yind <- data[[yindname]]
        }
      } else {
        yindname <- "yindex"
      }
      stopifnot(is.vector(yind), is.numeric(yind),
                length(yind) == nyindex)
    }
    #make sure it's a valid name
    if(length(yindname)>1) yindname <- "yindex"
    # make sure yind is sorted
    stopifnot(all.equal(order(yind), 1:nyindex))


    yindvec <- rep(yind, times = nobs)
    yindvecname <- as.symbol(paste(yindname,".vec",sep=""))
    assign(x=deparse(yindvecname), value=yindvec, envir=newfrmlenv)

    #assign response in _long_ format to newfrmlenv
    assign(x=deparse(responsename), value=as.vector(t(eval(responsename, envir=evalenv,
                                                           enclos=frmlenv))),
           envir=newfrmlenv)


    missingind <-  if(any(is.na(get(as.character(responsename), newfrmlenv)))){
      which(is.na(get(as.character(responsename), newfrmlenv)))
    } else NULL

    # repeat which row in <data> how many times
    stackpattern <- rep(1:nobs, each=nyindex)

  } else {
    # sparseOrNongrid:
    yindname <- "yindex"
    yindvec <- ydata$.index
    yindvecname <- as.symbol(paste(yindname,".vec",sep=""))
    assign(x=deparse(yindvecname), value=ydata$.index, envir=newfrmlenv)

    #assign response in _long_ format to newfrmlenv
    assign(x=deparse(responsename), value=ydata$.value, envir=newfrmlenv)

    missingind <- NULL

    # repeat which row in <data> how many times:
    stackpattern <- ydata$.obs
  }


  ##################################################################################
  #modify formula terms....
  newtrmstrings <- attr(tf, "term.labels")

  #if intercept, add \mu(yindex)
  if(attr(tf, "intercept")){
    # have to jump thru some hoops to get bs.yindex handed over properly
    # without having yind evaluated within the call
    arglist <- c(name="s", x = as.symbol(yindvecname), bs.int)
    intcall <- NULL
    assign(x= "intcall", value= do.call("call", arglist, envir=newfrmlenv), envir=newfrmlenv)
    newfrmlenv$intcall$x <- as.symbol(yindvecname)

    intstring <- deparse(newfrmlenv$intcall)
    rm(intcall, envir=newfrmlenv)

    newfrml <- paste(newfrml, intstring, sep=" ")
    addFint <- TRUE
    names(intstring) <- paste("Intercept(",yindname,")",sep="")
  } else{
    newfrml <-paste(newfrml, "0", sep="")
    addFint <- FALSE
  }

  #transform: c(foo) --> foo
  if(length(where.specials$c)){
    newtrmstrings[where.specials$c] <- sapply(trmstrings[where.specials$c], function(x){
      sub("\\)$", "", sub("^c\\(", "", x)) #c(BLA) --> BLA
    })
  }

  #prep function-on-function-terms
  if(length(c(where.specials$ff, where.specials$sff))){
    ffterms <- lapply(terms[c(where.specials$ff, where.specials$sff)], function(x){
      eval(x, envir=evalenv, enclos=frmlenv)
    })

    newtrmstrings[c(where.specials$ff, where.specials$sff)] <- sapply(ffterms, function(x) {
      safeDeparse(x$call)
    })

    #apply limits function and assign stacked data to newfrmlenv
    makeff <- function(x){
      tmat <- matrix(yindvec, nrow=length(yindvec), ncol=length(x$xind))
      smat <- matrix(x$xind, nrow=length(yindvec), ncol=length(x$xind),
                     byrow=TRUE)
      if(!is.null(x[["LX"]])){
        # for ff: stack weights * covariate
        LStacked <- x$LX[stackpattern,]
      } else {
        # for sff: stack weights, X separately
        LStacked <- x$L[stackpattern,]
        XStacked <- x$X[stackpattern, ]
      }
      if(!is.null(x$limits)){
        # find int-limits and set weights to 0 outside
        use <- x$limits(smat, tmat)
        LStacked <- LStacked * use

        # find indices for row-wise int-range & maximal occuring width:
        windows <- t(apply(use, 1, function(x){
          use_this <- which(x)
          # edge case: no integration
          if(!any(use_this)) return(c(1,1))
          range(use_this)
        }))
        windows <- cbind(windows, windows[,2]-windows[,1]+1)
        maxwidth <- max(windows[,3])
        # reduce size of matrix-covariates if possible:
        if(maxwidth < ncol(smat)){
          # all windows have to have same length, so modify windows:
          eff.windows <- t(apply(windows, 1, function(window,
                                                      maxw=maxwidth,
                                                      maxind=ncol(smat)){
            width <- window[3]
            if((window[2] + maxw - width) <= maxind){
              window[1] : (window[2] + maxw -width)
            } else {
              (window[1] + width - maxw) : window[2]
            }
          }))

          # extract relevant parts of each row and stack'em
          shift_and_shorten <- function(X, eff.windows){
            t(sapply(1:(nrow(X)),
                     function(i) X[i, eff.windows[i,]]))
          }
          smat <- shift_and_shorten(smat, eff.windows)
          tmat <- shift_and_shorten(tmat, eff.windows)
          LStacked <- shift_and_shorten(LStacked, eff.windows)
          if(is.null(x$LX)){ # sff
            XStacked <- shift_and_shorten(XStacked, eff.windows)
          }
        }
      }
      assign(x=x$yindname,
             value=tmat,
             envir=newfrmlenv)
      assign(x=x$xindname,
             value=smat,
             envir=newfrmlenv)

      assign(x=x$LXname,
             value=LStacked,
             envir=newfrmlenv)
      if(is.null(x[["LX"]])){ # sff
        assign(x=x$xname,
               value=XStacked,
               envir=newfrmlenv)
      }
      invisible(NULL)
    }
    lapply(ffterms, makeff)
  } else ffterms <- NULL

  if(length(where.specials$ffpc)){ ##TODO for sparse
    ffpcterms <- lapply(terms[where.specials$ffpc], function(x){
      eval(x, envir=evalenv, enclos=frmlenv)
    })

    lapply(ffpcterms, function(trm){
      lapply(colnames(trm$data), function(nm){
        assign(x=nm, value=trm$data[stackpattern, nm], envir=newfrmlenv)
        invisible(NULL)
      })
      invisible(NULL)
    })


    getFfpcFormula <- function(trm) {
      frmls <- lapply(colnames(trm$data), function(pc) {
        arglist <- c(name="s", x = as.symbol(yindvecname), by= as.symbol(pc),
                     id=trm$id, trm$splinepars)
        call <- do.call("call", arglist, envir=newfrmlenv)
        call$x <- as.symbol(yindvecname)
        call$by <- as.symbol(pc)
        safeDeparse(call)
      })
      return(paste(unlist(frmls), collapse=" + "))
    }
    newtrmstrings[where.specials$ffpc] <- sapply(ffpcterms, getFfpcFormula)

    ffpcterms <- lapply(ffpcterms, function(x) x[names(x)!="data"])
  } else ffpcterms <- NULL

  #prep PC-based random effects
  if(length(where.specials$pcre)){
    pcreterms <- lapply(terms[where.specials$pcre], function(x){
      eval(x, envir=evalenv, enclos=frmlenv)
    })
    #assign newly created data to newfrmlenv
    lapply(pcreterms, function(trm){
      if(!sparseOrNongrid && all(trm$yind==yind)){
        lapply(colnames(trm$efunctions), function(nm){
          assign(x=nm, value=trm$efunctions[rep(1:nyindex, times=nobs), nm],
                 envir=newfrmlenv)
          invisible(NULL)
        })
      } else {
        # don't ever extrapolate eigenfunctions:
        stopifnot(min(trm$yind)<=min(yind))
        stopifnot(max(trm$yind)>=max(yind))

        # interpolate given eigenfunctions to observed index values:
        lapply(colnames(trm$efunctions), function(nm){
          tmp <- approx(x=trm$yind,
                        y=trm$efunctions[, nm],
                        xout=yindvec,
                        method = "linear")$y
          assign(x=nm, value=tmp,
                 envir=newfrmlenv)
          invisible(NULL)
        })
      }
      assign(x=trm$idname, value=trm$id[stackpattern], envir=newfrmlenv)
      invisible(NULL)
    })

    newtrmstrings[where.specials$pcre] <- sapply(pcreterms, function(x) {
      safeDeparse(x$call)
    })
  }else pcreterms <- NULL


  #transform: s(x, ...), te(x, z,...), t2(x, z, ...) --> <ti|t2>(x, <z,> yindex, ..., <bs.yindex>)
  makeSTeT2 <- function(x){

    xnew <- x
    if(deparse(x[[1]]) %in%  c("te", "ti") && as.character(algorithm) == "gamm4") xnew[[1]] <- quote(t2)
    if(deparse(x[[1]]) == "s"){
      xnew[[1]] <- if(as.character(algorithm) != "gamm4") {
        tensortype
      } else quote(t2)
      #accomodate multivariate s()-terms
      xnew$d <- if(!is.null(names(xnew))){
        c(length(all.vars(xnew[names(xnew)==""])), 1)
      } else c(length(all.vars(xnew)), 1)
    } else {
      if("d" %in% names(x)){ #either expand given d...
        xnew$d <- c(eval(x$d), 1)
      } else {#.. or default to univariate marginal bases
        xnew$d <- rep(1, length(all.vars(x))+1)
      }
    }
    xnew[[length(xnew)+1]] <- yindvecname
    this.bs.yindex <- if("bs.yindex" %in% names(x)){
      eval(x$bs.yindex)
    } else bs.yindex
    xnew <- xnew[names(xnew) != "bs.yindex"]

    if(deparse(xnew[[1]]) == "ti"){
      # apply sum-to-zero constraints to marginal bases for covariate(s),
      #  but not to <yindex> to get terms with sum-to-zero-for-each-t constraints
      xnew$mc <- c(rep(TRUE, length(xnew$d)-1), FALSE)
    }

    xnew$bs <- if("bs" %in% names(x)){
      if("bs" %in% names(this.bs.yindex)){
        c(eval(x$bs), this.bs.yindex$bs)
      } else {
        c(xnew$bs, "tp")
      }
    } else {
      if("bs" %in% names(this.bs.yindex)){
        c(rep("tp", length(xnew$d)-1), this.bs.yindex$bs)
      } else {
        rep("tp", length(all.vars(x))+1)
      }
    }
    xnew$m <- if("m" %in% names(x)){
      if("m" %in% names(this.bs.yindex)){
        warning("overriding bs.yindex for m in ", deparse(x))
      }
      #TODO: adjust length if necessary, m can be a list for bs="ps","cp","ds"!
      x$m
    } else {
      if("m" %in% names(this.bs.yindex)){
        this.bs.yindex$m
      } else {
        NA
      }
    }
    #defaults to 8 basis functions
    xnew$k <- if("k" %in% names(x)){
      if("k" %in% names(this.bs.yindex)){
        c(eval(xnew$k), eval(this.bs.yindex$k))
      } else {
        c(eval(xnew$k), 8)
      }
    } else {
      if("k" %in% names(this.bs.yindex)){
        c(pmax(8, 5^head(eval(xnew$d), -1)), eval(this.bs.yindex$k))
      } else {
        pmax(8, 5^eval(xnew$d))
      }
    }
    xnew$k <- unlist(xnew$k)

    if("xt" %in% names(x)){
#       # xt has to be supplied as a list, with length(x$d) entries,
#       # each of which is a list or NULL:
#       stopifnot(x$xt[[1]]==as.symbol("list") &&
#                   # =length(x$d)+1, since first element in parse tree is 'list'
#                   all(sapply(2:length(x$xt), function(i)
#                     x$xt[[i]][[1]] == as.symbol("list") ||
#                       is.null(eval(x$xt[[i]][[1]])))))
      xnew$xt <- x$xt
    }

    ret <- safeDeparse(xnew)
    return(ret)
  }

  if(length(c(where.specials$s, where.specials$te, where.specials$t2))){
    newtrmstrings[c(where.specials$s, where.specials$te, where.specials$t2)] <-
      sapply(terms[c(where.specials$s, where.specials$te, where.specials$t2)],
             makeSTeT2)
  }

  #transform: x --> s(YINDEX, by=x)
  if(length(where.specials$par)){
    newtrmstrings[where.specials$par] <- sapply(terms[where.specials$par], function(x){
      xnew <- bs.yindex
      xnew <- as.call(c(quote(s), yindvecname, by=x, xnew))
      safeDeparse(xnew)
    })

  }

  #... & assign expanded/additional variables to newfrmlenv
  where.specials$notff <- c(where.specials$c, where.specials$par,
                            where.specials$s, where.specials$te, where.specials$t2)
  if(length(where.specials$notff)){
    # evalenv below used to be list2env(eval.parent(call$data)), frmlenv),
    # but that assigned everything in <data> to the global workspace if frmlenv was the global
    # workspace.
    evalenv <- if("data" %in% names(call))  {
      list2env(eval.parent(call$data))
    } else frmlenv
    lapply(terms[where.specials$notff],
           function(x){
             #nms <- all.vars(x)
             isC <- safeDeparse(x) %in% sapply(terms[where.specials$c], safeDeparse)
             if(isC) {
               # drop c()
               # FIXME: FUGLY!
               x <- formula(paste("~", gsub("\\)$", "",
                                            gsub("^c\\(", "", deparse(x)))))[[2]]
             }
             ## remove names in xt, k, bs,  information (such as variable names for MRF penalties etc)
             nms <- if(!is.null(names(x))){
               all.vars(x[names(x) %in% c("", "by")])
             }  else all.vars(x)


             sapply(nms, function(nm){
               var <- get(nm, envir=evalenv)
               if(is.matrix(var)){
                   stopifnot(!sparseOrNongrid || ncol(var) == nyindex)
                   assign(x=nm,
                          value=as.vector(t(var)),
                          envir=newfrmlenv)
               } else {
                   stopifnot(length(var) == nobs)
                   assign(x=nm,
                          value=var[stackpattern],
                          envir=newfrmlenv)
               }
               invisible(NULL)
             })
             invisible(NULL)
           })
  }

  newfrml <- formula(paste(c(newfrml, newtrmstrings), collapse="+"))
  environment(newfrml) <- newfrmlenv

  # variance formula for gaulss
  if(gaulss) {
    if(is.null(dots$varformula)) {
        dots$varformula <- formula(paste("~", safeDeparse(
          as.call(c(as.name("s"), x = as.symbol(yindvecname), bs.int)))))
    }
    environment(dots$varformula) <- newfrmlenv
    newfrml <- list(newfrml, dots$varformula)
  }
  pffrdata <- list2df(as.list(newfrmlenv))

  newcall <- expand.call(pffr, call)
  newcall$yind <- newcall$tensortype <- newcall$bs.int <-
    newcall$bs.yindex <- newcall$algorithm <- newcall$ydata <- NULL
  newcall$formula <- newfrml
  newcall$data <- quote(pffrdata)
  newcall[[1]] <- algorithm
  # make sure ...-args are taken from ..., not GlobalEnv:
  dotargs <- names(newcall)[names(newcall) %in% names(dots)]
  newcall[dotargs] <- dots[dotargs]
  if("subset" %in% dotargs){
    stop("<subset>-argument is not supported.")
  }
  if("weights" %in% dotargs){
    wtsdone <- FALSE
    if(length(dots$weights) == nobs){
        newcall$weights <- dots$weights[stackpattern]
        wtsdone <- TRUE
    }
    if (!is.null(dim(dots$weights)) &&
       all(dim(dots$weights) == c(nobs, nyindex))) {
        newcall$weights <- as.vector(t(dots$weights))
        wtsdone <- TRUE
    }
    if(!wtsdone){
        stop("weights have to be supplied as a vector with length=rows(data) or
               a matrix with the same dimensions as the response.")
    }
  }
  if("offset" %in% dotargs){
    ofstdone <- FALSE
    if(length(dots$offset) == nobs){
        newcall$offset <- dots$offset[stackpattern]
        ofstdone <- TRUE
    }
    if(!is.null(dim(dots$offset)) &&
       all(dim(dots$offset) == c(nobs, nyindex))){
        newcall$offset <- as.vector(t(dots$offset))
        ofstdone <- TRUE
    }
    if(!ofstdone){
        stop("offsets have to be supplied as a vector with length=rows(data) or
             a matrix with the same dimensions as the response.")
    }
  }
  if(as.character(algorithm) == "jagam"){
    newcall <- newcall[names(newcall) %in% c("", names(formals(jagam)))]
    if(is.null(newcall$file)) {
      newcall$file <- tempfile("pffr2jagam", tmpdir = getwd(), fileext = ".jags")
    }
  }
  # call algorithm to estimate model
  m <- eval(newcall)
  if(as.character(algorithm) == "jagam"){
    m$modelfile <- newcall$file
    message("JAGS/BUGS model code written to \n", m$modelfile, ",\n see ?jagam")
    return(m)
  }

  m.smooth <- if(as.character(algorithm) %in% c("gamm4","gamm")){
    m$gam$smooth
  } else m$smooth

  #return some more info s.t. custom predict/plot/summary will work
  trmmap <- newtrmstrings
  names(trmmap) <- names(terms)
  if(addFint) trmmap <- c(trmmap, intstring)

  # map labels to terms --
  # ffpc are associated with multiple smooths
  # parametric are associated with multiple smooths if covariate is a factor
  labelmap <- as.list(trmmap)
  lbls <- sapply(m.smooth, function(x) x$label)
  if(length(c(where.specials$par, where.specials$ffpc))){
    if(length(where.specials$par)){
      for(w in where.specials$par){
        # only combine if <by>-variable is a factor!
        if(is.factor(get(names(labelmap)[w], envir=newfrmlenv))){
          labelmap[[w]] <- {
            #covariates for parametric terms become by-variables:
            where <- sapply(m.smooth, function(x) x$by) == names(labelmap)[w]
            sapply(m.smooth[where], function(x) x$label)
          }
        } else {
          labelmap[[w]] <- paste0("s(",yindvecname,"):",names(labelmap)[w])
        }
      }
    }
    if(length(where.specials$ffpc)){
      ind <- 1
      for(w in where.specials$ffpc){
        labelmap[[w]] <- {
          #PCs for X become by-variables:
          where <- sapply(m.smooth, function(x) x$id) == ffpcterms[[ind]]$id
          sapply(m.smooth[where], function(x) x$label)
        }
        ind <- ind+1
      }
    }
    labelmap[-c(where.specials$par, where.specials$ffpc)] <- lbls[pmatch(
      sapply(labelmap[-c(where.specials$par, where.specials$ffpc)], function(x){
        ## FUGLY: check whether x is a function call of some sort
        ##  or simply a variable name.
        if(length(parse(text=x)[[1]]) != 1){
          tmp <- eval(parse(text=x))
          return(tmp$label)
        } else {
          return(x)
        }
      }), lbls)]
  } else{
    labelmap[1:length(labelmap)] <-  lbls[pmatch(
      sapply(labelmap[1:length(labelmap)], function(x){
        ## FUGLY: check whether x is a function call of some sort
        ##  or simply a variable name.
        if(length(parse(text=x)[[1]]) != 1){
          tmp <- eval(parse(text=x))
          return(tmp$label)
        } else {
          return(x)
        }
      }), lbls)]

  }
  # check whether any parametric terms were left out & add them
  nalbls <- sapply(labelmap,
                   function(x) {
                     any(is.null(x)) | any(is.na(x[!is.null(x)]))
                   })
  if (any(nalbls)) {
    labelmap[nalbls] <- trmmap[nalbls]
  }

  names(m.smooth) <- lbls
  if(as.character(algorithm) %in% c("gamm4","gamm")){
    m$gam$smooth <- m.smooth
  } else{
    m$smooth  <- m.smooth
  }

  ret <-  list(
    call=call,
    formula=formula,
    termmap=trmmap,
    labelmap=labelmap,
    responsename = responsename,
    nobs=nobs,
    nyindex=nyindex,
    yindname = yindname,
    yind=yind,
    where=where.specials,
    ff=ffterms,
    ffpc=ffpcterms,
    pcreterms=pcreterms,
    missingind = missingind,
    sparseOrNongrid=sparseOrNongrid,
    ydata=ydata)

  if(as.character(algorithm) %in% c("gamm4","gamm")){
    m$gam$pffr <- ret
    class(m$gam) <- c("pffr", class(m$gam))
  } else {
    m$pffr <- ret
    class(m) <- c("pffr", class(m))
  }
  return(m)
}# end pffr()
