
#### see http://adv-r.had.co.nz/S3.html
#### for the advices on best practices for a S3 class  

######### Create a construction method that checks the types of the input, 
# and returns a list with the correct class label. XXX <- function(...) {}

#' A S3 class for univariate functional data on a common grid
#' 
#' The hmatrix class represents data for a functional historical effect. 
#' The class is basically a matrix containing the time and the id for the observations of the 
#' functional response. The functional covariate is contained as attribute. 
#' @param time set of argument values of the response in long format, 
#' i.e. at which \code{t} the response curve is observed
#' @param id specify to which curve the point belongs to, id from 1, 2, ..., n.  
#' @param x matrix of functional covariate, each trajectory is in one row 
#' @param argvals set of argument values, i.e., the common gird at which the functional covariate 
#' is observed, by default \code{1:ncol(x)}
#' @param timeLab name of the time axis, by default \code{t}
#' @param idLab name of the id variable, by default \code{wideIndex}
#' @param xLab name of the functional variable, by default NULL
#' @param argvalsLab name of the argument for the covariate by default \code{s}
#' 
#' @details In the hmatrix class the id has to run from i=1, 2, ..., n including all integers from 1 to n. 
#' The rows of the functional covariate x correspond to those observations. 
#' 
#' @seealso \code{\link{getTime.hmatrix}} to extract attributes, 
#' and ?"[.hmatrix" for the extract method. 
#'
#' @examples 
#' ## Example for a hmatrix object
#' t1 <- rep((1:5)/2, each = 3)
#' id1 <- rep(1:3, 5)
#' x1 <- matrix(1:15, ncol = 5) 
#' s1 <- (1:5)/2 
#' myhmatrix <- hmatrix(time = t1, id = id1, x = x1, argvals = s1, 
#'                      timeLab = "t1", argvalsLab = "s1", xLab = "test")
#' 
#' # extract with [ keeps attributes 
#' # select observations of subjects 2 and 3
#' myhmatrixSub <- myhmatrix[id1 %in% c(2, 3), ]  
#' str(myhmatrixSub)
#' getX(myhmatrixSub)
#' getX(myhmatrix)
#' 
#' # get time
#' myhmatrix[ , 1] # as column matrix as drop = FALSE
#' getTime(myhmatrix) # as vector
#' 
#' # get id
#' myhmatrix[ , 2] # as column matrix as drop = FALSE
#' getId(myhmatrix) # as vector
#' 
#' # subset hmatrix on the basis of an index, which is defined on the curve level
#' reweightData(data = list(hmat = myhmatrix), vars = "hmat", index = c(1, 1, 2))
#' # this keeps only the unique x values in attr(,'x') but multiplies the corresponding
#' # ids and times in the time id matrix 
#' # for bhistx baselearner, there may be an additional id variable for the tensor product
#' newdat <- reweightData(data = list(hmat = myhmatrix, 
#'   repIDx = rep(1:nrow(attr(myhmatrix,'x')), length(attr(myhmatrix,"argvals")))), 
#'   vars = "hmat", index = c(1,1,2), idvars="repIDx")
#' length(newdat$repIDx) 
#' 
#' ## use hmatrix within a data.frame
#' mydat <- data.frame(I(myhmatrix), z=rnorm(3)[id1])
#' str(mydat)
#' str(mydat[id1 %in% c(2, 3), ])
#' str(myhmatrix[id1 %in% c(2, 3), ])
#'  
#' @export
hmatrix <- function(time, id, x, argvals=1:ncol(x), 
                    timeLab="t", idLab="wideIndex", xLab="x", argvalsLab="s"){
   
  ## check that id is integer valued containing 1, 2, 3, ..., n 
  ## and that x has n rows
  stopifnot( all(sort(unique(id)) == 1:nrow(x)) )  
  stopifnot(length(time)==length(id))
    
  # convert x to a matrix, especially if x is of class AsIs
  x <- matrix(x, ncol=ncol(x), nrow=nrow(x))  
   
  #### check argvals and x
  if( any(duplicated(argvals)) ){
    stop("argvals contains duplicates.")
  } 
  if( is.unsorted(argvals) ){
    stop("argvals is not sorted.")
  }
  
  if (ncol(x)!=length(argvals)) {
    stop(quote(x), " must have same number of columns as the length of ", quote(s), ".")
  } 
  
  ret <- matrix(c(time, id), ncol=2)
  colnames(ret) <- c("time","id")
  ## ret <- data.frame(time=time, id=id) # use matrix to use hmatrix within a data.frame
  
  attr(ret, "x") <- x
  attr(ret, "argvals") <- argvals
  attr(ret, "timeLab") <- timeLab
  attr(ret, "idLab") <- idLab 
  attr(ret, "xLab") <- xLab
  attr(ret, "argvalsLab") <- argvalsLab
  class(ret) <- c("hmatrix", class(ret)) 
  ret  
}


### Define the generic methods
#' Generic functions to asses attributes of functional data objects
#' 
#' Extract attributes of an object.  
#' @param object an R-object, currently implemented for hmatrix and fmatrix
#' 
#' @details Extract the time variable \code{getTime}, the id\code{getId}, 
#' the functional covariate \code{getX}, its argument values \code{getArgvals}. 
#' Or the names of the different variables \code{getTimeLab}, 
#' \code{getIdLab}, \code{getXLab}, \code{getArgvalsLab}. 
#' 
#' @seealso \code{\link{hmatrix}} for the h.atrix class. 
#' 
#' @aliases getId getX getArgvals getTimeLab getIdLab getXLab getArgvalsLab
#'
#' @export
getTime <- function(object) { UseMethod("getTime", object) }

#' @rdname getTime
#' @export
getId <- function(object) { UseMethod("getId", object) }

#' @rdname getTime
#' @export
getX <- function(object) { UseMethod("getX", object) }

#' @rdname getTime
#' @export
getArgvals <- function(object) { UseMethod("getArgvals", object) }

#' @rdname getTime
#' @export
getTimeLab <- function(object) { UseMethod("getTimeLab", object) }

#' @rdname getTime
#' @export
getIdLab <- function(object) { UseMethod("getIdLab", object) }

#' @rdname getTime
#' @export
getXLab <- function(object) { UseMethod("getXLab", object) }

#' @rdname getTime
#' @export
getArgvalsLab <- function(object) { UseMethod("getArgvalsLab", object) }



#' Extract attributes of hmatrix
#' 
#' Extract attributes of an object of class \code{hmatrix}.  
#' @param object object of class hmatrix
#' 
#' @details Extract the time variable \code{getTime}, the id\code{getId}, 
#' the functional covariate \code{getX}, its argument values \code{getArgvals}. 
#' Or the names of the different variables \code{getTimeLab}, 
#' \code{getIdLab}, \code{getXLab}, \code{getArgvalsLab} for an object of class \code{hmatrix}.  
#' 
#' @aliases getId.hmatrix getX.hmatrix getArgvals.hmatrix getTimeLab.hmatrix getXLab.hmatrix getArgvalsLab.hmatrix
#'
#' @export
getTime.hmatrix <- function(object) object[ , 1, drop=TRUE]

#' @rdname getTime.hmatrix
#' @export
getId.hmatrix <- function(object)  object[ , 2, drop=TRUE]

#' @rdname getTime.hmatrix
#' @export
getX.hmatrix <- function(object) attr(object, "x")

#' @rdname getTime.hmatrix
#' @export
getArgvals.hmatrix <- function(object) attr(object, "argvals")

#' @rdname getTime.hmatrix
#' @export
getTimeLab.hmatrix <- function(object) attr(object, "timeLab")

#' @rdname getTime.hmatrix
#' @export
getIdLab.hmatrix <- function(object) attr(object, "idLab")

#' @rdname getTime.hmatrix
#' @export
getXLab.hmatrix <- function(object) attr(object, "xLab")

#' @rdname getTime.hmatrix
#' @export
getArgvalsLab.hmatrix <- function(object) attr(object, "argvalsLab")

######### Write a function to check if an object is of your class: 
# is.XXX <- function(x) inherits(x, "XXX")
#' Test to class of hmatrix
#' 
#' is.hmatrix tests if its argument is an object of class hmatrix.   
#' @param object object of class hmatrix
#'
#' @export
is.hmatrix <- function(object){
  inherits(object, "hmatrix")
}

######### When implementing a vector class, you should implement these methods: 
# length, [, [<-, [[, [[<-, c. (If [ is implemented rev, head, and tail should all work).

#' Extract or replace parts of a hmatrix-object
#' 
#' Operator acting on hmatrix preserving the attributes when rows are extracted.  
#' @param x object from which to extract element(s) or in which to replace element(s).
#' @param i,j indices specifying elements to extract or replace. Indices are numeric 
#' vectors or empty (missing) or NULL. Numeric values are coerced to integer as by as.integer 
#' (and hence truncated towards zero). 
#' @param ... not used
#' @param drop  If \code{TRUE} the result is coerced to the lowest possible dimension 
#' (or just a matrix). This only works for extracting elements, not for the 
#' replacement, defaults to \code{FALSE}.
#' 
#' @details If used on columns or rows/columns a matrix is returned. 
#' If used on rows only, i.e. x[i,] an object of class hmatrix is returned. 
#' The id is changed so that it runs from 1, ..., nNew, where nNew is the number of different 
#' id values in the new hmatrix-object. 
#' From the functional covariate \code{x} rows are selected accordingly.
#'  
#' @seealso ?"["
#'
#' @export 
`[.hmatrix` <- function(x, i, j, ..., drop=FALSE) {
  
  # number of arguments without drop
  Narg <- nargs() - (!missing(drop)) 
  
  # save attributs of x
  xAttr <- attributes(x) 
  
  ## use "[" method as for a matrix
  r <- NextMethod("[", drop=drop)
  class(r) <- class(r)[class(r)!="hmatrix"] 
  
  ## x[i] return column i 
  if(Narg == 2){
    return(r)
  }
  
  ## x[i,] whole rows are selected 
  ## the is.symbol(j) is used if hmatrix is part of a data.frame using I()
  if(missing(j) || is.symbol(j)){ 
    
    tempId <- r[ ,2] # get the id of the corresponding rows
    tempId <- (1:length(unique(tempId)))[factor(tempId)]  # transform the id to 1, 2, 3, ...        
    
    return( hmatrix(time=r[ ,1], id=tempId, 
                  x=xAttr$x[unique(r[ ,2]), , drop=FALSE], argvals = xAttr$argvals, 
                  timeLab = xAttr$timeLab, idLab = xAttr$idLab, xLab = xAttr$xLab, argvalsLab = xAttr$argvalsLab) )
  }
    
  # x[i,j] select on rows and colums, or only columns x[,j]
  return(r) 
}

#' Transform id and time of wide format into long format
#' 
#' Transform id and time from wide format into long format, i.e., time and id are 
#' repeated accordingly so that two vectors of the same length are returned. 
#' @param time the observation points
#' @param id the id for the curve 
#'
#' @export
wide2long <- function(time, id){
  newtime <- rep(time, each=length(unique(id)))
  newid <- rep(id, length(time))
  return(list(time=newtime, id=newid))
}

#' Subsets hmatrix according to an index
#' 
#' @param x hmatix object that should be subsetted 
#' @param index integer vector with (possibly duplicated) indices
#' for each curve to select
#' @param compress logical, defaults to \code{TRUE}. Only used to force a meaningful
#' behaviour of \code{applyFolds} with hmatrix objects when using nested resampling.
#' 
#' @details This methods is primary useful when subsetting repeatedly.  
#' @examples 
#' t1 <- rep((1:5)/2, each = 3)
#' id1 <- rep(1:3, 5)
#' x1 <- matrix(1:15, ncol = 5) 
#' s1 <- (1:5)/2 
#' hmat <- hmatrix(time = t1, id = id1, x = x1, argvals = s1, timeLab = "t1", 
#'                 argvalsLab = "s1", xLab = "test")
#' 
#' index1 <- c(1, 1, 3)
#' index2 <- c(2, 3, 3)
#' resMat <- subset_hmatrix(hmat, index = index1)
#' try(resMat2 <- subset_hmatrix(resMat, index = index2))
#' resMat <- subset_hmatrix(hmat, index = index1, compress = FALSE)
#' try(resMat2 <- subset_hmatrix(resMat, index = index2))
#'
#' @export
subset_hmatrix <- function(x, index, compress = TRUE)
{
  
  ## get attributes
  attrTemp <- attributes(x)
  
  # save time and id variable of hmatrix-object as ordinary matrix
  # otherwise [ on a hmatrix-object behaves unexpectedly 
  tempMat <- cbind(x[, 1], x[, 2])
  
  # create new matrix for results
  resMat <- matrix(ncol=3)
  
  # for all unique time points t do
  for(t in unique(x[, 1])){ 
    
    # check whether the id exists for this time point 
    idInT <- index %in% tempMat[tempMat[,1] == t, 2]
    # add rows for observations selected by index for time t
    resMat <- rbind(resMat, 
                    matrix(c(rep(t, sum(idInT)), # for time points in hmatrix
                             index[idInT], # for id in hmatrix
                             (1:length(index))[idInT]), # for idvars 
                           ncol=3))
    
  }
  
  # drop first row with NAs
  resMat <- resMat[-1,]
  
  if(compress)
  {
    # id with duplicates
    idvars <- c(factor(resMat[,2])) 
    # correct ordering
    idvars <- (1:length(unique(idvars)))[factor(idvars)] 
    
    # rewrite index for actual matrix
    index <- unique(index)
    
  }else{
    # id with unique values
    idvars <- resMat[,3]
  }  
  
  new_time <- resMat[,1]
  
  newHmat <- hmatrix(time = new_time, 
                     id = idvars, 
                     x = attrTemp$x[index, , drop=FALSE], 
                     argvals = attrTemp$argvals, 
                     timeLab = attrTemp$timeLab, 
                     idLab = attrTemp$idLab, 
                     xLab = attrTemp$xLab, 
                     argvalsLab = attrTemp$argvalsLab)
  
  return(newHmat)
  
}
