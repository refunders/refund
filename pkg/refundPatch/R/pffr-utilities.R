# Utility functions for pffr:

safeDeparse <- function(expr){
	# turn an expression into a _single_ string, regardless of the expression's length 	
	ret <- paste(deparse(expr), collapse="")
	#rm whitespace
	gsub("[[:space:]][[:space:]]+", " ", ret)
}

#' Return call with all possible arguments
#' 
#' Return a call in which all of the arguments which were supplied or have presets are specified by their full names and their supplied or default values. 
#' 
#' @param definition a function. See \code{\link[base]{match.call}}.
#' @param call an unevaluated call to the function specified by definition. See \code{\link[base]{match.call}}.
#' @param expand.dots logical. Should arguments matching ... in the call be included or left as a ... argument? See \code{\link[base]{match.call}}.
#' @return An object of mode "\code{\link[base]{call}}". 
#' @author Fabian Scheipl
#' @export
#' @seealso \code{\link[base]{match.call}}
expand.call <-
        function(definition=NULL, call=sys.call(sys.parent(1)), expand.dots = TRUE)
{
    call <- match.call(definition, call, expand.dots)
    #given args:
    ans <- as.list(call)
    
    #possible args:
    frmls <- formals(safeDeparse(ans[[1]]))
    #remove formal args with no presets:
    frmls <- frmls[!sapply(frmls, is.symbol)]
    
    add <- which(!(names(frmls) %in% names(ans)))
    return(as.call(c(ans, frmls[add])))
}



list2df <- function(l){
# make a list into a dataframe -- matrices are left as matrices!
    nrows <- sapply(l, function(x) nrow(as.matrix(x)))
    stopifnot(length(unique(nrows)) == 1)
    ret <- data.frame(rep(NA, nrows[1]))
    for(i in 1:length(l)) ret[[i]] <- l[[i]]
    names(ret) <- names(l)
    return(ret)
}

## TODO: this does not always yield unique labels, e.g. if you have
##   s(g, bs="re") + s(g, bs="mrf", xt=somepenalty) 
getShrtlbls <- function(object){
# make short labels for display/coef-list, etc...
    
    labelmap <- object$pffr$labelmap
    
    ret <- sapply(names(unlist(labelmap)),
            function(x){
                #make a parseable expression for ffpc terms
                if(grepl("^ffpc", x)){
                    ffpcnumber <- gsub("(^.+))([0-9]+$)","\\2", x)
                    x <- gsub(")[0-9]+",")",x)
                   }
                exprx <- parse(text=x)
                
                #remove argument names
                x <- gsub("((?:[A-Za-z]+))(\\s)(=)(\\s)", "", x, perl=TRUE)
                
                #remove whitespace
                x <- gsub("\\s", "", x)
                
                #remove everything after and including first quoted argument
                if(any(regexpr(",\".*",x)[[1]]>0)) {
                    x <- gsub("([^\"]*)(,[c\\(]*\".*)(\\)$)", "\\1\\3", x, perl=TRUE) 
                }
                
                #remove everything after last variable:
                lstvrbl <- tail(all.vars(exprx),1)                
                x <- gsub(paste("(^.*?(?=",lstvrbl,")",lstvrbl,")(.*$)",sep=""), "\\1", x, perl=TRUE)
                
                #match braces
                openbr <- sum(grepl("\\(", strsplit(x, "")[[1]]))
                closebr <- sum(grepl("\\)", strsplit(x, "")[[1]]))
                if(openbr>closebr) x <- paste(c(x, rep(")", openbr-closebr)), sep="",collapse="") 
                
                #add number of PC for ffpc terms
                if(grepl("^ffpc", x)){
                    x <- paste(x, ffpcnumber, sep="")
                }
                return(x)
            })
    # correct labels for factor variables:
    if(any(sapply(labelmap, length) > 1 & !sapply(names(labelmap), function(x) grepl("^ffpc", x)))){
        which <- which(sapply(labelmap, length) > 1)
        inds <- c(0, cumsum(sapply(labelmap, length)))
        for(w in which){
            ret[(inds[w]+1):inds[w+1]] <- {
                lbls <- labelmap[[w]]
                bylevels <- sapply(object$smooth[lbls], function(x) x$by.level)
                by <- object$smooth[[lbls[1]]]$by
                paste(by, bylevels, "(", object$pffr$yindname, ")", sep="")
            }
        }
    }        
    #append labels for varying coefficient terms
    if(any(!grepl("\\(", ret))){
        which <- which(!grepl("\\(", ret))
        ret[which] <- paste(ret[which],"(", object$pffr$yindname, ")", sep="")
    }
    
    return(ret)
}
