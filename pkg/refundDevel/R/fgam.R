fgam <- function(formula,fitter=NA,tensortype=c('te','t2'),...){
  
  call <- match.call()
  tensortype <- as.symbol(match.arg(tensortype))
  dots <- list(...)
  if (length(dots)) {
    validDots <- if (!is.na(fitter) && fitter == "gamm4") {
      c(names(formals(gamm4)), names(formals(lmer)))
    }
    else {
      c(names(formals(gam)), names(formals(gam.fit)))
    }
    notUsed <- names(dots)[!(names(dots) %in% validDots)]
    if (length(notUsed)) 
      warning("Arguments <", paste(notUsed, collapse = ", "), 
              "> supplied but not used.")
  }
  tf <- terms.formula(formula, specials = c("s", "te", "t2", "lf", "af"))
  trmstrings <- attr(tf, "term.labels")
  terms <- sapply(trmstrings, function(trm) as.call(parse(text = trm))[[1]], 
                  simplify = FALSE)
  frmlenv <- environment(formula)
  where.af <- attr(tf, "specials")$af - 1
  where.lf <- attr(tf, "specials")$lf - 1
  where.s <- attr(tf, "specials")$s - 1
  where.te <- attr(tf, "specials")$te - 1
  where.t2 <- attr(tf, "specials")$t2 - 1
  
  
  if (length(trmstrings)) {
    where.par <- which(!(1:length(trmstrings) %in% c(where.af, where.lf, where.s,where.te, where.t2)))
  }else where.par <- numeric(0)
  
  responsename <- attr(tf, "variables")[2][[1]]
  newfrml <- paste(responsename, "~", sep = "")
  newfrmlenv <- new.env()
  evalenv <- if ("data" %in% names(call)) 
    eval(call$data)
  else NULL
  nobs <- length(eval(responsename, envir = evalenv, enclos = frmlenv))
  
  if (missing(fitter) || is.na(fitter)) {
    fitter <- ifelse(nobs > 1e+05, "bam", "gam")
  }

  fitter <- as.symbol(fitter)
  if (as.character(fitter) == "bam" && !("chunk.size" %in% 
    names(call))) {
    call$chunk.size <- max(nobs/5, 10000)
  }
  if (as.character(fitter) == "gamm4") 
    stopifnot(length(where.te) < 1)
  
  assign(x = deparse(responsename), value = as.vector(t(eval(responsename, 
                                      envir = evalenv, enclos = frmlenv))), envir = newfrmlenv)

  newtrmstrings <- attr(tf, "term.labels")
  if (!attr(tf, "intercept")) {
    newfrml <- paste(newfrml, "0", sep = "")
  }
  
  if (length(c(where.af, where.lf))) {
    fterms <- lapply(terms[c(where.af, where.lf)], function(x) {
      eval(x, envir = evalenv, enclos = frmlenv)
    })
    newtrmstrings[c(where.af, where.lf)] <- sapply(fterms, 
                                                    function(x) {
                                                      safeDeparse(x$call)
                                                    })
    lapply(fterms, function(x) {
      lapply(names(x$data), function(nm) {
        assign(x = nm, value = x$data[[nm]], envir = newfrmlenv)
        invisible(NULL)
      })
      invisible(NULL)
    })
    fterms <- lapply(fterms, function(x) x[names(x) != 
      "data"])
  }
  else fterms <- NULL
  
  where.notf <- c(where.par,where.s,where.te,where.t2)
  if (length(where.notf)) {
    if ("data" %in% names(call)) 
      frmlenv <- list2env(eval(call$data), frmlenv)
    lapply(terms[where.notf], function(x) {

      nms <- if (!is.null(names(x))) {
        all.vars(x[names(x) == ""])
      }
      else all.vars(x)
      sapply(nms, function(nm) {
        stopifnot(length(get(nm, envir = frmlenv)) == 
          nobs)
        assign(x = nm, value = get(nm, envir = frmlenv), envir = newfrmlenv)
        invisible(NULL)
      })
      invisible(NULL)
    })
  }
  
  newfrml <- formula(paste(c(newfrml, newtrmstrings), collapse = "+"))
  environment(newfrml) <- newfrmlenv
  fgamdata <- list2df(as.list(newfrmlenv))
  datameans <- sapply(as.list(newfrmlenv),mean)
  newcall <- expand.call(fgam, call)
  newcall$fitter <- type <- newcall$bs.int <- newcall$bs.yindex <- newcall$fitter <- NULL
  newcall$formula <- newfrml
  newcall$data <- quote(fgamdata)
  newcall$fitter <- newcall$tensortype <- NULL
  newcall[[1]] <- fitter
  
  
  res <- eval(newcall)
  
  res.smooth <- if (as.character(fitter) %in% c("gamm4", "gamm")) {
    res$gam$smooth
  }else res$smooth
  
  trmmap <- newtrmstrings
  names(trmmap) <- names(terms)
  labelmap <- as.list(trmmap)
  names(trmmap) <- names(terms)

    
  lbls <- sapply(res.smooth, function(x) x$label)
  if (length(where.par)) {
    for (w in where.par) labelmap[[w]] <- {
      where <- sapply(res.smooth, function(x) x$by) == 
        names(labelmap)[w]
      sapply(res.smooth[where], function(x) x$label)
    }
  labelmap[-c(where.par)] <- lbls[pmatch(sapply(labelmap[-c(where.par)], function(x) {
                                                                tmp <- eval(parse(text = x))
                                                                if (is.list(tmp)) {
                                                                  return(tmp$label)
                                                                }else {
                                                                  return(x)
                                                                }
                                                              }), lbls)]
  }else{
    labelmap[1:length(labelmap)] <- lbls[pmatch(sapply(labelmap, 
                                                       function(x) {
                                                         tmp <- eval(parse(text = x))
                                                         if (is.list(tmp)) {
                                                           return(tmp$label)
                                                         }
                                                         else {
                                                           return(x)
                                                         }
                                                       }), lbls)]
  }
  if (any(nalbls <- sapply(labelmap, function(x) any(is.na(x))))) {
    labelmap[nalbls] <- trmmap[nalbls]
  }
  names(res.smooth) <- lbls
  if (as.character(fitter) %in% c("gamm4", "gamm")) {
    res$gam$smooth <- res.smooth
  }else {
    res$smooth <- res.smooth
  }
  ret <- list(formula = formula, termmap = trmmap, labelmap = labelmap, 
              responsename = responsename, nobs = nobs,
              where = list(where.af = where.af, where.lf = where.lf, 
                              where.s = where.s, where.te = where.te, where.t2 = where.t2, 
                               where.par = where.par), datameans=datameans,ft = fterms)
  if (as.character(fitter) %in% c("gamm4", "gamm")) {
    res$gam$fgam <- ret
    class(res$gam) <- c("fgam", class(res$gam))
  }
  else {
    res$fgam <- ret
    class(res) <- c("fgam", class(res))
  }

	return(res)
}
