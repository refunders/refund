## 04_robustness.R — robustness of the corrected Canadian-Weather pffr result.
## Shows that, once the two scalar coordinates are entered correctly, the matched
## pffr ISE is essentially invariant to the remaining modelling choices; the only
## lever is the scalar-covariate treatment itself (see 01_canadianweather.R).
##
## Variants reproduced here (LOO ISE over 35 stations; cf. Section 4.3 of the Comment):
##   basis family      : cr vs cubic B-spline (ps) + 2nd-derivative penalty
##   smoothing criterion: exact REML vs GCV.Cp
##   fitting engine     : gam (exact REML) vs bam (fast REML)  [see 01_* for bam as-is]
##   covariate pre-smoothing: raw temperature vs a 100-basis B-spline projection
##                            (exactly PenFFR's my_functrep step)
## NOTE: several integral fits are involved; a full run takes ~1-2 hours.
source(file.path("comment-on-penffr","analysis","00_setup.R"))
load_pffr(); suppressMessages(library(MASS))     # MASS::Null needed by ff(limits=...)
CW <- get_canadian_weather()
raw.temp <- t(CW[["dailyAv"]][,,1]); raw.prec <- t(CW[["dailyAv"]][,,3])
coord <- as.matrix(CW[["coordinates"]])
n <- 35; m <- 365; obs.grid <- ((1:m)-1)/(m-1)

## PenFFR-style pre-smoothing: project each temperature curve onto 100 cubic B-splines
## by unpenalized least squares (the no-argvals branch of my_functrep), computed once.
Bpre <- splines::bs(obs.grid, df = 100, intercept = TRUE)
temp.sm <- t(apply(raw.temp, 1, function(y) as.vector(Bpre %*% coef(lm.fit(Bpre, y)))))
cat(sprintf("pre-smoothing reconstruction: max|temp.sm-raw.temp| = %.3g (temp range %.1f)\n",
            max(abs(temp.sm - raw.temp)), diff(range(raw.temp))))

run_loo <- function(f,label){ t0<-Sys.time()
  e <- sapply(1:n, function(i){ x<-try(suppressWarnings(f(i)),silent=TRUE); if(inherits(x,"try-error")) NA_real_ else x })
  cat(sprintf("%-52s ISE mean=%.2f (sd %.2f) | %d/%d | %.0fs\n", label, mean(e,na.rm=TRUE),
      sd(e,na.rm=TRUE), sum(!is.na(e)), n, as.numeric(Sys.time()-t0,units="secs"))); e }

## --- concurrent: basis family (authors' as-is settings, cr -> cubic B-spline) ---
conc_asis_bspline <- function(i){
  tmp.prec<-raw.prec[-i,]; tmp.temp<-raw.temp[-i,]
  lat<-matrix(rep(coord[-i,1],m),ncol=m); lon<-matrix(rep(coord[-i,2],m),ncol=m)
  mod<-pffr(tmp.prec~tmp.temp+lat+lon, yind=obs.grid, algorithm="bam",
            bs.int=list(bs="ps",m=c(2,2),k=50), bs.yindex=list(bs="ps",m=c(2,2),k=50), sandwich="none")
  latP<-matrix(rep(coord[i,1],m),nrow=1); lonP<-matrix(rep(coord[i,2],m),nrow=1)
  ISE_curve(raw.prec[i,], predict(mod,type="response",
            newdata=list(tmp.temp=matrix(raw.temp[i,],nrow=1),lat=latP,lon=lonP))) }

## --- concurrent: matched + GCV instead of REML ---
conc_matched_gcv <- function(i){
  tmp.prec<-raw.prec[-i,]; tmp.temp<-raw.temp[-i,]; lat<-coord[-i,1]; lon<-coord[-i,2]
  mod<-pffr(tmp.prec~tmp.temp+c(lat)+c(lon), yind=obs.grid, algorithm="gam", method="GCV.Cp",
            bs.int=list(bs="ps",m=c(2,2),k=40), bs.yindex=list(bs="ps",m=c(2,2),k=40), sandwich="none")
  ISE_curve(raw.prec[i,], predict(mod,type="response",
            newdata=list(tmp.temp=matrix(raw.temp[i,],nrow=1),lat=coord[i,1],lon=coord[i,2]))) }

## --- concurrent: matched + PenFFR-style pre-smoothed covariate ---
conc_matched_presm <- function(i){
  tmp.prec<-raw.prec[-i,]; tmp.temp<-temp.sm[-i,]; lat<-coord[-i,1]; lon<-coord[-i,2]
  mod<-pffr(tmp.prec~tmp.temp+c(lat)+c(lon), yind=obs.grid, algorithm="gam", method="REML",
            bs.int=list(bs="ps",m=c(2,2),k=40), bs.yindex=list(bs="ps",m=c(2,2),k=40), sandwich="none")
  ISE_curve(raw.prec[i,], predict(mod,type="response",
            newdata=list(tmp.temp=matrix(temp.sm[i,],nrow=1),lat=coord[i,1],lon=coord[i,2]))) }

## --- integral (historical): matched, basis family (cubic B-spline -> cr) ---
inte_matched_cr <- function(i){
  tmp.prec<-raw.prec[-i,]; tmp.temp<-raw.temp[-i,]; lat<-coord[-i,1]; lon<-coord[-i,2]
  mod<-pffr(tmp.prec~ff(tmp.temp,xind=obs.grid,limits="s<=t",
            splinepars=list(bs="cr",k=c(10,10)))+c(lat)+c(lon), yind=obs.grid,
            algorithm="gam", method="REML", bs.int=list(bs="cr",k=10),
            bs.yindex=list(bs="cr",k=10), sandwich="none")
  ISE_curve(raw.prec[i,], predict(mod,type="response",
            newdata=list(tmp.temp=matrix(raw.temp[i,],nrow=1),lat=coord[i,1],lon=coord[i,2]))) }

## --- integral (historical): matched + GCV instead of REML ---
inte_matched_gcv <- function(i){
  tmp.prec<-raw.prec[-i,]; tmp.temp<-raw.temp[-i,]; lat<-coord[-i,1]; lon<-coord[-i,2]
  mod<-pffr(tmp.prec~ff(tmp.temp,xind=obs.grid,limits="s<=t",
            splinepars=list(bs="ps",m=list(c(2,2),c(2,2)),k=c(10,10)))+c(lat)+c(lon), yind=obs.grid,
            algorithm="gam", method="GCV.Cp", bs.int=list(bs="ps",m=c(2,2),k=10),
            bs.yindex=list(bs="ps",m=c(2,2),k=10), sandwich="none")
  ISE_curve(raw.prec[i,], predict(mod,type="response",
            newdata=list(tmp.temp=matrix(raw.temp[i,],nrow=1),lat=coord[i,1],lon=coord[i,2]))) }

res <- list(
  conc_asis_bspline  = run_loo(conc_asis_bspline,  "Concurrent, as-is, cubic B-spline + 2nd-deriv:"),
  conc_matched_gcv   = run_loo(conc_matched_gcv,   "Concurrent, matched, GCV.Cp:"),
  conc_matched_presm = run_loo(conc_matched_presm, "Concurrent, matched, pre-smoothed covariate:"),
  inte_matched_cr    = run_loo(inte_matched_cr,    "Integral (hist), matched, cr basis:"),
  inte_matched_gcv   = run_loo(inte_matched_gcv,   "Integral (hist), matched, GCV.Cp:"))
out <- data.frame(model=names(res), ISE_mean=sapply(res,mean,na.rm=TRUE), ISE_sd=sapply(res,sd,na.rm=TRUE))
write.csv(out, file.path("comment-on-penffr","results","tables","canadianweather_robustness.csv"), row.names=FALSE)

cat("\nReference (matched, from 01_canadianweather.R): concurrent 64.80 / integral 37.98.\n",
    "Expected here: as-is B-spline ~89.9; matched GCV ~64.9; pre-smoothed ~65.0;",
    "integral cr ~39.6; integral GCV ~38.7. Only the scalar-covariate fix (89.5->64.8) moves the result.\n")
