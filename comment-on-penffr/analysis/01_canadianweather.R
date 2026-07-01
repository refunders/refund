## 01_canadianweather.R — re-run pffr on Canadian Weather (concurrent + integral).
##  * "as-is"   reproduces the authors' own call (algorithm="bam", cr bases,
##              lat/lon entered as n x 365 matrices => time-varying k=50 effects).
##  * "matched" is specified as close to PenFFR as possible: cubic B-spline (ps)
##              bases with a 2nd-derivative penalty (m=c(2,2)); coefficient bases
##              matched to PenFFR (40 concurrent, 10x10 integral); scalars entered
##              as CONSTANT effects c(lat)+c(lon) (as PenFFR does); historical
##              integration limits="s<=t"; and algorithm="gam" with exact REML
##              (bam's fREML approximation tends to give slightly worse fits).
## Metric + protocol are the authors' own (leave-one-out ISE over 35 stations).
source(file.path("comment-on-penffr","analysis","00_setup.R"))
load_pffr(); suppressMessages(library(MASS))  # MASS::Null needed by ff(limits=...)
CW <- get_canadian_weather()
raw.temp <- t(CW[["dailyAv"]][,,1]); raw.prec <- t(CW[["dailyAv"]][,,3])  # [,,3] = log10 precip
coord <- as.matrix(CW[["coordinates"]])
n <- 35; m <- 365; obs.grid <- ((1:m)-1)/(m-1)
run_loo <- function(f,label){ t0<-Sys.time()
  e <- sapply(1:n, function(i){ x<-try(suppressWarnings(f(i)),silent=TRUE); if(inherits(x,"try-error")) NA_real_ else x })
  cat(sprintf("%-34s ISE mean=%.2f (sd %.2f) | %d/%d | %.0fs\n", label, mean(e,na.rm=TRUE),
      sd(e,na.rm=TRUE), sum(!is.na(e)), n, as.numeric(Sys.time()-t0,units="secs"))); e }

conc_asis <- function(i){                              # authors' call (bam, lat/lon as matrices)
  tmp.prec<-raw.prec[-i,]; tmp.temp<-raw.temp[-i,]
  lat<-matrix(rep(coord[-i,1],m),ncol=m); lon<-matrix(rep(coord[-i,2],m),ncol=m)
  mod<-pffr(tmp.prec~tmp.temp+lat+lon, yind=obs.grid, algorithm="bam",
            bs.int=list(bs="cr",k=50), bs.yindex=list(bs="cr",k=50), sandwich="none")
  latP<-matrix(rep(coord[i,1],m),nrow=1); lonP<-matrix(rep(coord[i,2],m),nrow=1)
  ISE_curve(raw.prec[i,], predict(mod,type="response",
            newdata=list(tmp.temp=matrix(raw.temp[i,],nrow=1),lat=latP,lon=lonP))) }
conc_matched <- function(i){                           # matched to PenFFR (gam/REML, c() scalars)
  tmp.prec<-raw.prec[-i,]; tmp.temp<-raw.temp[-i,]; lat<-coord[-i,1]; lon<-coord[-i,2]
  mod<-pffr(tmp.prec~tmp.temp+c(lat)+c(lon), yind=obs.grid, algorithm="gam", method="REML",
            bs.int=list(bs="ps",m=c(2,2),k=40), bs.yindex=list(bs="ps",m=c(2,2),k=40), sandwich="none")
  ISE_curve(raw.prec[i,], predict(mod,type="response",
            newdata=list(tmp.temp=matrix(raw.temp[i,],nrow=1),lat=coord[i,1],lon=coord[i,2]))) }
inte_asis <- function(i){                              # authors' call (bam, full-range ff, lat/lon matrices)
  tmp.prec<-raw.prec[-i,]; tmp.temp<-raw.temp[-i,]
  lat<-matrix(rep(coord[-i,1],m),ncol=m); lon<-matrix(rep(coord[-i,2],m),ncol=m)
  mod<-pffr(tmp.prec~ff(tmp.temp,xind=obs.grid,check.ident=FALSE,basistype="te",
            integration="rieman",splinepars=list(bs="cr",k=10))+lat+lon, yind=obs.grid,
            algorithm="bam", bs.int=list(bs="cr",k=50), bs.yindex=list(bs="cr",k=50), sandwich="none")
  latP<-matrix(rep(coord[i,1],m),nrow=1); lonP<-matrix(rep(coord[i,2],m),nrow=1)
  ISE_curve(raw.prec[i,], predict(mod,type="response",
            newdata=list(tmp.temp=matrix(raw.temp[i,],nrow=1),lat=latP,lon=lonP))) }
inte_matched <- function(i){                           # matched: historical, gam/REML, c() scalars
  tmp.prec<-raw.prec[-i,]; tmp.temp<-raw.temp[-i,]; lat<-coord[-i,1]; lon<-coord[-i,2]
  mod<-pffr(tmp.prec~ff(tmp.temp,xind=obs.grid,limits="s<=t",
            splinepars=list(bs="ps",m=list(c(2,2),c(2,2)),k=c(10,10)))+c(lat)+c(lon), yind=obs.grid,
            algorithm="gam", method="REML", bs.int=list(bs="ps",m=c(2,2),k=10),
            bs.yindex=list(bs="ps",m=c(2,2),k=10), sandwich="none")
  ISE_curve(raw.prec[i,], predict(mod,type="response",
            newdata=list(tmp.temp=matrix(raw.temp[i,],nrow=1),lat=coord[i,1],lon=coord[i,2]))) }

res <- list(
  conc_asis    = run_loo(conc_asis,    "Concurrent pffr (authors as-is, bam):"),
  conc_matched = run_loo(conc_matched, "Concurrent pffr (matched, gam/REML):"),
  inte_asis    = run_loo(inte_asis,    "Integral pffr (authors as-is, bam):"),
  inte_matched = run_loo(inte_matched, "Integral pffr (matched, gam/REML, hist):"))
out <- data.frame(model=names(res), ISE_mean=sapply(res,mean,na.rm=TRUE), ISE_sd=sapply(res,sd,na.rm=TRUE))
write.csv(out, file.path("comment-on-penffr","results","tables","canadianweather_pffr.csv"), row.names=FALSE)
cat("\nPaper Table 4: Concurrent pffr 89.31 (52.03); Integral pffr 41.37 (48.91);",
    "Concurrent PenFFR 36.40; Integral PenFFR 33.66 (best).\n")
