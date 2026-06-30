## 01_canadianweather.R — re-run pffr on Canadian Weather (concurrent + integral),
## as the authors call it vs. with the lat/lon mis-specification corrected.
## Reproduces the authors' LOO ISE protocol and metric.
source(file.path("comment-on-penffr","analysis","00_setup.R"))  # adjust path if run from analysis/
load_pffr()
CW <- get_canadian_weather()
raw.temp <- t(CW[["dailyAv"]][,,1]); raw.prec <- t(CW[["dailyAv"]][,,3])  # [,,3] = log10 precip
coord <- as.matrix(CW[["coordinates"]])
n <- 35; m <- 365; obs.grid <- ((1:m)-1)/(m-1)

run_loo <- function(fitfun, label) {
  t0 <- Sys.time()
  e <- sapply(1:n, function(i){ x <- try(fitfun(i), silent=TRUE)
    if (inherits(x,"try-error")) NA_real_ else x })
  cat(sprintf("%-34s ISE mean=%.2f (sd %.2f) | %d/%d ok | %.0fs\n",
      label, mean(e,na.rm=TRUE), sd(e,na.rm=TRUE), sum(!is.na(e)), n,
      as.numeric(Sys.time()-t0, units="secs"))); e
}

## ---------- concurrent ----------
conc_asis <- function(i){
  tmp.prec<-raw.prec[-i,]; tmp.temp<-raw.temp[-i,]
  lat<-matrix(rep(coord[-i,1],m),ncol=m); lon<-matrix(rep(coord[-i,2],m),ncol=m)
  mod<-pffr(tmp.prec~tmp.temp+lat+lon, yind=obs.grid, algorithm="bam",
            bs.int=list(bs="cr",k=50), bs.yindex=list(bs="cr",k=50), sandwich="none")
  latP<-matrix(rep(coord[i,1],m),nrow=1); lonP<-matrix(rep(coord[i,2],m),nrow=1)
  ISE_curve(raw.prec[i,], predict(mod,type="response",
            newdata=list(tmp.temp=matrix(raw.temp[i,],nrow=1),lat=latP,lon=lonP)))
}
conc_fair <- function(i){
  tmp.prec<-raw.prec[-i,]; tmp.temp<-raw.temp[-i,]; lat<-coord[-i,1]; lon<-coord[-i,2]
  mod<-pffr(tmp.prec~tmp.temp+c(lat)+c(lon), yind=obs.grid, algorithm="bam",
            bs.int=list(bs="ps",k=40), bs.yindex=list(bs="ps",k=40), sandwich="none")
  ISE_curve(raw.prec[i,], predict(mod,type="response",
            newdata=list(tmp.temp=matrix(raw.temp[i,],nrow=1),lat=coord[i,1],lon=coord[i,2])))
}
## ---------- integral ----------
inte_asis <- function(i){
  tmp.prec<-raw.prec[-i,]; tmp.temp<-raw.temp[-i,]
  lat<-matrix(rep(coord[-i,1],m),ncol=m); lon<-matrix(rep(coord[-i,2],m),ncol=m)
  mod<-pffr(tmp.prec~ff(tmp.temp,xind=obs.grid,check.ident=FALSE,basistype="te",
            integration="rieman",splinepars=list(bs="cr",k=10))+lat+lon, yind=obs.grid,
            algorithm="bam", bs.int=list(bs="cr",k=50), bs.yindex=list(bs="cr",k=50), sandwich="none")
  latP<-matrix(rep(coord[i,1],m),nrow=1); lonP<-matrix(rep(coord[i,2],m),nrow=1)
  ISE_curve(raw.prec[i,], predict(mod,type="response",
            newdata=list(tmp.temp=matrix(raw.temp[i,],nrow=1),lat=latP,lon=lonP)))
}
inte_fair <- function(i){
  tmp.prec<-raw.prec[-i,]; tmp.temp<-raw.temp[-i,]; lat<-coord[-i,1]; lon<-coord[-i,2]
  mod<-pffr(tmp.prec~ff(tmp.temp,xind=obs.grid,check.ident=FALSE,basistype="te",
            integration="rieman",splinepars=list(bs="ps",k=10))+c(lat)+c(lon), yind=obs.grid,
            algorithm="bam", bs.int=list(bs="ps",k=40), bs.yindex=list(bs="ps",k=40), sandwich="none")
  ISE_curve(raw.prec[i,], predict(mod,type="response",
            newdata=list(tmp.temp=matrix(raw.temp[i,],nrow=1),lat=coord[i,1],lon=coord[i,2])))
}

res <- list(
  conc_asis = run_loo(conc_asis, "Concurrent pffr (authors as-is):"),
  conc_fair = run_loo(conc_fair, "Concurrent pffr (lat/lon fixed):"),
  inte_asis = run_loo(inte_asis, "Integral pffr (authors as-is):"),
  inte_fair = run_loo(inte_fair, "Integral pffr (lat/lon fixed):"))
out <- data.frame(model=names(res),
                  ISE_mean=sapply(res,mean,na.rm=TRUE), ISE_sd=sapply(res,sd,na.rm=TRUE))
write.csv(out, file.path("comment-on-penffr","results","tables","canadianweather_pffr.csv"), row.names=FALSE)
cat("\nPaper Table 4: Concurrent pffr 89.31 (52.03); Integral pffr 41.37 (48.91);",
    "Concurrent PenFFR 36.40; Integral PenFFR 33.66 (best).\n")
