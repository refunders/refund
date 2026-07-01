## 02_hawaii_ocean.R — re-run pffr on the Hawaii Ocean data (concurrent + integral).
## Predict Salinity from Potential.density, Temperature, Oxygen, Chloropigment.
## Train = first 50 curves, test = remaining 66 (as in the paper). Metric: ISE x 100.
source(file.path("comment-on-penffr","analysis","00_setup.R"))
load_pffr()
ocean <- get_ocean()
Sal<-ocean$Salinity; Dens<-ocean$Potential.density; Temp<-ocean$Temperature
Oxy<-ocean$Oxygen; Chl<-ocean$Chloropigment
m<-101; depth<-seq(0,1,length.out=m); tr<-1:50; te<-51:116
ISE_set <- function(Yact,Yhat) rowSums((Yact-Yhat)^2)
dens<-Dens[tr,]; temp<-Temp[tr,]; oxy<-Oxy[tr,]; chl<-Chl[tr,]; sal<-Sal[tr,]

## concurrent pffr (paper Table 5 reports this as the BEST method: 0.52 x10^-2)
mc <- pffr(sal~dens+temp+oxy+chl, yind=depth, algorithm="gam", method="REML",
           bs.yindex=list(bs="ps",k=20), bs.int=list(bs="ps",k=20), sandwich="none")
pc <- predict(mc, type="response", newdata=list(dens=Dens[te,],temp=Temp[te,],oxy=Oxy[te,],chl=Chl[te,]))
ec <- ISE_set(Sal[te,], pc)

## integral pffr
mi <- pffr(sal~ff(dens,xind=depth,check.ident=FALSE)+ff(temp,xind=depth,check.ident=FALSE)+
                ff(oxy,xind=depth,check.ident=FALSE)+ff(chl,xind=depth,check.ident=FALSE),
           yind=depth, algorithm="gam", method="REML", bs.yindex=list(bs="ps",k=20), bs.int=list(bs="ps",k=20), sandwich="none")
pii <- predict(mi, type="response", newdata=list(dens=Dens[te,],temp=Temp[te,],oxy=Oxy[te,],chl=Chl[te,]))
ei <- ISE_set(Sal[te,], pii)

out <- data.frame(model=c("Concurrent pffr","Integral pffr"),
                  ISEx100_mean=c(100*mean(ec),100*mean(ei)), ISEx100_sd=c(100*sd(ec),100*sd(ei)))
print(out)
write.csv(out, file.path("comment-on-penffr","results","tables","hawaii_pffr.csv"), row.names=FALSE)
cat("\nPaper Table 5 (x10^2): Integral PenFFR 0.57; Concurrent PenFFR 1.83; Integral pffr 2.37;",
    "Concurrent pffr 0.52 (BEST, the authors' own bold); wSigcomp 4.79.\n")
