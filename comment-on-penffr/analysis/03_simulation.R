## 03_simulation.R — basis-count (in)sensitivity of penalized pffr.
## Faithful reconstruction of the paper's Sect 5.1 CONCURRENT design (p=5 functional
## covariates, smooth time-varying coefficients, Gaussian noise) with the paper's MRPE.
## Shows: (a) penalized pffr is flat in the basis dim k; (b) its error tracks the noise.
source(file.path("comment-on-penffr","analysis","00_setup.R"))
load_pffr()
suppressMessages(library(splines))

m <- 80; tg <- seq(0,1,length.out=m)
B <- cbind(2.5*sin(2*pi*tg), cos(2*pi*tg), -1.5-1.5*cos(pi*tg), 0.5*sin(2*pi*tg), 0.3*cos(2*pi*tg))
b0 <- 0.7 + 0.2*tg
simX <- function(n){ P<-6; ph<-sapply(1:P,function(r) sin(2*pi*r*tg)); ph2<-sapply(1:P,function(r) cos(2*pi*r*tg))
  lapply(1:5, function(l){ a<-matrix(rnorm(n*P,sd=1/(1:P)),n,P,byrow=TRUE); a2<-matrix(rnorm(n*P,sd=1/(1:P)),n,P,byrow=TRUE)
    (a%*%t(ph)+a2%*%t(ph2)) + l*0.2 }) }
genY <- function(Xl,s2){ n<-nrow(Xl[[1]]); Y<-matrix(rep(b0,n),n,m,byrow=TRUE)
  for(l in 1:5) Y<-Y+Xl[[l]]*matrix(rep(B[,l],n),n,m,byrow=TRUE); Y+matrix(rnorm(n*m,sd=sqrt(s2)),n,m) }
MRPE <- function(Ya,Yh) mean(rowSums((Ya-Yh)^2)/rowSums(Ya^2))
fit_pffr <- function(Xtr,Ytr,Xte,k){ X1<-Xtr[[1]];X2<-Xtr[[2]];X3<-Xtr[[3]];X4<-Xtr[[4]];X5<-Xtr[[5]]
  mod<-pffr(Ytr~X1+X2+X3+X4+X5, yind=tg, algorithm="bam",
            bs.yindex=list(bs="ps",k=k), bs.int=list(bs="ps",k=k), sandwich="none")
  predict(mod,type="response",newdata=list(X1=Xte[[1]],X2=Xte[[2]],X3=Xte[[3]],X4=Xte[[4]],X5=Xte[[5]])) }

set.seed(2026); scen<-list(c(200,1),c(200,4),c(500,1),c(500,4)); ks<-c(5,10,20,40); reps<-5
rows<-list()
for(s in scen){ n<-s[1]; s2<-s[2]; for(k in ks){ mp<-numeric(reps)
  for(r in 1:reps){ Xtr<-simX(n);Ytr<-genY(Xtr,s2);Xte<-simX(200);Yte<-genY(Xte,s2); mp[r]<-MRPE(Yte,fit_pffr(Xtr,Ytr,Xte,k)) }
  rows[[length(rows)+1]]<-data.frame(n=n,sigma2=s2,k=k,pffr=mean(mp),pffr_sd=sd(mp))
  cat(sprintf("n=%d s2=%d k=%2d | pffr MRPE=%.4f (sd %.4f)\n",n,s2,k,mean(mp),sd(mp))) }}
out<-do.call(rbind,rows)
write.csv(out, file.path("comment-on-penffr","results","tables","sim_results.csv"), row.names=FALSE)

png(file.path("comment-on-penffr","results","figures","sim_basis_insensitivity.png"),width=1500,height=650,res=150)
par(mfrow=c(1,2),mar=c(4,4,3,1)); cols<-c("#1b9e77","#d95f02","#1b9e77","#d95f02"); lt<-c(1,1,2,2)
plot(NA,xlim=range(out$k),ylim=range(out$pffr),xlab="basis dim k",ylab="test MRPE",main="(a) penalized pffr is flat in k")
j<-0; for(s in scen){j<-j+1; sub<-out[out$n==s[1]&out$sigma2==s[2],]; lines(sub$k,sub$pffr,col=cols[j],lty=lt[j],lwd=2); points(sub$k,sub$pffr,col=cols[j],pch=19)}
legend("right",c("n=200,s2=1","n=200,s2=4","n=500,s2=1","n=500,s2=4"),col=cols,lty=lt,lwd=2,bty="n",cex=.8)
sub<-out[out$k==20,]; bp<-barplot(sub$pffr,names.arg=paste0("n",sub$n,"\ns2=",sub$sigma2),col=ifelse(sub$sigma2==1,"#1b9e77","#d95f02"),ylab="test MRPE (k=20)",main="(b) pffr error tracks noise",ylim=c(0,.35)); text(bp,sub$pffr+.012,sprintf("%.3f",sub$pffr),cex=.85)
dev.off()
cat("\nPaper Table 2: pffr MRPE x10^3 ~ 91.5 (constant across all 4 scenarios); FFR/PenFFR ~54.\n")
