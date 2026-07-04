## Clean WITHIN-PARTITION-FAMILY comparison (the paper's evidence foundation).
## EF (equal-quantile, chi2_{G-2}) vs HL, Pigeon-Heyse, Tsiatis, Xie. NO Stukel. Same G=10.
suppressWarnings(suppressMessages({ library(ResourceSelection); library(MASS); library(dplyr) }))
set.seed(20260703); lg<-function(p)log(p/(1-p)); PROJ<-"c:/Users/ebrah/.cursor-tutor/projects"
sw<-function(x) suppressWarnings(suppressMessages(x))
for(f in c("pigeonheyse.R","Tsiatis.R","Xie.R")) try(sw(source(file.path(PROJ,f))), silent=TRUE)
cat("available -> PH:",exists("pigeon_heyse_test")," Tsiatis:",exists("score_gof_clustering")," Xie:",exists("XieGoodnessOfFitTest"),"\n")
sc<-function(e){v<-tryCatch(sw(e),error=function(x)NA_real_); if(is.null(v)||length(v)!=1||!is.finite(as.numeric(v)))NA_real_ else as.numeric(v)}
sqb<-function(J)solve(rbind(c(1,-1.5,2.25),c(1,3,9),c(1,-3,9)),c(lg(.05),lg(.95),lg(J)))
sib<-function(I){e0<-lg(.1);e1<-lg(.2);e2<-lg(.2+I);b0<-(e0+e1)/2;b1<-(e1-e0)/6;b3<-(e2-(b0+3*b1))/6;c(b0,b1,3*b3,b3)}
ef_full<-function(y,ph,G=10){ph<-pmin(pmax(ph,1e-6),1-1e-6);n<-length(y)
  g<-pmin(ceiling(rank(ph,ties.method="first")/(n/G)),G)
  o<-tapply(y,g,sum);e<-tapply(ph,g,sum);ng<-tapply(y,g,length);pb<-as.numeric(tapply(ph,g,mean))
  V<-ng*pb*(1-pb);oe<-as.numeric(o-e);w<-1-2*pb; 1-pchisq(sum(oe^2/V)-sum(w*oe/V),G-2)}
gen<-function(scn,n){x<-runif(n,-3,3);d<-rbinom(n,1,.5);eta<-0.6*x+0.5*d
  if(scn=="null")   return(list(D=data.frame(x=x,d=d,y=rbinom(n,1,plogis(eta))),f=y~x+d))
  if(scn=="cloglog")return(list(D=data.frame(x=x,d=d,y=rbinom(n,1,1-exp(-exp(eta)))),f=y~x+d))
  if(scn=="probit") return(list(D=data.frame(x=x,d=d,y=rbinom(n,1,pnorm(eta))),f=y~x+d))
  if(scn=="quad"){b<-sqb(0.03);return(list(D=data.frame(x=x,y=rbinom(n,1,plogis(b[1]+b[2]*x+b[3]*x^2))),f=y~x))}
  if(scn=="int_bin"){b<-sib(0.5);return(list(D=data.frame(x=x,d=d,y=rbinom(n,1,plogis(b[1]+b[2]*x+b[3]*d+b[4]*x*d))),f=y~x+d))}}
one<-function(scn,n){G0<-gen(scn,n);dat<-G0$D;fit<-sw(glm(G0$f,data=dat,family=binomial()));ph<-as.numeric(fitted(fit));y<-dat$y
  f2<-fit;f2$predicted_probs<-ph
  c(EF=sc(ef_full(y,ph,10)), HL=sc(hoslem.test(y,ph,g=10)$p.value),
    PH=sc(pigeon_heyse_test(data.frame(y=y),fit,g=10)$p_value),
    Tsiatis=sc(score_gof_clustering(fit,num_groups=10,y=y)$p_value),
    Xie=sc(as.numeric(XieGoodnessOfFitTest(dat,f2))))}
cat("\nDEBUG one rep (cloglog,1000):\n"); print(one("cloglog",1000))
runf<-function(scn,n,reps){ tests<-c("EF","HL","PH","Tsiatis","Xie")
  m<-sapply(seq_len(reps), function(i) one(scn,n)); m<-matrix(as.numeric(m), nrow=length(tests), dimnames=list(tests,NULL))
  rej<-round(rowMeans(m<0.05, na.rm=TRUE),3); na<-round(rowMeans(is.na(m)),2)
  cat(sprintf("  [%s n=%d] NA-rates: %s\n", scn, n, paste(sprintf("%s=%.2f",tests,na),collapse=" ")))
  data.frame(scenario=scn,n=n, as.list(rej)) }
R<-400
cat("\n============ TYPE I (should be ~0.05; PH known conservative) ============\n")
tab <- rbind(runf("null",500,R), runf("null",1000,R))
print(tab, row.names=FALSE)
cat("\n============ POWER (n=1000) ============\n")
pw <- rbind(runf("cloglog",1000,R), runf("probit",1000,R), runf("quad",1000,R), runf("int_bin",1000,R))
print(pw, row.names=FALSE)
write.csv(rbind(tab,pw), "within_family_results.csv", row.names=FALSE)
cat("\nHEADLINE CHECK: on cloglog, is EF the max of {EF,HL,PH,Tsiatis,Xie}?\n")
cl <- pw[pw$scenario=="cloglog",c("EF","HL","PH","Tsiatis","Xie")]
cat("  cloglog powers:", paste(names(cl), unlist(cl), sep="="), "\n")
cat("  EF is best-in-family on cloglog:", which.max(unlist(cl))==1, "\n")
