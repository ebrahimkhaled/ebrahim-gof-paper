## DEFINITIVE simulation for the EF paper -- PARALLEL (24-core), K=5000.
## EF vs HL. Size (g x covariate x n) + quadratic/interaction sweeps + the money LINK-ASYMMETRY
## sweep via the Aranda-Ordaz family (alpha=1 -> logit/symmetric ; alpha->0 -> cloglog/asymmetric).
## Also A(delta) per asymmetry level, for the mechanism figure.
suppressMessages(library(parallel)); lg<-function(p)log(p/(1-p))
## ---- tests (self-contained; no external scripts -> no segfault) ----
grp<-function(ph,G){ph<-pmin(pmax(ph,1e-6),1-1e-6);n<-length(ph);pmin(ceiling(rank(ph,ties.method="first")/(n/G)),G)}
ef_chi<-function(y,ph,G){g<-grp(ph,G);o<-tapply(y,g,sum);e<-tapply(ph,g,sum);ng<-tapply(y,g,length);pb<-as.numeric(tapply(ph,g,mean))
  V<-ng*pb*(1-pb);oe<-as.numeric(o-e);1-pchisq(sum(oe^2/V)-sum((1-2*pb)*oe/V),max(G-2,1))}
hl_chi<-function(y,ph,G){g<-grp(ph,G);o<-tapply(y,g,sum);e<-tapply(ph,g,sum);ng<-tapply(y,g,length);pb<-as.numeric(tapply(ph,g,mean))
  V<-ng*pb*(1-pb);oe<-as.numeric(o-e);1-pchisq(sum(oe^2/V),max(G-2,1))}
## ---- Aranda-Ordaz asymmetric link: pi(eta;alpha)=1-(1+alpha*e^eta)^(-1/alpha) ----
ao<-function(eta,alpha) 1-(1+alpha*exp(eta))^(-1/alpha)   # alpha=1 logit ; alpha->0 cloglog
sqb<-function(J)solve(rbind(c(1,-1.5,2.25),c(1,3,9),c(1,-3,9)),c(lg(.05),lg(.95),lg(J)))
sib<-function(I){e0<-lg(.1);e1<-lg(.2);e2<-lg(.2+I);b0<-(e0+e1)/2;b1<-(e1-e0)/6;b3<-(e2-(b0+3*b1))/6;c(b0,b1,3*b3,b3)}
## ---- one replication: returns EF,HL p-values for a scenario ----
gen<-function(fam,par,n,skew=FALSE){
  if(skew){x<-rchisq(n,4);eta<- -2.0+0.5*x;X<-data.frame(x=x);f<-y~x}
  else {x<-runif(n,-3,3);d<-rbinom(n,1,.5);eta<-0.6*x+0.5*d;X<-data.frame(x=x,d=d);f<-y~x+d}
  if(fam=="null") p<-plogis(eta)
  else if(fam=="ao") p<-ao(eta,par)
  else if(fam=="quad"){b<-sqb(par);x<-runif(n,-3,3);X<-data.frame(x=x);f<-y~x;p<-plogis(b[1]+b[2]*x+b[3]*x^2)}
  else if(fam=="int"){b<-sib(par);x<-runif(n,-3,3);d<-rbinom(n,1,.5);X<-data.frame(x=x,d=d);f<-y~x+d;p<-plogis(b[1]+b[2]*x+b[3]*d+b[4]*x*d)}
  X$y<-rbinom(n,1,p); list(D=X,f=f)
}
one<-function(fam,par,n,G,skew=FALSE){G0<-gen(fam,par,n,skew);fit<-tryCatch(suppressWarnings(glm(G0$f,data=G0$D,family=binomial())),error=function(e)NULL)
  if(is.null(fit))return(c(EF=NA,HL=NA));ph<-as.numeric(fitted(fit));y<-G0$D$y;c(EF=ef_chi(y,ph,G),HL=hl_chi(y,ph,G))}
## A(delta): large-n alignment functional for AO(alpha)
Adelta<-function(alpha,n=40000,runs=3,G=10){v<-numeric(runs)
  for(r in 1:runs){x<-runif(n,-3,3);d<-rbinom(n,1,.5);eta<-0.6*x+0.5*d;y<-rbinom(n,1,ao(eta,alpha))
    fit<-suppressWarnings(glm(y~x+d,family=binomial()));ph<-pmin(pmax(fitted(fit),1e-6),1-1e-6)
    g<-grp(ph,G);o<-tapply(y,g,sum);e<-tapply(ph,g,sum);ng<-tapply(y,g,length);pb<-as.numeric(tapply(ph,g,mean))
    V<-ng*pb*(1-pb);v[r]<-sum((1-2*pb)*as.numeric(o-e)/V)};mean(v)}

K<-5000; ncores<-max(1L,detectCores()-1L)
cl<-makeCluster(ncores); clusterSetRNGStream(cl,20260703)
clusterExport(cl,c("grp","ef_chi","hl_chi","ao","sqb","sib","gen","one","lg"))
run_cell<-function(fam,par,n,G,skew=FALSE){force(fam);force(par);force(n);force(G);force(skew)
  m<-parSapply(cl,1:K,function(i)one(fam,par,n,G,skew));c(EF=mean(m["EF",]<0.05,na.rm=TRUE),HL=mean(m["HL",]<0.05,na.rm=TRUE))}
res<-list(); t0<-Sys.time()
## ---- SIZE: g x covariate x n ----
for(G in c(5,10,20)) for(n in c(100,500,1000,5000)) for(sk in c(FALSE,TRUE)){
  r<-run_cell("null",NA,n,G,sk); res[[length(res)+1]]<-data.frame(block="size",fam="null",par=NA,n=n,G=G,skew=sk,EF=r["EF"],HL=r["HL"]) }
cat("size done",format(Sys.time()-t0),"\n")
## ---- POWER sweeps at G=10 ----
for(J in c(0.01,0.02,0.03,0.05,0.1)) for(n in c(500,1000,2000)){r<-run_cell("quad",J,n,10);res[[length(res)+1]]<-data.frame(block="power",fam="quad",par=J,n=n,G=10,skew=FALSE,EF=r["EF"],HL=r["HL"])}
for(I in c(0.1,0.3,0.5,0.7)) for(n in c(500,1000,2000)){r<-run_cell("int",I,n,10);res[[length(res)+1]]<-data.frame(block="power",fam="int",par=I,n=n,G=10,skew=FALSE,EF=r["EF"],HL=r["HL"])}
## ---- MONEY: Aranda-Ordaz asymmetry sweep (alpha x n) ----
for(al in c(0.9,0.7,0.5,0.35,0.2,0.1)) for(n in c(500,1000,2000,5000)){r<-run_cell("ao",al,n,10);res[[length(res)+1]]<-data.frame(block="power",fam="ao",par=al,n=n,G=10,skew=FALSE,EF=r["EF"],HL=r["HL"])}
cat("power done",format(Sys.time()-t0),"\n")
stopCluster(cl)
out<-do.call(rbind,res); rownames(out)<-NULL
write.csv(out,"definitive_results.csv",row.names=FALSE)
## A(delta) per alpha
ad<-data.frame(alpha=c(0.9,0.7,0.5,0.35,0.2,0.1)); ad$A<-sapply(ad$alpha,Adelta)
write.csv(ad,"definitive_Adelta.csv",row.names=FALSE)
cat("DONE",format(Sys.time()-t0),"-> definitive_results.csv, definitive_Adelta.csv\n")
print(out[out$fam=="ao" & out$n==1000,])
