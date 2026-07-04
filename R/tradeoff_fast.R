con<-file("tradeoff_fast.out","w"); w<-function(...) {cat(...,"\n",file=con);flush(con)}
set.seed(20260704)
grp<-function(ph,G=10){n<-length(ph);pmin(ceiling(rank(ph,ties.method="first")/(n/G)),G)}
rp<-function(y,ph,ef=TRUE){ph<-pmin(pmax(ph,1e-6),1-1e-6);g<-grp(ph);o<-tapply(y,g,sum);e<-tapply(ph,g,sum);ng<-tapply(y,g,length);pb<-as.numeric(tapply(ph,g,mean));V<-ng*pb*(1-pb);oe<-as.numeric(o-e);s<-sum(oe^2/V);if(ef)s<-s-sum((1-2*pb)*oe/V);1-pchisq(s,8)}
stukel<-function(D,f){eta<-predict(f);z1<-ifelse(eta>=0,eta^2/2,0);z2<-ifelse(eta<0,-eta^2/2,0)
  a<-tryCatch(suppressWarnings(glm(y~x+d+z1+z2,data=cbind(D,z1=z1,z2=z2),family=binomial)),error=function(e)NULL)
  if(is.null(a)||!a$converged)return(NA);1-pchisq(deviance(f)-deviance(a),2)}
DGP<-list(
  null=function(n){x<-runif(n,-3,3);d<-rbinom(n,1,.5);data.frame(x=x,d=d,y=rbinom(n,1,plogis(0.6*x+0.5*d)))},
  cloglog=function(n){x<-runif(n,-3,3);d<-rbinom(n,1,.5);data.frame(x=x,d=d,y=rbinom(n,1,1-exp(-exp(0.6*x+0.5*d))))},
  oscillatory=function(n){x<-runif(n,-3,3);d<-rbinom(n,1,.5);data.frame(x=x,d=d,y=rbinom(n,1,plogis(0.6*x+0.5*d+1.5*sin(5*x))))},
  interaction=function(n){x<-runif(n,-3,3);d<-rbinom(n,1,.5);data.frame(x=x,d=d,y=rbinom(n,1,plogis(0.6*x+0.5*d+0.5*x*d)))})
pvals<-function(scn,B,n=1000){m<-matrix(NA,B,3);colnames(m)<-c("EF","HL","Stukel")
  for(i in 1:B){D<-DGP[[scn]](n);f<-suppressWarnings(glm(y~x+d,data=D,family=binomial));ph<-as.numeric(fitted(f))
    m[i,]<-c(rp(D$y,ph,TRUE),rp(D$y,ph,FALSE),stukel(D,f))};m}
B<-3000
nullp<-pvals("null",B); pstar<-apply(nullp,2,quantile,0.05,na.rm=TRUE)
w(sprintf("raw size @0.05: EF=%.3f HL=%.3f Stukel=%.3f | 5%% p-thresholds: EF=%.4f HL=%.4f Stukel=%.4f",
  mean(nullp[,1]<.05),mean(nullp[,2]<.05),mean(nullp[,3]<.05,na.rm=T),pstar[1],pstar[2],pstar[3]))
w("\nscenario         | RAW power (a=.05)      | SIZE-ADJUSTED power")
w("                 | EF    HL    Stukel     | EF    HL    Stukel")
for(scn in c("cloglog","oscillatory","interaction")){P<-pvals(scn,B)
  raw<-colMeans(P<.05,na.rm=T); adj<-sapply(1:3,function(j)mean(P[,j]<pstar[j],na.rm=T))
  w(sprintf("%-16s | %.2f  %.2f  %.2f     | %.2f  %.2f  %.2f",scn,raw[1],raw[2],raw[3],adj[1],adj[2],adj[3]))}
## Stukel refit-fail vs sample size (sparsity fragility)
w("\nStukel refit-failure rate by n (cloglog data):")
for(n in c(1000,200,100,60,40)){f2<-0;B2<-400
  for(i in 1:B2){D<-DGP[["cloglog"]](n);ff<-suppressWarnings(glm(y~x+d,data=D,family=binomial));if(is.na(stukel(D,ff)))f2<-f2+1}
  w(sprintf("  n=%4d : refit-fail=%.1f%%",n,100*f2/B2))}
w("DONE"); close(con)
