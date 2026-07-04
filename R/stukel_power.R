con<-file("_stukel2.out","w"); set.seed(20260703)
grp<-function(ph,G=10){n<-length(ph);pmin(ceiling(rank(ph,ties.method="first")/(n/G)),G)}
raw_ef<-function(y,ph){ph<-pmin(pmax(ph,1e-6),1-1e-6);g<-grp(ph);o<-tapply(y,g,sum);e<-tapply(ph,g,sum);ng<-tapply(y,g,length);pb<-as.numeric(tapply(ph,g,mean));V<-ng*pb*(1-pb);oe<-as.numeric(o-e);1-pchisq(sum(oe^2/V)-sum((1-2*pb)*oe/V),8)}
raw_hl<-function(y,ph){ph<-pmin(pmax(ph,1e-6),1-1e-6);g<-grp(ph);o<-tapply(y,g,sum);e<-tapply(ph,g,sum);ng<-tapply(y,g,length);pb<-as.numeric(tapply(ph,g,mean));V<-ng*pb*(1-pb);1-pchisq(sum((o-e)^2/V),8)}
stukel<-function(D,base){eta<-predict(base);z1<-ifelse(eta>=0,eta^2/2,0);z2<-ifelse(eta<0,-eta^2/2,0)
  aug<-tryCatch(suppressWarnings(glm(y~x+d+z1+z2,data=cbind(D,z1=z1,z2=z2),family=binomial)),error=function(e)NULL)
  if(is.null(aug)||!aug$converged) return(NA); 1-pchisq(deviance(base)-deviance(aug),2)}
# EXACT complementary log-log truth, fit logit (matches tab:family "cloglog link")
gen<-function(n=1000){x<-runif(n,-3,3);d<-rbinom(n,1,.5);eta<-0.6*x+0.5*d;y<-rbinom(n,1,1-exp(-exp(eta)));data.frame(x=x,d=d,y=y)}
K<-2000; ef<-hl<-st<-numeric(K);fail<-0
for(i in 1:K){D<-gen();f<-suppressWarnings(glm(y~x+d,data=D,family=binomial));ph<-as.numeric(fitted(f))
  ef[i]<-raw_ef(D$y,ph);hl[i]<-raw_hl(D$y,ph);sp<-stukel(D,f);if(is.na(sp))fail<-fail+1;st[i]<-sp}
cat(sprintf("EXACT cloglog n=1000 (K=%d): EF=%.3f HL=%.3f Stukel=%.3f (Stukel refit-fail=%.3f)\n",K,mean(ef<.05),mean(hl<.05),mean(st<.05,na.rm=TRUE),fail/K),file=con)
close(con)
