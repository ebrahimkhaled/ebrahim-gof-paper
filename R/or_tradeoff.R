con<-file("or_tradeoff.out","w"); w<-function(...) {cat(...,"\n",file=con);flush(con)}
suppressMessages(library(ebrahim.gof)); set.seed(20260704)
DGP<-list(
  null=function(n){x<-runif(n,-3,3);d<-rbinom(n,1,.5);data.frame(x=x,d=d,y=rbinom(n,1,plogis(0.6*x+0.5*d)))},
  cloglog=function(n){x<-runif(n,-3,3);d<-rbinom(n,1,.5);data.frame(x=x,d=d,y=rbinom(n,1,1-exp(-exp(0.6*x+0.5*d))))},
  oscillatory=function(n){x<-runif(n,-3,3);d<-rbinom(n,1,.5);data.frame(x=x,d=d,y=rbinom(n,1,plogis(0.6*x+0.5*d+1.5*sin(5*x))))},
  interaction=function(n){x<-runif(n,-3,3);d<-rbinom(n,1,.5);data.frame(x=x,d=d,y=rbinom(n,1,plogis(0.6*x+0.5*d+0.5*x*d)))})
getp<-function(bat,name){B<-as.data.frame(bat);i<-grep(name,B$Test,ignore.case=TRUE)[1];if(is.na(i))return(NA);as.numeric(B$p_value[i])}
orp<-function(scn,B,n=1000){v<-numeric(B)
  for(i in 1:B){D<-DGP[[scn]](n);f<-tryCatch(suppressWarnings(glm(y~x+d,data=D,family=binomial)),error=function(e)NULL);if(is.null(f)){v[i]<-NA;next}
    bat<-tryCatch(run.all.gof(f,include_slow=FALSE),error=function(e)NULL);v[i]<-if(is.null(bat))NA else getp(bat,"Osius")};v}
B<-500
nu<-orp("null",B);pstar<-quantile(nu,0.05,na.rm=TRUE)
w(sprintf("Ungrouped Osius-Rojek: raw size=%.3f  5%%-threshold p*=%.4f",mean(nu<.05,na.rm=T),pstar))
for(scn in c("cloglog","oscillatory","interaction")){p<-orp(scn,B)
  w(sprintf("  %-13s : raw power=%.3f  size-adj power=%.3f",scn,mean(p<.05,na.rm=T),mean(p<pstar,na.rm=T)))}
w("DONE");close(con)
