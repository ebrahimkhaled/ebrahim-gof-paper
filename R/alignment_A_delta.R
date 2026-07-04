con<-file("_ll.out","w"); set.seed(20260703)
grp<-function(ph,G=10){n<-length(ph);pmin(ceiling(rank(ph,ties.method="first")/(n/G)),G)}
Ad<-function(linkinv,n=100000){x<-runif(n,-3,3);d<-rbinom(n,1,.5);eta<-0.6*x+0.5*d;y<-rbinom(n,1,linkinv(eta))
  ph<-pmin(pmax(suppressWarnings(fitted(glm(y~x+d,family=binomial))),1e-6),1-1e-6);g<-grp(ph,10)
  o<-tapply(y,g,sum);e<-tapply(ph,g,sum);ng<-tapply(y,g,length);pb<-as.numeric(tapply(ph,g,mean));V<-ng*pb*(1-pb)
  sum((1-2*pb)*as.numeric(o-e)/V)}
cat(sprintf("cloglog A=%+.3f | loglog A=%+.3f\n", Ad(function(e)1-exp(-exp(e))), Ad(function(e)exp(-exp(-e)))),file=con)
close(con)
con<-file("_theory_check.out","w"); w<-function(...) {cat(...,"\n",file=con);flush(con)}
set.seed(20260703)
grp<-function(ph,G=10){n<-length(ph);pmin(ceiling(rank(ph,ties.method="first")/(n/G)),G)}
Adelta<-function(genf,n=200000,G=10){G0<-genf(n);y<-G0$y;ph<-pmin(pmax(G0$ph,1e-6),1-1e-6)
  g<-grp(ph,G);o<-tapply(y,g,sum);e<-tapply(ph,g,sum);ng<-tapply(y,g,length);pb<-as.numeric(tapply(ph,g,mean))
  V<-ng*pb*(1-pb);list(A=sum((1-2*pb)*as.numeric(o-e)/V), meanpi=mean(ph), fracAbovehalf=mean(pb>0.5))}
## MAIN design: truth 0.6x+0.5d + departure ; fit y~x+d
quad<-function(g) function(n){x<-runif(n,-3,3);d<-rbinom(n,1,.5);eta<-0.6*x+0.5*d+g*x^2;y<-rbinom(n,1,plogis(eta))
  list(y=y,ph=suppressWarnings(fitted(glm(y~x+d,family=binomial))))}
intr<-function(g) function(n){x<-runif(n,-3,3);d<-rbinom(n,1,.5);eta<-0.6*x+0.5*d+g*x*d;y<-rbinom(n,1,plogis(eta))
  list(y=y,ph=suppressWarnings(fitted(glm(y~x+d,family=binomial))))}
nullf<-function(n){x<-runif(n,-3,3);d<-rbinom(n,1,.5);eta<-0.6*x+0.5*d;y<-rbinom(n,1,plogis(eta))
  list(y=y,ph=suppressWarnings(fitted(glm(y~x+d,family=binomial))))}
w("=== A(delta) for QUADRATIC misfit (main design), sign predicts EF vs HL (A>0 => EF LOSES) ===")
for(g in c(0.05,0.1,0.2,0.3)){r<-Adelta(quad(g));w(sprintf("  quad gamma=%.2f: A(delta)=%+.3f  (mean pi=%.2f, %.0f%% deciles above 0.5)",g,r$A,r$meanpi,100*r$fracAbovehalf))}
w("=== A(delta) for INTERACTION misfit (main design) ===")
for(g in c(0.3,0.5,0.7)){r<-Adelta(intr(g));w(sprintf("  int gamma=%.1f: A(delta)=%+.3f",g,r$A))}
## SD(C) under the NULL at G=10, n=1000 (checks whether C=o_p(1): should be small, ~0.5-0.7)
w("=== null SD(C) and sum w_g^2/V_g (Var(C) leading term) at G=10 ===")
Cnull<-function(n,B=2000){cc<-numeric(B);vv<-numeric(B)
  for(i in 1:B){G0<-nullf(n);y<-G0$y;ph<-pmin(pmax(G0$ph,1e-6),1-1e-6);g<-grp(ph,10)
    o<-tapply(y,g,sum);e<-tapply(ph,g,sum);ng<-tapply(y,g,length);pb<-as.numeric(tapply(ph,g,mean));Vg<-ng*pb*(1-pb)
    cc[i]<-sum((1-2*pb)*as.numeric(o-e)/Vg); vv[i]<-sum((1-2*pb)^2/Vg)}
  c(sdC=sd(cc),meanC=mean(cc),sumw2V=mean(vv))}
for(n in c(500,1000,5000)){r<-Cnull(n);w(sprintf("  n=%4d: SD(C)=%.3f meanC=%.3f  sum w^2/V=%.3f (=Var(C) leading)",n,r["sdC"],r["meanC"],r["sumw2V"]))}
close(con)
