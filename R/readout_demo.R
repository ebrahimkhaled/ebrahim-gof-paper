con<-file("readout_demo.out","w"); w<-function(...) {cat(...,"\n",file=con);flush(con)}
set.seed(20260703)
ao<-function(eta,alpha) 1-(1+alpha*exp(eta))^(-1/alpha)   # paper's Aranda-Ordaz; alpha->0 = cloglog
grp<-function(ph,G=10){n<-length(ph);pmin(ceiling(rank(ph,ties.method="first")/(n/G)),G)}
ro<-function(y,ph,G=10){ph<-pmin(pmax(ph,1e-6),1-1e-6);g<-grp(ph,G)
  o<-as.numeric(tapply(y,g,sum));e<-as.numeric(tapply(ph,g,sum));ng<-as.numeric(tapply(y,g,length))
  pb<-as.numeric(tapply(ph,g,mean));V<-ng*pb*(1-pb);cg<-(1-2*pb)*(o-e)/V
  list(g=1:G,pbar=pb,ng=ng,o=o,e=e,cg=cg,HL=sum((o-e)^2/V),EF=sum((o-e)^2/V)-sum(cg),C=sum(cg))}
pv<-function(s,df=8)1-pchisq(s,df)
gen<-function(n=1000,alpha=0.1){x<-runif(n,-3,3);d<-rbinom(n,1,.5);eta<-0.6*x+0.5*d;y<-rbinom(n,1,ao(eta,alpha))
  ph<-suppressWarnings(fitted(glm(y~x+d,family=binomial)));list(y=y,ph=ph)}
K<-3000; efp<-hzp<-numeric(K); cgmat<-matrix(NA,K,10); cand<-list()
for(i in 1:K){d<-gen();r<-ro(d$y,d$ph);efp[i]<-pv(r$EF);hzp[i]<-pv(r$HL);cgmat[i,]<-r$cg
  if(efp[i]<0.05 && hzp[i]>0.05) cand[[length(cand)+1]]<-list(i=i,r=r,ef=efp[i],hz=hzp[i])}
mc<-colMeans(cgmat)
w(sprintf("AO alpha=0.1 (cloglog) misfit, n=1000, K=%d: EF rejects %.1f%%, HL rejects %.1f%% (Delta=%+.1fpp)",K,100*mean(efp<.05),100*mean(hzp<.05),100*(mean(efp<.05)-mean(hzp<.05))))
w(sprintf("mean C = %+.3f (negative => EF systematically > HL, matching the niche)",mean(rowSums(cgmat))))
w("mean c_g by decile (SYSTEMATIC read-out tilt):")
w(paste(sprintf("g%d:%+.3f",1:10,mc),collapse="  "))
## pick the candidate (EF<0.05<HL) whose c_g pattern is CLOSEST to the mean tilt (most representative, not extreme)
if(length(cand)>0){
  dist<-sapply(cand,function(c) sum((c$r$cg-mc)^2)); sel<-cand[[which.min(dist)]]; r<-sel$r
  w(sprintf("\n#candidates EF<.05<HL = %d ; REPRESENTATIVE = rep %d (closest c_g to mean): EF p=%.4f, HL p=%.4f, C=%+.2f",length(cand),sel$i,sel$ef,sel$hz,r$C))
  for(k in 1:10) w(sprintf("  g%2d pbar=%.3f n=%3d o=%3d e=%5.1f cg=%+.3f",k,r$pbar[k],r$ng[k],r$o[k],r$e[k],r$cg[k]))
  write.csv(data.frame(g=r$g,pbar=r$pbar,ng=r$ng,o=r$o,e=r$e,cg=r$cg,meancg=mc),"readout_demo.csv",row.names=FALSE)
  w("wrote readout_demo.csv")
}
close(con)
