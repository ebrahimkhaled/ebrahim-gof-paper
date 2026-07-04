con<-file("_s6b.out","w"); w<-function(...) {cat(...,"\n",file=con);flush(con)}
suppressMessages(library(MASS))
grp<-function(ph,G=10){n<-length(ph);pmin(ceiling(rank(ph,ties.method="first")/(n/G)),G)}
info<-function(y,ph,G=10){ph<-pmin(pmax(ph,1e-6),1-1e-6);g<-grp(ph,G)
  o<-as.numeric(tapply(y,g,sum));e<-as.numeric(tapply(ph,g,sum));ng<-as.numeric(tapply(y,g,length))
  pb<-as.numeric(tapply(ph,g,mean));V<-ng*pb*(1-pb);cvar<-as.numeric(tapply(ph*(1-ph),g,sum));Gr<-length(o)
  HL<-sum((o-e)^2/V);C<-sum((1-2*pb)*(o-e)/V);EF<-HL-C;PH<-sum((o-e)^2/cvar)
  list(EF_p=1-pchisq(EF,Gr-2),HL_p=1-pchisq(HL,Gr-2),PH_p=1-pchisq(PH,Gr-1),C=C,Gr=Gr,
       tab=data.frame(g=1:Gr,pbar=round(pb,3),ng=ng,o=o,e=round(e,1),cg=round((1-2*pb)*(o-e)/V,3)))}
rep<-function(tag,y,ph,G=10){s<-info(y,ph,G);w(sprintf("%-32s n=%3d ev=%.2f npat=%4d | EF=%.3f HL=%.3f PH=%.3f C=%+.2f (G=%d)",
  tag,length(y),mean(y),length(unique(round(ph,8))),s$EF_p,s$HL_p,s$PH_p,s$C,s$Gr));invisible(s)}
data(birthwt);d<-birthwt;d$low<-as.integer(d$low);y<-d$low
w("=== continuous-covariate benchmarks (proper decile regime) ===")
s_lbw<-rep("LBW additive",y,fitted(glm(low~age+lwt+smoke+ptl+ht+ui+factor(race),data=d,family=binomial)))
rep("LBW +age:lwt,smoke:lwt [fix]",y,fitted(glm(low~age+lwt+smoke+ptl+ht+ui+factor(race)+age:lwt+smoke:lwt,data=d,family=binomial)))
try({data(kyphosis,package="rpart");k<-kyphosis;ky<-as.integer(k$Kyphosis=="present")
  rep("Kyphosis additive",ky,fitted(glm(ky~Age+Number+Start,data=k,family=binomial)))
  rep("Kyphosis +poly(Age,2) [fix]",ky,fitted(glm(ky~poly(Age,2)+Number+Start,data=k,family=binomial)))},silent=TRUE)
try({ic<-aplore3::icu;iy<-as.integer(ic$sta=="Died")
  rep("ICU age+sys+typ+loc",iy,fitted(glm(iy~age+sys+typ+loc,data=ic,family=binomial)))},silent=TRUE)
try({gl<-aplore3::glow500;gy<-as.integer(gl$fracture=="Yes")
  rep("GLOW age+weight+priorfrac..",gy,fitted(glm(gy~age+weight+priorfrac+premeno+raterisk,data=gl,family=binomial)))},silent=TRUE)
## Pima external validation (parsimonious, well-specified) for the application
data(Pima.tr);data(Pima.te);Pima.tr$yy<-as.integer(Pima.tr$type=="Yes");Pima.te$yy<-as.integer(Pima.te$type=="Yes")
dev<-glm(yy~npreg+glu+bp+skin+bmi+ped+age,data=Pima.tr,family=binomial)
p.te<-plogis(predict(dev,newdata=Pima.te))
rep("Pima external (parsimonious)",Pima.te$yy,p.te)
w(""); w("=== LBW additive per-decile read-out (c_g) ===")
for(i in 1:nrow(s_lbw$tab)){r<-s_lbw$tab[i,];w(sprintf("  g%2d pbar=%.3f n=%2d o=%2d e=%4.1f cg=%+.3f",r$g,r$pbar,r$ng,r$o,r$e,r$cg))}
write.csv(s_lbw$tab,"lbw_readout.csv",row.names=FALSE)
w("wrote lbw_readout.csv")
close(con)
