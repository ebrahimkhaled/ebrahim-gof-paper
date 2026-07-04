#!/usr/bin/env python
"""The 'crazy' figures from the definitive Aranda-Ordaz sweep:
 fig7 = EF-advantage LANDSCAPE contour over (link asymmetry x sample size);
 fig8 = money power small-multiples (EF vs HL, gap widening with asymmetry);
 fig4 = A(delta) predicts the advantage (clean sweep). Austin/Demler aesthetic."""
import os, csv, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import cm
CODE=r"C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/EF_test_version_2/code"
OUT=[r"C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/EF_test_version_2/figures",
     r"C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/EF_test_version_2/wiley_statistics_in_medicine/figures"]
for d in OUT: os.makedirs(d,exist_ok=True)
C=dict(EF="#0072B2",HL="#D55E00")
plt.rcParams.update({"font.family":"serif","mathtext.fontset":"dejavuserif","font.size":10,
    "axes.titlesize":11,"axes.labelsize":10,"axes.linewidth":.8,"axes.spines.top":False,
    "axes.spines.right":False,"savefig.dpi":400,"savefig.bbox":"tight","figure.dpi":150})
def save(fig,name):
    for d in OUT: fig.savefig(os.path.join(d,name+".pdf")); fig.savefig(os.path.join(d,name+".png"))
    plt.close(fig); print("wrote",name)
ao=[r for r in csv.DictReader(open(os.path.join(CODE,"definitive_results.csv"))) if r["fam"]=="ao"]
for r in ao: r["alpha"]=float(r["par"]); r["n"]=int(r["n"]); r["EF"]=float(r["EF"]); r["HL"]=float(r["HL"]); r["d"]=r["EF"]-r["HL"]
Ad={float(r["alpha"]):float(r["A"]) for r in csv.DictReader(open(os.path.join(CODE,"definitive_Adelta.csv")))}
alphas=sorted(set(r["alpha"] for r in ao), reverse=True)   # 0.9..0.1
ns=sorted(set(r["n"] for r in ao))                          # 500..5000
def cell(al,n):
    for r in ao:
        if r["alpha"]==al and r["n"]==n: return r
asym=[round(1-a,2) for a in alphas]                          # 0.1 (symmetric) .. 0.9 (cloglog)

# ================= FIG 7: EF-advantage LANDSCAPE contour =================
Z=np.array([[cell(a,n)["d"] for a in alphas] for n in ns])   # rows=n, cols=alpha
X=np.array([1-a for a in alphas]); Y=np.array(ns,dtype=float)
# smooth interpolation
try:
    from scipy.interpolate import griddata
    gx,gy=np.meshgrid(np.linspace(X.min(),X.max(),120), np.linspace(np.log10(Y.min()),np.log10(Y.max()),120))
    pts=np.array([(1-r["alpha"], np.log10(r["n"])) for r in ao]); vals=np.array([r["d"] for r in ao])
    gz=griddata(pts, vals, (gx,gy), method="cubic"); gz=np.clip(gz,0,None); yl=10**gy
    smooth=True
except Exception:
    smooth=False
fig,ax=plt.subplots(figsize=(6.4,4.8))
if smooth:
    cf=ax.contourf(gx,yl,gz,levels=np.linspace(0,0.07,15),cmap="YlGnBu",extend="max")
    cl=ax.contour(gx,yl,gz,levels=[0.01,0.02,0.03,0.04,0.05,0.06],colors="white",linewidths=.7)
    ax.clabel(cl,fmt="%.2f",fontsize=7,inline=True)
else:
    cf=ax.contourf(X,Y,Z,levels=15,cmap="YlGnBu")
ax.set_yscale("log"); ax.set_yticks(ns); ax.set_yticklabels([str(n) for n in ns])
ax.scatter([1-r["alpha"] for r in ao],[r["n"] for r in ao],s=10,c="black",alpha=.35,zorder=5)
ax.set_xlabel(r"link asymmetry  $1-\alpha$   (0 = logit $\;\to\;$ 1 = complementary log--log)")
ax.set_ylabel("sample size  $n$")
cb=fig.colorbar(cf,ax=ax,pad=.02); cb.set_label("power gain  power(EF) $-$ power(HL)")
ax.annotate("EF's advantage\nregion",(0.55,1500),fontsize=9,
    color="#08306b",ha="center",fontweight="bold")
ax.set_title("The Ebrahim–Farrington advantage over Hosmer–Lemeshow,\nmapped across link asymmetry and sample size",fontsize=10)
save(fig,"fig7_landscape")

# ================= FIG 8: money power small-multiples =================
fig,axes=plt.subplots(1,len(ns),figsize=(11,3.0),sharey=True)
for ax,n in zip(axes,ns):
    ef=[cell(a,n)["EF"] for a in alphas]; hl=[cell(a,n)["HL"] for a in alphas]
    ax.fill_between(asym,hl,ef,color=C["EF"],alpha=.15,zorder=1)
    ax.plot(asym,hl,color=C["HL"],marker="s",ls="--",label="HL",zorder=3)
    ax.plot(asym,ef,color=C["EF"],marker="o",label="EF",zorder=4)
    ax.axhline(0.05,color="#888",lw=.7,ls=":")
    ax.set_title(f"$n = {n}$",fontsize=10); ax.set_xlabel(r"asymmetry $1-\alpha$"); ax.set_ylim(0,1.02)
axes[0].set_ylabel("power ($\\alpha=0.05$)"); axes[0].legend(frameon=False,loc="upper left")
fig.suptitle("Power to detect a misspecified (Aranda–Ordaz) link: the EF advantage (shaded) grows with asymmetry, then recedes as $n\\to\\infty$",y=1.03,fontsize=10)
save(fig,"fig8_sweep")

# ================= FIG 4 (update): A(delta) predicts the advantage =================
fig,ax=plt.subplots(figsize=(5.0,4.2))
cmap=cm.get_cmap("viridis");
for i,n in enumerate(ns):
    xs=[-Ad[a] for a in alphas]; ys=[cell(a,n)["d"] for a in alphas]
    ax.plot(xs,ys,"-o",color=cmap(i/(len(ns)-1)),label=f"n={n}",ms=5,lw=1.3)
ax.axhline(0,color="#888",lw=.7)
ax.set_xlabel(r"$-A(\delta)$   (alignment functional; larger = more directional misfit)")
ax.set_ylabel("power(EF) − power(HL)")
ax.legend(frameon=False,title="sample size",fontsize=8.5)
ax.set_title(r"A single functional $A(\delta)$ predicts EF's advantage:"+"\n"+r"the gain rises monotonically with $-A(\delta)$",fontsize=9.5)
save(fig,"fig4_alignment")
print("done")
