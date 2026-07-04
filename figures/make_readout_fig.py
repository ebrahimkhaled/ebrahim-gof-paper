#!/usr/bin/env python
"""fig_readout: the per-group directional READ-OUT demonstrated on a representative
cloglog-misfit dataset (Aranda-Ordaz alpha=0.1, n=1000; EF p=0.039<0.05<HL p=0.052).
Left  = signed per-group contributions c_g (bars) + the systematic mean tilt (line);
Right = observed-vs-expected decile calibration, one-sided bow = the directional misfit.
Same Okabe-Ito palette as the other figures."""
import os, csv, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
CODE=r"C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/EF_test_version_2/code"
OUT=[r"C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/EF_test_version_2/figures",
     r"C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/EF_test_version_2/wiley_statistics_in_medicine/figures"]
for d in OUT: os.makedirs(d,exist_ok=True)
C=dict(EF="#0072B2",HL="#D55E00"); POS="#009E73"; NEG="#D55E00"
plt.rcParams.update({"font.family":"serif","mathtext.fontset":"dejavuserif","font.size":10,
    "axes.titlesize":11,"axes.labelsize":10,"axes.linewidth":.8,"axes.spines.top":False,
    "axes.spines.right":False,"savefig.dpi":400,"savefig.bbox":"tight","figure.dpi":150})
def save(fig,name):
    for d in OUT: fig.savefig(os.path.join(d,name+".pdf")); fig.savefig(os.path.join(d,name+".png"))
    plt.close(fig); print("wrote",name)
R=[r for r in csv.DictReader(open(os.path.join(CODE,"readout_demo.csv")))]
g=[int(r["g"]) for r in R]; pbar=[float(r["pbar"]) for r in R]; o=[float(r["o"]) for r in R]
e=[float(r["e"]) for r in R]; ng=[float(r["ng"]) for r in R]; cg=[float(r["cg"]) for r in R]; mcg=[float(r["meancg"]) for r in R]

fig,(axL,axR)=plt.subplots(1,2,figsize=(9.4,3.9))
# ---- LEFT: signed c_g by decile + systematic mean tilt ----
cols=[POS if v>=0 else NEG for v in cg]
axL.bar(g,cg,color=cols,alpha=.85,width=.66,zorder=3,label="this dataset")
axL.plot(g,mcg,color="#333",lw=1.6,marker="o",ms=3.5,zorder=4,label="mean over 3000 reps")
axL.axhline(0,color="#555",lw=.9)
axL.set_xlabel("decile of fitted risk (low $\\to$ high)")
axL.set_ylabel(r"contribution $c_g=(1-2\bar\pi_g)(o_g-e_g)/V_g$")
axL.set_xticks(g)
axL.set_title("Directional read-out: a systematic tilt",fontsize=10)
axL.legend(frameon=False,fontsize=8.2,loc="lower left")
axL.annotate("under-prediction\nin low-risk groups",(1.6,0.30),fontsize=7.6,color="#00695c",ha="left")
axL.annotate("over-weighted\nhigh-risk misfit",(9.4,-0.55),fontsize=7.6,color="#a13d00",ha="right")
# ---- RIGHT: observed vs expected decile calibration ----
obs=[o[i]/ng[i] for i in range(len(g))]; exp=[e[i]/ng[i] for i in range(len(g))]
axR.plot([0,1],[0,1],color="#888",ls="--",lw=1,zorder=1)
axR.plot(exp,obs,color=C["EF"],marker="o",ms=5,lw=1.3,zorder=3)
axR.set_xlim(0,1); axR.set_ylim(0,1); axR.set_aspect("equal")
axR.set_xlabel(r"expected proportion $\bar\pi_g$ (logit fit)")
axR.set_ylabel("observed proportion $o_g/n_g$")
axR.set_title("Observed vs expected by decile",fontsize=10)
axR.text(0.04,0.93,"EF $p=0.039$\nHL $p=0.052$",fontsize=9,va="top",
    bbox=dict(boxstyle="round,pad=0.3",fc="white",ec="#bbb",lw=.7))
fig.suptitle(r"The read-out localises an asymmetric-link departure that HL's global statistic averages away "
             r"(cloglog misfit, $n=1000$)",y=1.02,fontsize=10)
save(fig,"fig_readout")
print("done")
