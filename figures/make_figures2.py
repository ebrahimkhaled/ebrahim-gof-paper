#!/usr/bin/env python
"""Upgraded figures: (1) fig1 refocused on EF's winning scenarios + real-world labels,
(2) fig3 heatmap yellow->green, (3) a striking diverging EF-HL 'who wins' bar chart."""
import os, csv, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
SIM=r"C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/simulations"
OUT=[r"C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/EF_test_version_2/figures",
     r"C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/EF_test_version_2/wiley_statistics_in_medicine/figures"]
for d in OUT: os.makedirs(d,exist_ok=True)
C=dict(EF="#0072B2",HL="#D55E00",win="#1b9e77",tie="#999999",loss="#d95f02")
plt.rcParams.update({"font.family":"serif","mathtext.fontset":"dejavuserif","font.size":10,
    "axes.titlesize":11,"axes.labelsize":10,"axes.linewidth":.8,"axes.spines.top":False,
    "axes.spines.right":False,"savefig.dpi":400,"savefig.bbox":"tight","figure.dpi":150})
def save(fig,name):
    for d in OUT: fig.savefig(os.path.join(d,name+".pdf")); fig.savefig(os.path.join(d,name+".png"))
    plt.close(fig); print("wrote",name)
# ---- load Master_full (19 scenarios, n=1000) ----
rows=list(csv.reader(open(os.path.join(SIM,"Master_full.csv"),encoding="utf-8")))
hdr=rows[0]; idx={h:i for i,h in enumerate(hdr)}
def val(r,t):
    v=r[idx[t]] if idx.get(t) is not None else "";
    return float(v) if v not in ("",None) else np.nan
data={}
for r in rows[1:]:
    scn=r[0]
    if scn in ("null",""): continue
    data[scn]=(val(r,"EF"),val(r,"HL"))
# pretty names + real-world context for the winners
NICE={"cloglog":"cloglog link","loglog":"loglog link","cubic":"cubic","threshold":"threshold",
 "sawtooth":"sawtooth","osc4":"oscillatory","bump":"bump","skew":"skewed","stukel_asym":"Stukel asym.",
 "quad":"quadratic","int_bin":"interaction","int_cont":"cont. interaction","crossover":"crossover","int_binbin":"interaction (bin.×bin.)",
 "probit":"probit","cauchit":"cauchit","robit_t4":"robit","scobit2":"scobit","scobit_half":"scobit-½",
 "stukel_heavy":"Stukel heavy","stukel_light":"Stukel light","joint":"joint","corr":"correlated","logx":"log-x"}
REAL={"cloglog":"dose–response toxicity;\ndiscrete-time survival; rare events",
      "loglog":"left-skewed risk\n(complementary asymmetry)"}

# ================= FIG 1 v2: EF wins where the link is asymmetric =================
fig,ax=plt.subplots(figsize=(5.6,5.2))
ax.fill_between([0,1],[0,1],[1,1],color=C["EF"],alpha=.06,zorder=0)
ax.plot([0,1],[0,1],color="black",lw=1.0,ls=(0,(4,3)),zorder=1)
ax.text(.30,.42,"EF more powerful",color=C["EF"],fontsize=8.5,rotation=45,alpha=.65)
ax.text(.64,.57,"EF = Hosmer–Lemeshow",color=C["tie"],fontsize=8,rotation=45,alpha=.85)
for scn,(ef,hl) in data.items():
    if np.isnan(ef) or np.isnan(hl): continue
    d=ef-hl
    if d>0.025:  ax.scatter(hl,ef,s=150,marker="*",color=C["win"],edgecolor="black",lw=.7,zorder=6)
    elif d<-0.03: ax.scatter(hl,ef,s=34,color=C["loss"],edgecolor="white",lw=.4,zorder=4)
    else: ax.scatter(hl,ef,s=30,color=C["tie"],edgecolor="white",lw=.4,alpha=.75,zorder=3)
# label the two winners compactly + one real-world box
for scn in ("cloglog","loglog"):
    if scn in data:
        ef,hl=data[scn]
        ax.annotate(NICE.get(scn,scn),(hl,ef),xytext=(9,4),textcoords="offset points",
            fontsize=8,color=C["win"],fontweight="bold")
ax.text(0.03,0.79,"EF's niche: asymmetric-link misfit\n • dose–response toxicity\n • discrete-time survival / hazards\n • rare, skewed events",
    fontsize=7.6,color=C["win"],va="top",
    bbox=dict(boxstyle="round,pad=0.4",fc="white",ec=C["win"],lw=.8,alpha=.92))
ax.set_xlim(0,1);ax.set_ylim(0,1);ax.set_aspect("equal")
ax.set_xlabel("Hosmer–Lemeshow power");ax.set_ylabel("Ebrahim–Farrington power")
leg=[plt.Line2D([],[],marker="*",ls="",ms=12,mfc=C["win"],mec="black",label="EF more powerful (asymmetric link)"),
     plt.Line2D([],[],marker="o",ls="",ms=6,mfc=C["tie"],label="EF ≈ HL (symmetric misfit)"),
     plt.Line2D([],[],marker="o",ls="",ms=6,mfc=C["loss"],label="EF < HL")]
ax.legend(handles=leg,loc="lower right",frameon=False,fontsize=7.8)
ax.set_title("Where the Ebrahim–Farrington test wins:\nasymmetric-link misfit, common in survival and dose–response",fontsize=9.5)
save(fig,"fig1_agreement")

# ================= FIG 3 v2: heatmap yellow -> green =================
fam=["EF","HL","PH","Tsiatis","Xie"]
order=["cloglog","loglog","scobit2","scobit_half","stukel_asym","cubic","threshold","sawtooth","osc4","bump",
       "quad","int_bin","int_cont","crossover","probit","cauchit","robit_t4","stukel_heavy","stukel_light"]
avail=[s for s in order if s in data]
def rowvals(r,scn):
    for rr in rows[1:]:
        if rr[0]==scn: return [ (float(rr[idx[t]]) if idx.get(t) and rr[idx[t]] not in("",None) else np.nan) for t in fam]
    return [np.nan]*len(fam)
M=np.array([rowvals(rows,s) for s in avail])
fig,ax=plt.subplots(figsize=(5.2,6.6))
im=ax.imshow(M,aspect="auto",cmap="YlGn",vmin=0,vmax=1)
ax.set_xticks(range(len(fam)));ax.set_xticklabels(fam,rotation=30,ha="right")
ax.set_yticks(range(len(avail)));ax.set_yticklabels([NICE.get(s,s) for s in avail],fontsize=8)
for i in range(len(avail)):
    jb=int(np.nanargmax(M[i]))
    for j in range(len(fam)):
        if not np.isnan(M[i,j]):
            ax.text(j,i,f"{M[i,j]:.2f}",ha="center",va="center",fontsize=6.4,
                    color=("black" if M[i,j]>0.4 else "#333"),fontweight=("bold" if j==jb else "normal"))
    if fam[jb]=="EF": ax.add_patch(plt.Rectangle((-.5,i-.5),.999,1,fill=False,edgecolor=C["EF"],lw=1.8))
cb=fig.colorbar(im,ax=ax,fraction=.046,pad=.04);cb.set_label("power (n = 1000)")
ax.set_title("Power across the misspecification space\n(blue box = EF best in family)",fontsize=9.5)
save(fig,"fig3_heatmap")

# ================= FIG 6: diverging 'who wins' bar chart (crazy + honest) =================
deltas=sorted([(s,data[s][0]-data[s][1]) for s in data if not np.isnan(data[s][0]) and not np.isnan(data[s][1])],
              key=lambda kv:kv[1])
labels=[NICE.get(s,s) for s,_ in deltas]; d=[v for _,v in deltas]
fig,ax=plt.subplots(figsize=(6.4,6.2))
cols=[C["win"] if v>0.02 else (C["loss"] if v<-0.02 else C["tie"]) for v in d]
ax.barh(range(len(d)),d,color=cols,edgecolor="white",height=.72)
ax.axvline(0,color="black",lw=.9)
ax.set_yticks(range(len(d)));ax.set_yticklabels(labels,fontsize=8)
ax.set_xlabel("power(EF) − power(HL)   (n = 1000)")
for i,(s,v) in enumerate(deltas):
    if abs(v)>0.02: ax.text(v+0.004,i,f"{v:+.02f}",va="center",ha="left",fontsize=7,color="#222")
ax.text(0.045,len(d)-3.2,"EF wins\n(asymmetric links)",color=C["win"],fontsize=9,fontweight="bold")
ax.text(-0.155,3.2,"EF behind\n(symmetric curvature,\ncovariate-space)",color=C["loss"],fontsize=8)
ax.set_xlim(-0.18,0.14)
ax.set_title("Honest ledger: EF equals HL almost everywhere,\nand wins precisely where the link is asymmetric",fontsize=9.5)
save(fig,"fig6_whowins")
print("done")
