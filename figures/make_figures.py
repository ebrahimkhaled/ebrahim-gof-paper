#!/usr/bin/env python
"""Publication-quality figures for the EF paper (Statistics in Medicine style).
Honest + compelling: shows EF tracks HL almost everywhere, with a visible edge only
for asymmetric-link misfit. Colourblind-safe (Okabe-Ito). Saves PDF + PNG."""
import os, csv, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

SIM = r"C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/simulations"
CODE = r"C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/EF_test_version_2/code"
OUT = [r"C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/EF_test_version_2/figures",
       r"C:/Users/ebrah/.gemini/Projects/PDFs/Paper_Ebrahim_Frangiton/EF_test_version_2/wiley_statistics_in_medicine/figures"]
for d in OUT: os.makedirs(d, exist_ok=True)

# Okabe-Ito colourblind-safe palette
C = dict(EF="#0072B2", HL="#D55E00", PH="#009E73", Tsiatis="#CC79A7", Xie="#56B4E9",
         grey="#888888", ok="#009E73", bad="#D55E00")
plt.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"], "mathtext.fontset": "dejavuserif",
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10, "legend.fontsize": 9,
    "xtick.labelsize": 9, "ytick.labelsize": 9, "axes.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 150,
    "savefig.dpi": 400, "savefig.bbox": "tight", "lines.linewidth": 1.8, "lines.markersize": 5,
})
def save(fig, name):
    for d in OUT:
        fig.savefig(os.path.join(d, name+".pdf")); fig.savefig(os.path.join(d, name+".png"))
    plt.close(fig); print("wrote", name)
def rd(path):
    return list(csv.DictReader(open(path, encoding="utf-8")))

# ---------- load ----------
pw = rd(os.path.join(SIM, "Power_chisq_EF_HL_PH.csv"))     # Scenario,x,n,EF,HL,PH
for r in pw:
    for k in ("EF","HL","PH"): r[k] = float(r[k])
    r["n"] = int(r["n"])
mf = rd(os.path.join(SIM, "Master_full.csv"))              # ,DEF..,EF,HL,HLeqw,PH,Tsiatis,Xie,PR,Stukel
fam = ["EF","HL","PH","Tsiatis","Xie"]
mrows = []
for r in mf:
    scn = r[list(r.keys())[0]] or r.get("")
    try: mrows.append((scn, {t: float(r[t]) for t in fam if r.get(t) not in (None,"")}))
    except Exception: pass
try:
    v1 = rd(os.path.join(CODE, "validate_V1.csv"))         # scn,G,A_delta,...,dPower,...
    for r in v1:
        r["A_delta"]=float(r["A_delta"]); r["dPower"]=float(r["dPower"])
except Exception: v1=[]

# ================= FIGURE 1: EF-vs-HL agreement scatter (the flagship) =================
fig, ax = plt.subplots(figsize=(4.6,4.4))
groups = {"Omitted quadratic":"Quad","Binary interaction":"Interact","Continuous interaction":"ContInteract"}
gcol = {"Omitted quadratic":C["grey"], "Binary interaction":"#66a3c2", "Continuous interaction":"#b0a0c0"}
ax.plot([0,1],[0,1], color="black", lw=1.0, ls=(0,(4,3)), zorder=1)
ax.fill_between([0,1],[0,1],[1,1], color=C["EF"], alpha=0.05, zorder=0)   # EF>HL region
ax.text(0.15,0.55,"EF more powerful", color=C["EF"], fontsize=8, rotation=45, alpha=.8)
for lab,scn in groups.items():
    xs=[r["HL"] for r in pw if r["Scenario"]==scn]; ys=[r["EF"] for r in pw if r["Scenario"]==scn]
    ax.scatter(xs,ys,s=26,color=gcol[lab],edgecolor="white",linewidth=.4,label=lab,zorder=3)
# link points: cloglog (the niche) highlighted, probit
cl=[r for r in pw if r["Scenario"]=="Link" and r["x"]=="cloglog"]
pr=[r for r in pw if r["Scenario"]=="Link" and r["x"]=="probit"]
ax.scatter([r["HL"] for r in cl],[r["EF"] for r in cl],s=95,marker="*",color=C["EF"],
           edgecolor="black",linewidth=.6,zorder=5,label="Asymmetric link (cloglog)")
ax.scatter([r["HL"] for r in pr],[r["EF"] for r in pr],s=34,marker="D",color=C["PH"],
           edgecolor="white",linewidth=.4,zorder=4,label="Symmetric link (probit)")
ax.set_xlim(0,1); ax.set_ylim(0,1); ax.set_aspect("equal")
ax.set_xlabel("Hosmer–Lemeshow power"); ax.set_ylabel("Ebrahim–Farrington power")
ax.legend(loc="lower right", frameon=False, fontsize=7.5, handletextpad=.4)
ax.set_title("EF tracks Hosmer–Lemeshow, except for asymmetric links", fontsize=10)
save(fig, "fig1_agreement")

# ================= FIGURE 2: power-profile small multiples =================
fig, axes = plt.subplots(1,3, figsize=(9.6,3.2), sharey=True)
panels=[("Quad","Omitted quadratic (severity $J$)"),("Interact","Binary interaction (severity $I$)"),
        ("ContInteract","Continuous interaction (severity $K$)")]
for ax,(scn,title) in zip(axes,panels):
    rows=sorted([r for r in pw if r["Scenario"]==scn and r["n"]==1000], key=lambda r:float(r["x"]))
    X=[float(r["x"]) for r in rows]
    ax.plot(X,[r["HL"] for r in rows], color=C["HL"], marker="s", ls="--", label="HL")
    ax.plot(X,[r["EF"] for r in rows], color=C["EF"], marker="o", label="EF")
    ax.plot(X,[r["PH"] for r in rows], color=C["PH"], marker="^", ls=":", label="PH", alpha=.9)
    ax.axhline(0.05, color=C["grey"], lw=.7, ls=":")
    ax.set_title(title, fontsize=9.5); ax.set_xlabel("severity"); ax.set_ylim(0,1.02)
axes[0].set_ylabel("power (n = 1000)")
axes[0].legend(loc="upper left", frameon=False)
axes[1].text(.5,.06,"EF and HL nearly coincide", ha="center", fontsize=8.5, color=C["grey"], transform=axes[1].transAxes)
fig.suptitle("EF is competitive with HL across nonlinear and interaction departures", y=1.02, fontsize=10.5)
save(fig, "fig2_power_profiles")

# ================= FIGURE 3: comprehensive comparison heatmap =================
order=["cloglog","loglog","scobit2","scobit_half","stukel_asym","cubic","threshold","sawtooth","osc4","bump",
       "quad","int_bin","int_cont","crossover","probit","cauchit","robit_t4","stukel_heavy","stukel_light"]
avail=[s for s in order if any(scn==s for scn,_ in mrows)]
M=np.array([[dict(mrows).get(s,{}).get(t,np.nan) for t in fam] for s in avail])
fig, ax = plt.subplots(figsize=(5.2, 6.4))
im=ax.imshow(M, aspect="auto", cmap="magma", vmin=0, vmax=1)
ax.set_xticks(range(len(fam))); ax.set_xticklabels(fam, rotation=30, ha="right")
ax.set_yticks(range(len(avail))); ax.set_yticklabels(avail, fontsize=8)
for i in range(len(avail)):
    jbest=int(np.nanargmax(M[i]));
    for j in range(len(fam)):
        if not np.isnan(M[i,j]):
            ax.text(j,i,f"{M[i,j]:.2f}",ha="center",va="center",fontsize=6.4,
                    color=("white" if M[i,j]<0.6 else "black"),
                    fontweight=("bold" if j==jbest else "normal"))
    if fam[jbest]=="EF": ax.add_patch(plt.Rectangle((-.5,i-.5),0.999,1,fill=False,edgecolor=C["EF"],lw=1.6))
cb=fig.colorbar(im,ax=ax,fraction=.046,pad=.04); cb.set_label("power (n = 1000, α = 0.05)")
ax.set_title("Partition-family power across the misspecification space\n(EF leads on asymmetric links; boxed = EF best)", fontsize=9.5)
save(fig, "fig3_heatmap")

# ================= FIGURE 4: A(delta) predicts the (small) EF-HL gap =================
if v1:
    fig, ax = plt.subplots(figsize=(4.8,4.0))
    A=np.array([r["A_delta"] for r in v1]); dP=np.array([r["dPower"] for r in v1])
    ax.axhline(0,color=C["grey"],lw=.7); ax.axvline(0,color=C["grey"],lw=.7)
    ax.fill_betweenx([-1,1],-8,0,color=C["EF"],alpha=.05);
    col=[C["EF"] if a<0 else C["HL"] for a in A]
    ax.scatter(-A,dP,s=34,c=col,edgecolor="white",linewidth=.4,zorder=3)
    for r in v1:
        if abs(r["dPower"])>0.05: ax.annotate(r["scn"],(-r["A_delta"],r["dPower"]),fontsize=6.5,
            xytext=(3,3),textcoords="offset points",color=C["grey"])
    ax.set_xlabel(r"$-A(\delta)=-\sum_g (1-2\bar\pi_g)\,\delta_g/V_g$")
    ax.set_ylabel("observed power(EF) − power(HL)")
    ax.set_title("A single functional predicts where EF beats HL\n(EF gains iff $A(\\delta)<0$: asymmetric departures)", fontsize=9)
    save(fig, "fig4_alignment")

# ================= FIGURE 5: Type-I calibration (EF = HL) =================
t1 = [("n=1000\nU(-3,3)",0.043,0.043),("n=5000\nU(-3,3)",0.043,0.045),("n=20000\nU(-3,3)",0.048,0.048),
      ("n=1000\nchi2(4)",0.065,0.072),("n=5000\nchi2(4)",0.075,0.072),("n=20000\nchi2(4)",0.055,0.055)]
fig, ax = plt.subplots(figsize=(6.0,3.2))
x=np.arange(len(t1)); w=.36
ax.bar(x-w/2,[r[1] for r in t1],w,color=C["EF"],label="EF")
ax.bar(x+w/2,[r[2] for r in t1],w,color=C["HL"],label="HL")
ax.axhline(0.05,color="black",lw=.9,ls="--"); ax.text(len(t1)-.5,0.052,"nominal 0.05",fontsize=7,ha="right")
ax.set_xticks(x); ax.set_xticklabels([r[0] for r in t1],fontsize=7.5)
ax.set_ylabel("Type-I error"); ax.set_ylim(0,0.09); ax.legend(frameon=False,loc="upper left")
ax.set_title("EF and HL are indistinguishable in size, at every sample size", fontsize=9.5)
save(fig, "fig5_typeI")
print("ALL FIGURES DONE ->", OUT[0])
