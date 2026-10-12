"""make_esm_tables.py -- the tables of Online Resource 1, written from the result files (never typed):
S1-S4 from the markdown tables of theory/SPARSE_GROUP_THEOREM.md (written by theory/sparse_sim), S5 from
theory/sparse_sim/out/corrupt_summary.csv, S6 from results/realdata_screen.csv (beetle row replaced by the
one-group-per-dose analysis in results/beetle.rds, read through Rscript).

    python make_esm_tables.py
"""
import csv
import os
import re
import subprocess

ROOT = r"."
OUT = os.path.join(ROOT, "sp", "esm", "tables")
os.makedirs(OUT, exist_ok=True)
md = open(os.path.join(ROOT, "theory", "SPARSE_GROUP_THEOREM.md"), encoding="utf-8").read()


def md_table(tag):
    """rows of the markdown table under the heading '### <tag>.'"""
    block = md.split("### " + tag + ".", 1)[1]
    rows, started = [], False
    for l in block.splitlines():
        if l.startswith("|"):
            rows.append(l); started = True
        elif started:
            break
    out = []
    for l in rows[2:]:
        out.append([c.strip() for c in l.strip().strip("|").split("|")])
    return out


def num(s):
    s = s.strip()
    return re.sub(r"(^|[ /])-([0-9])", r"\1$-$\2", s)


def write(name, colspec, header, body):
    lines = ["\\begin{tabular}{@{}" + colspec + "@{}}", "\\toprule", header + " \\\\", "\\midrule"]
    lines += body
    lines += ["\\bottomrule", "\\end{tabular}", ""]
    open(os.path.join(OUT, name), "w", encoding="utf-8").write("\n".join(lines))
    print("wrote", name, len(body), "rows")


# S1: null, fitted-risk grouping
body, last = [], None
for r in md_table("T1"):
    n, m, G, sdEF, sdHL, mEF, mC, zEF, zHL, azHL, hlc, efc, sk = r
    if last is not None and n != last:
        body.append("\\addlinespace")
    last = n
    body.append(" & ".join(num(x) for x in [n, m, G, sdEF, sdHL, mEF, mC, zEF, zHL, hlc, efc]) + " \\\\")
write("tab_S1.tex", "rrrcccrrrrr",
      "$n$ & $m$ & $G$ & sd $\\EF$ sim/pred & sd $\\HL$ sim/pred & mean $\\EF-G$ sim/pred & mean $C$ & $Z_{\\EF}$ & $Z_{\\HL}$ & $\\HL$, $\\chi^2_{G-2}$ & $\\EF$, $\\chi^2_{G-2}$",
      body)

# S2: null, true-risk grouping
body, last = [], None
for r in md_table("T1o"):
    if last is not None and r[0] != last:
        body.append("\\addlinespace")
    last = r[0]
    body.append(" & ".join(num(x) for x in r) + " \\\\")
write("tab_S2.tex", "rrccrr",
      "$n$ & $m$ & sd $\\EF$ sim/pred & sd $\\HL$ sim/pred & $Z_{\\EF}$ & $Z_{\\HL}$", body)

# S3: variance ratio, wide risk range
body = []
for r in md_table("T2"):
    n, m, R, emp, plug, form, zEF, zHL, hlc = r
    body.append(" & ".join(num(x) for x in [m, emp, plug, form, zEF, zHL, hlc]) + " \\\\")
write("tab_S3.tex", "rrrrrrr",
      "$m$ & simulated & plug-in & $1+R/\\{2(m-1)\\}$ & $Z_{\\EF}$ & $Z_{\\HL}$ & $\\HL$, $\\chi^2_{G-2}$", body)

# S4: power
lab = {"cloglog": "cloglog", "loglog": "log-log", "quadneg": "square ($-$)", "quadpos": "square ($+$)"}
body, last = [], None
for r in md_table("T3"):
    alt, n, m, D, A, ef, hl, ahl, hlc, zl, fis = r
    if last is not None and alt != last:
        body.append("\\addlinespace")
    first = alt != last
    last = alt
    body.append(" & ".join([lab[alt] if first else ""] + [num(x) for x in [n, m, D, A, ef, hl, hlc, zl]]) + " \\\\")
write("tab_S4.tex", "lrrrrccrr",
      "Departure & $n$ & $m$ & $\\Delta_n$ & $A_n$ & $\\EF$ sim/pred & $\\HL$ sim/pred & $\\HL$, $\\chi^2_{G-2}$ & $Z_L$", body)

# S5: corruption with many small groups
body, last = [], None
for r in csv.DictReader(open(os.path.join(ROOT, "theory", "sparse_sim", "out", "corrupt_summary.csv"))):
    if last is not None and r["m"] != last:
        body.append("\\addlinespace")
    last = r["m"]
    err = "none" if r["k"] == "0" else ("$\\times4$" if float(r["mult"]) > 0 else "$\\times(-4)$")
    body.append(" & ".join([r["m"], err, r["k"]] + ["%.3f" % float(r[c]) for c in ("EF", "HL", "ZL_two_sided")]) + " \\\\")
write("tab_S5.tex", "rlrrrr", "$m$ & Error & $k$ & $Z_{\\EF}$ & $Z_{\\HL}$ & $Z_L$ (two-sided)", body)

# S6: public data sets; beetle by dose
beetle = subprocess.run(["Rscript", "-e",
                         "B <- readRDS('" + os.path.join(ROOT, "results", "beetle.rds").replace("\\", "/") + "');"
                         "cat(B$t[['HL']], B$t[['C']], B$t[['EF']], B$t[['p_HL_chisq']], B$t[['p_EF_chisq']])"],
                        capture_output=True, text=True).stdout.split()
body = []
for r in csv.DictReader(open(os.path.join(ROOT, "results", "realdata_screen.csv"), encoding="utf-8")):
    HL, C, EF, pHL, pEF = (float(x) for x in (beetle if r["data"] == "beetle" else
                                              (r["HL"], r["C"], r["EF"], r["p_HL"], r["p_EF"])))
    name = r["data"] + (" (by dose)" if r["data"] == "beetle" else "")
    src = r["source"].split(",")[-1].strip().replace("::", "::")
    body.append(" & ".join([name, r["n"], r["events"], num("%.2f" % float(r["AIC_logit_minus_cloglog"])),
                            "%.2f" % HL, num("%.2f" % C), "%.2f" % EF, "%.3f" % pHL, "%.3f" % pEF,
                            "%.3f" % float(r["p_Stukel"]), "\\texttt{" + src.replace("_", "\\_") + "}"]) + " \\\\")
write("tab_S6.tex", "lrrrrrrrrrl",
      "Data & $n$ & events & $\\Delta$AIC & $\\HL$ & $C$ & $\\EF$ & $p_{\\HL}$ & $p_{\\EF}$ & $p$ Stukel & source", body)
