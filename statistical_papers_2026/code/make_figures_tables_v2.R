## make_figures_tables_v2.R -- figures and tables of the revised Statistical Papers manuscript, from the result files.
## Figures: cairo PDF, Arial, 8-10 pt at final size (text width 131 mm), no titles inside the artwork.
## Fig1 map of the second-order gain with the simulated cells marked; Fig2 robustness screen (stage 1);
## Fig3 second-order prediction against simulation; Fig4 many small groups (power against m); Fig5 beetle read-out.
ROOT <- "."
suppressMessages(library(data.table))
FIG <- file.path(ROOT, "sp", "figures"); TAB <- file.path(ROOT, "sp", "tables"); RES <- file.path(ROOT, "results")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE); dir.create(TAB, showWarnings = FALSE, recursive = TRUE)
fig <- function(f, w_mm, h_mm) {
  grDevices::cairo_pdf(file.path(FIG, f), width = w_mm / 25.4, height = h_mm / 25.4, family = "Arial", pointsize = 9)
  par(cex = 1, mgp = c(1.9, 0.5, 0), tcl = -0.25, las = 1)
}
## layout()/mfrow silently shrink text to 66-83%; every multi-panel figure resets cex = 1 after them, and no text
## uses cex below 0.9, so all lettering is at least 8.1 pt at print size (journal rule: 8-12 pt)
## colours by family: risk-ordered partition (greens), covariate clusters (blues), record-level directed (browns)
col <- c(EF = "#01665E", HL = "#35978F", PH = "#80CDC1", Tsiatis = "#2166AC", Xie = "#67A9CF",
         Stukel = "#8C510A", StukelW = "#543005", Cubic = "#BF812D", GiViTI = "#DFC27D")
lty <- c(EF = 1, HL = 1, PH = 1, Tsiatis = 2, Xie = 2, Stukel = 3, StukelW = 4, Cubic = 3, GiViTI = 3)
pch <- c(EF = 16, HL = 1, PH = 2, Tsiatis = 15, Xie = 0, Stukel = 17, StukelW = 18, Cubic = 6, GiViTI = 4)
f2 <- function(v) { s <- sprintf("%+.2f", v); s[s %in% c("+0.00", "-0.00")] <- "0.00"; s }
f3 <- function(v) sprintf("%.3f", v)
mm <- function(s) gsub("(^|[ (])-([0-9])", "\\1$-$\\2", s)          # typographic minus in table cells

## ------------------------------------------------------------------------------------------- Fig. 1 (map)
M  <- fread(file.path(RES, "N2_sign_map.csv"))
S6 <- fread(file.path(RES, "study6_summary.csv"))
pal <- colorRampPalette(c("#8C510A", "#D8B365", "#F5F5F5", "#5AB4AC", "#01665E"))(41); zlim <- c(-0.35, 0.35)
fig("Fig1.pdf", 131, 112)
layout(matrix(c(1:6, 7, 7, 7), 3, 3, byrow = TRUE), heights = c(1, 1, 0.26))
par(mar = c(2.5, 3.0, 1.3, 0.3), oma = c(0, 0.6, 0, 0), cex = 1)
for (sd_ in c("cloglog", "loglog")) for (sp_ in c(0.5, 1, 2)) {
  S <- M[M$side == sd_ & M$s == sp_]
  Z <- rbind(as.matrix(dcast(S, lam ~ mu, value.var = "gain_n1000")[, -1]), 0)   # alpha = 1: the logit itself, gain 0
  mus <- sort(unique(S$mu)); lams <- c(sort(unique(S$lam)), 1)
  image(seq_along(mus), seq_along(lams), t(Z), col = pal, zlim = zlim, axes = FALSE, xlab = "", ylab = "")
  axis(1, seq_along(mus), mus); axis(2, seq_along(lams), c(head(lams, -1), "1"))
  for (i in seq_along(mus)) for (j in seq_along(lams))
    text(i, j + 0.12, sub("^([+-]?)0\\.", "\\1.", f2(Z[j, i])), cex = 0.9, col = if (abs(Z[j, i]) > 0.2) "white" else "black")
  ## simulated cells: outlined, simulated gain printed under the prediction (size-adjusted; at alpha = 1, where
  ## the size-adjusted gain is zero by construction, the raw difference)
  for (r in seq_len(nrow(S6))) {
    x <- S6[r]; if (x$side != sd_ || x$s != sp_) next
    i <- match(x$mu, mus); j <- match(x$lam, lams); if (is.na(i) || is.na(j)) next
    rect(i - 0.5, j - 0.5, i + 0.5, j + 0.5, border = "black", lwd = 1.4)
    sim <- if (x$lam == 1 || is.na(x$sim_adj)) x$sim_raw else x$sim_adj
    text(i, j - 0.27, sprintf("(%s)", sub("^([+-]?)0\\.", "\\1.", f2(sim))), cex = 0.9,
         col = if (abs(Z[j, i]) > 0.2) "white" else "black")
  }
  box()
  mtext(sprintf("(%s)", letters[(sd_ == "loglog") * 3 + match(sp_, c(0.5, 1, 2))]), side = 3, line = 0.2, adj = 0)
  if (sp_ == 0.5) mtext(expression(alpha), side = 2, line = 2.0, las = 1)
  mtext(expression(mu), side = 1, line = 1.5)
}
par(mar = c(1.9, 9, 1.1, 9))
image(seq(zlim[1], zlim[2], length.out = 41), 1, matrix(seq(zlim[1], zlim[2], length.out = 41)), col = pal,
      axes = FALSE, xlab = "", ylab = "")
axis(1, at = seq(-0.3, 0.3, 0.1), labels = f2(seq(-0.3, 0.3, 0.1)))
mtext("Gain in power of EF over HL", side = 3, line = 0.15)
dev.off()

## ------------------------------------------------------------------------------------------- Fig. 2 (stage 1)
S8 <- fread(file.path(RES, "study8_summary.csv"))
tests <- c("EF", "HL", "PH", "Tsiatis", "Xie", "Stukel", "StukelW", "Cubic", "GiViTI")
lab <- c(EF = "EF", HL = "HL", PH = "Pigeon-Heyse", Tsiatis = "Tsiatis", Xie = "Xie", Stukel = "Stukel",
         StukelW = "Stukel-W", Cubic = "cubic LR", GiViTI = "GiViTI")
fig("Fig2.pdf", 131, 112)
par(mfrow = c(2, 2), mar = c(3.0, 3.1, 1.3, 0.4), cex = 1)
for (er in c("x4", "x(-4)", "sign", "y flipped")) {
  S <- rbind(S8[design == "D1" & error == "none"], S8[design == "D1" & error == er])[order(k)]
  x <- seq_len(nrow(S))
  plot(NA, xlim = range(x), ylim = if (er %in% c("sign", "y flipped")) c(0, 0.15) else c(0, 1), xaxt = "n",
       xlab = if (er == "y flipped") "Misclassified responses in 1000"
       else "Corrupted records in 1000", ylab = "False-alarm rate at 5%")
  abline(h = 0.10, lty = 3, col = "grey55"); abline(h = 0.05, lty = 2, col = "grey55")
  for (t in rev(tests)) lines(x, S[[t]], type = "b", col = col[t], lty = lty[t], pch = pch[t], cex = 0.8, lwd = 1.2)
  axis(1, x, S$k)
  mtext(sprintf("(%s)", letters[match(er, c("x4", "x(-4)", "sign", "y flipped"))]), side = 3, line = 0.2, adj = 0)
  if (er == "x4") legend("topleft", lab[tests], col = col[tests], lty = lty[tests], pch = pch[tests], bty = "n",
                         cex = 0.9, ncol = 2, lwd = 1.2)
}
dev.off()

## ------------------------------------------------------------------------------------------- Fig. 3 (second order)
V <- fread(file.path(RES, "verify_N1_second_order.csv"))
se <- sqrt(2 * 0.25 / 2000)
fig("Fig3.pdf", 84, 86)
par(mar = c(3.6, 3.6, 0.5, 0.5), mgp = c(2.3, 0.5, 0))
pd <- c(D1 = 16, D3 = 17, D4 = 1); cd <- c(D1 = col[["EF"]], D3 = col[["Tsiatis"]], D4 = "grey40")
plot(V$second_order, V$simulated, xlim = c(-0.06, 0.20), ylim = c(-0.06, 0.12), type = "n",
     xlab = "Predicted gain (Theorem 2)", ylab = "Simulated gain, size-adjusted")
abline(h = 0, v = 0, col = "grey85"); abline(0, 1, lty = 2, col = "grey40")
segments(V$second_order, V$simulated - 1.96 * se, V$second_order, V$simulated + 1.96 * se, col = adjustcolor("grey50", .45))
points(V$second_order, V$simulated, pch = pd[V$design], col = cd[V$design], cex = 0.85)
legend("bottomright", c("balanced (D1)", "skewed covariate (D3)", "strong discrimination (D4)"), pch = pd, col = cd,
       bty = "n", cex = 0.9)
dev.off()

## ------------------------------------------------------------------------------------------- Fig. 4 (many small groups)
P <- fread(file.path(ROOT, "theory", "sparse_sim", "out", "power_summary.csv"))
P <- P[n == 4000]
fig("Fig4.pdf", 131, 62)
par(mfrow = c(1, 3), mar = c(3.0, 3.1, 1.3, 0.4), cex = 1)
for (sc in c("cloglog", "loglog", "quadneg")) {
  S <- P[scen == sc][order(m)]
  plot(NA, xlim = range(S$m), ylim = c(0, 1), log = "x", xlab = "Records per group, m", ylab = "Power",
       xaxt = "n")
  axis(1, S$m)
  lines(S$m, S$pred_EF, col = col[["EF"]], lwd = 1.3); points(S$m, S$pow_zEF, col = col[["EF"]], pch = 16, cex = 1.25)
  lines(S$m, S$pred_HL, col = col[["Stukel"]], lwd = 1.3, lty = 2); points(S$m, S$pow_zHL, col = col[["Stukel"]], pch = 1, cex = 1.25)
  abline(h = 0.05, lty = 3, col = "grey55")
  A <- sign(S$A[1])
  mtext(sprintf("(%s)", letters[match(sc, c("cloglog", "loglog", "quadneg"))]), side = 3, line = 0.2, adj = 0)
  if (sc == "cloglog") legend("topleft", c("EF", "HL"), col = col[c("EF", "Stukel")], pch = c(16, 1), lty = c(1, 2),
                               bty = "n", cex = 0.9)
}
dev.off()

## ------------------------------------------------------------------------------------------- Fig. 5 (Paul's rule)
S12 <- fread(file.path(RES, "study12_summary.csv"))
fig("Fig5.pdf", 131, 62)
par(mfrow = c(1, 3), mar = c(3.0, 3.1, 1.3, 0.4), cex = 1)
for (pn in list(c("D1", "cloglog"), c("D2", "loglog"), c("D4", "quad"))) {
  S <- S12[design == pn[1] & truth == pn[2]][order(n)]
  plot(NA, xlim = range(S$n), ylim = c(0, 1), log = "x", xaxt = "n", xlab = "Sample size, n", ylab = "Power, size-adjusted")
  axis(1, S$n, c("2000", "5000", "10000", "20000"))
  lines(S$n, S$adj_p_EF_norm_P, col = col[["EF"]], lwd = 1.3); points(S$n, S$adj_p_EF_norm_P, col = col[["EF"]], pch = 16)
  lines(S$n, S$adj_p_HL_chisq_P, col = col[["Stukel"]], lwd = 1.3, lty = 2); points(S$n, S$adj_p_HL_chisq_P, col = col[["Stukel"]], pch = 1)
  abline(h = 0.05, lty = 3, col = "grey55")
  mtext(sprintf("(%s)", letters[match(pn[1], c("D1", "D2", "D4"))]), side = 3, line = 0.2, adj = 0)
  if (pn[1] == "D1") legend("topright", c("EF, normal", "HL, chi-square"), col = col[c("EF", "Stukel")], pch = c(16, 1),
                            lty = c(1, 2), bty = "n", cex = 0.9)
}
dev.off()

## Table 5: size under Paul's rule (paper); Study 12 power and Studies 11, 13 in Online Resource 1
nl <- S12[truth == "logit"][order(design, n)]
rows <- nl[, sprintf("%s & %d & %d & %s & %s & %s & %s \\\\", ifelse(n == 2000, design, ""), n, as.integer(G_paul),
                     f3(raw_p_HL_chisq_P), f3(raw_p_EF_chisq_P), f3(raw_p_HL_norm_P), f3(raw_p_EF_norm_P))]
rows <- unlist(lapply(split(rows, rep(1:3, each = 4)), function(r) c(r, "\\addlinespace")))
writeLines(c("\\begin{tabular}{@{}lrrrrrr@{}}", "\\toprule",
             " & & & \\multicolumn{2}{c}{$\\chi^2_{G-2}$} & \\multicolumn{2}{c}{normal (Theorem 3)} \\\\",
             "\\cmidrule(lr){4-5}\\cmidrule(l){6-7}", "Design & $n$ & $G$ & HL & EF & HL & EF \\\\", "\\midrule",
             head(rows, -1), "\\bottomrule", "\\end{tabular}"), file.path(TAB, "tab_large.tex"))
tl12 <- c(cloglog = "cloglog", loglog = "log-log", quad = "square")
al <- S12[truth != "logit"][, truth := factor(truth, names(tl12))][order(design, truth, n)]
rows <- al[, sprintf("%s & %s & %d & %s & %s & %s & %s & %s & %s & %s & %s \\\\", design, tl12[as.character(truth)], n, mm(sprintf("%.2f", A)),
                     f3(adj_p_EF_norm_P), f3(adj_p_HL_chisq_P), f3(adj_p_HL_norm_P), f3(adj_p_comb_P), f3(adj_p_EF_norm_25),
                     f3(adj_EF10), f3(adj_StukelW))]
writeLines(c("\\begin{tabular}{@{}llrrrrrrrrr@{}}", "\\toprule",
             " & & & & \\multicolumn{4}{c}{Paul's $G$} & 25 per group & $G=10$ & \\\\",
             "\\cmidrule(lr){5-8}",
             "Design & Truth & $n$ & $A$ & EF-N & HL-$\\chi^2$ & HL-N & combined & EF-N & EF & Stukel-W \\\\", "\\midrule",
             rows, "\\bottomrule", "\\end{tabular}"), file.path(ROOT, "sp", "esm", "tables", "tab_S14.tex"))
S11 <- fread(file.path(RES, "study11_summary.csv"))[truth == "shared"][order(design, tau)]
S13 <- fread(file.path(RES, "study13_summary.csv"))[order(k, -mult)]
rowsA <- S11[, sprintf("%s & %.1f & %s & %s & %s & %s & %s & %s & %s \\\\", ifelse(design == "W", "wide", design), tau, f3(adj_p_EF_norm_P),
                       f3(adj_p_HL_norm_P), f3(adj_p_EF_norm_25), f3(adj_EF10), f3(adj_HL10), f3(adj_Stukel), f3(adj_StukelW))]
rowsB <- S13[, sprintf("%s & %d & %s & %s & %s & %s & %s & %s & %s \\\\", ifelse(k == 0, "none", ifelse(mult > 0, "$\\times4$", "$\\times(-4)$")), k,
                       f3(raw_p_EF_norm_P), f3(raw_p_HL_norm_P), f3(raw_p_L_P), f3(raw_p_comb_P), f3(raw_EF10), f3(raw_Stukel), f3(raw_StukelW))]
writeLines(c("\\begin{tabular}{@{}lrrrrrrrr@{}}", "\\toprule",
             "\\multicolumn{9}{@{}l}{(a) Shared deviations within risk bins, $n=4000$: size-adjusted power} \\\\",
             "Design & $\\tau$ & EF-N & HL-N & EF-N, 25 & EF, $G=10$ & HL, $G=10$ & Stukel & Stukel-W \\\\", "\\midrule", rowsA, "\\midrule",
             "\\multicolumn{9}{@{}l}{(b) Gross errors, design D1, $n=5000$, Paul's $G$: raw rate at 5\\%} \\\\",
             "Error & $k$ & EF-N & HL-N & $L$ & combined & EF, $G=10$ & Stukel & Stukel-W \\\\", "\\midrule", rowsB,
             "\\bottomrule", "\\end{tabular}"), file.path(ROOT, "sp", "esm", "tables", "tab_S15.tex"))

## ------------------------------------------------------------------------------------------- Fig. 6 (beetle)
B <- readRDS(file.path(RES, "beetle.rds"))
fig("Fig6.pdf", 84, 70)
par(mar = c(4.6, 3.1, 0.5, 0.5))
bp <- barplot(rbind(B$z$r, B$z$b * B$z$r), beside = TRUE, col = c("grey78", col[["EF"]]), border = NA,
              names.arg = rep("", 8), ylim = c(-2, 2.2), xlab = "", ylab = "Contribution")
mid <- colMeans(bp)                                    # two label rows: log dose, then mean fitted risk
mtext(sprintf("%.2f", B$logdose), side = 1, at = mid, line = 0.3, cex = 0.9)
mtext(sprintf("%.2f", B$z$pbar), side = 1, at = mid, line = 1.3, cex = 0.9)
mtext("log dose (top) and mean fitted risk (bottom)", side = 1, line = 2.6)
abline(h = 0, col = "grey40")
legend("topright", c(expression(r[g]), expression(c[g] == b[g] * r[g])), fill = c("grey78", col[["EF"]]), border = NA,
       bty = "n", cex = 0.9)
dev.off()

## ------------------------------------------------------------------------------------------- Table 1 (stage 1)
## Study 8 screen; D1 in the paper (tab_screen), D5 and the capped plasmode in Online Resource 1 (tab_screen_esm)
elab <- c(none = "none", x4 = "$\\times4$", `x(-4)` = "$\\times(-4)$", sign = "sign", `y flipped` = "$y$ flipped",
          `tbsa x4, capped` = "tbsa $\\times4$, $\\le100$")
screen_rows <- function(d, errs) {
  rows <- c()
  for (er in errs) {
    S <- S8[design == d & error == er][order(k)]
    for (r in seq_len(nrow(S))) {
      x <- S[r]
      cells <- sapply(tests, function(t) { v <- x[[t]]; if (v > 0.10) paste0("\\textbf{", f3(v), "}") else f3(v) })
      rows <- c(rows, sprintf("%s & %d & %s \\\\", if (r == 1) elab[[er]] else "", x$k, paste(cells, collapse = " & ")))
    }
    rows <- c(rows, "\\addlinespace")
  }
  head(rows, -1)
}
screen_head <- c("\\toprule",
  " & & \\multicolumn{3}{c}{risk-ordered groups} & \\multicolumn{2}{c}{covariate clusters} & \\multicolumn{4}{c}{record-level directed} \\\\",
  "\\cmidrule(lr){3-5}\\cmidrule(lr){6-7}\\cmidrule(l){8-11}",
  "Error & $k$ & EF & HL & P--H & Tsiatis & Xie & Stukel & Stukel-W & Cubic & GiViTI \\\\", "\\midrule")
writeLines(c("\\begin{tabular}{@{}lrrrrrrrrrr@{}}", screen_head,
             screen_rows("D1", c("none", "x4", "x(-4)", "sign", "y flipped")), "\\bottomrule", "\\end{tabular}"),
           file.path(TAB, "tab_screen.tex"))
writeLines(c("\\begin{tabular}{@{}lrrrrrrrrrr@{}}", screen_head,
             "\\multicolumn{11}{@{}l}{(a) Design D5, six covariates} \\\\",
             screen_rows("D5", c("none", "x4", "x(-4)", "sign", "y flipped")), "\\midrule",
             "\\multicolumn{11}{@{}l}{(b) Burn plasmode, logistic truth} \\\\",
             screen_rows("burn", c("tbsa x4, capped")), "\\bottomrule", "\\end{tabular}"),
           file.path(ROOT, "sp", "esm", "tables", "tab_S9.tex"))

## ------------------------------------------------------------------------------------------- Table 2 (size)
T1 <- fread(file.path(RES, "sparse_threshold.csv"))[order(design, n)]
dl <- c(D1 = "D1 balanced", D2 = "D2 rare events", D3 = "D3 skewed covariate", D4 = "D4 strong discrimination", D5 = "D5 six covariates")
rows <- T1[, sprintf("%s & %d & %.2f & %s & %s \\\\", ifelse(n == 100, dl[design], ""), n, min_expected,
                     ifelse(size_HL < 0.0354 | size_HL > 0.0646, paste0("\\textbf{", f3(size_HL), "}"), f3(size_HL)),
                     ifelse(size_EF < 0.0354 | size_EF > 0.0646, paste0("\\textbf{", f3(size_EF), "}"), f3(size_EF)))]
rows <- unlist(lapply(split(rows, rep(1:5, each = 5)), function(r) c(r, "\\addlinespace")))
writeLines(c("\\begin{tabular}{@{}lrrrr@{}}", "\\toprule", "Design & $n$ & $e_{\\min}$ & HL & EF \\\\", "\\midrule",
             head(rows, -1), "\\bottomrule", "\\end{tabular}"), file.path(TAB, "tab_size.tex"))

## ------------------------------------------------------------------------------------------- Table 3 (stage 2)
## Study 9 (event rate held fixed) at n = 1000 in the paper; all n in Online Resource 1; Study 5 (event rate free) there too
tl <- c(cloglog = "cloglog", loglog = "log-log", probit = "probit", quad = "square", inter = "interaction")
power_rows <- function(S, with_n) {
  S <- copy(S)[truth != "logit"]; S[, truth := factor(truth, names(tl))]; setorder(S, design, truth, n)
  rows <- S[, {
    best <- c("EF", "HL", "PH")[which.max(c(adj_EF, adj_HL, adj_PH))]
    bf <- function(t, v) if (t == best && round(abs(adj_EF - adj_HL), 3) > 0.010) paste0("\\textbf{", f3(v), "}") else f3(v)
    lead <- if (with_n) sprintf("%s & %s & %d", if (truth == "cloglog" && n == min(S$n)) design else "",
                                 if (n == min(S$n)) tl[as.character(truth)] else "", n)
            else sprintf("%s & %s", if (truth == "cloglog") design else "", tl[as.character(truth)])
    sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\", lead, mm(sprintf("%.2f", A)),
            bf("EF", adj_EF), bf("HL", adj_HL), bf("PH", adj_PH),
            if (!is.null(S$se_gain)) mm(sprintf("%.3f (%.3f)", adj_EF - adj_HL, se_gain)) else mm(sprintf("%.3f", adj_EF - adj_HL)),
            f3(raw_EF), f3(raw_HL), f3(adj_Tsiatis), f3(if (!is.null(S$adj_StukelW)) adj_StukelW else adj_Stukel)) }, by = .(design, truth, n)]$V1
  per <- if (with_n) 15 else 5
  head(unlist(lapply(split(rows, rep(seq_len(length(rows) / per), each = per)), function(r) c(r, "\\addlinespace"))), -1)
}
power_head <- function(with_n, last) c("\\toprule",
  sprintf("%s & & \\multicolumn{3}{c}{size-adjusted} & & \\multicolumn{2}{c}{raw} & \\multicolumn{2}{c}{reference (adjusted)} \\\\",
          if (with_n) " & &" else " &"),
  if (with_n) "\\cmidrule(lr){5-7}\\cmidrule(lr){9-10}\\cmidrule(l){11-12}" else "\\cmidrule(lr){4-6}\\cmidrule(lr){8-9}\\cmidrule(l){10-11}",
  sprintf("Design & Truth%s & $A$ & EF & HL & P--H & EF $-$ HL (SE) & EF & HL & Tsiatis & %s \\\\", if (with_n) " & $n$" else "", last),
  "\\midrule")
S9 <- fread(file.path(RES, "study9_summary.csv"))
writeLines(c("\\begin{tabular}{@{}llrrrrrrrrr@{}}", power_head(FALSE, "Stukel-W"), power_rows(S9[n == 1000], FALSE),
             "\\bottomrule", "\\end{tabular}"), file.path(TAB, "tab_power.tex"))
writeLines(c("\\begin{tabular}{@{}lllrrrrrrrrr@{}}", power_head(TRUE, "Stukel-W"), power_rows(S9, TRUE),
             "\\bottomrule", "\\end{tabular}"), file.path(ROOT, "sp", "esm", "tables", "tab_S10.tex"))
S5 <- fread(file.path(RES, "study5_summary.csv"))
writeLines(c("\\begin{tabular}{@{}lllrrrrrrrrr@{}}", power_head(TRUE, "Stukel"), power_rows(S5, TRUE),
             "\\bottomrule", "\\end{tabular}"), file.path(ROOT, "sp", "esm", "tables", "tab_S11.tex"))

## map cells (Study 6): predicted and simulated power of each test, Online Resource 1
M6 <- merge(fread(file.path(RES, "study6_summary.csv")), M[, .(side, lam, mu, s, power_HL)], by = c("side", "lam", "mu", "s"), all.x = TRUE)
setorder(M6, side, -lam, mu)
rows <- M6[, sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\", side, sub("\\.?0+$", "", sprintf("%.2f", lam)), mm(sprintf("%d", as.integer(mu))), s,
                     mm(sprintf("%.3f", predicted)), ifelse(is.na(power_HL), "---", f3(power_HL + predicted)), ifelse(is.na(power_HL), "---", f3(power_HL)),
                     ifelse(is.na(adj_EF), f3(raw_EF), f3(adj_EF)), ifelse(is.na(adj_HL), f3(raw_HL), f3(adj_HL)),
                     mm(sprintf("%.3f", ifelse(lam == 1 | is.na(sim_adj), sim_raw, sim_adj))))]
writeLines(c("\\begin{tabular}{@{}lrrrrrrrrr@{}}", "\\toprule",
             " & & & & \\multicolumn{3}{c}{predicted (Theorem 2)} & \\multicolumn{3}{c}{simulated} \\\\",
             "\\cmidrule(lr){5-7}\\cmidrule(l){8-10}",
             "Side & $\\alpha$ & $\\mu$ & $s$ & gain & EF & HL & EF & HL & gain \\\\", "\\midrule",
             rows, "\\bottomrule", "\\end{tabular}"), file.path(ROOT, "sp", "esm", "tables", "tab_S12.tex"))

## the claims ledger, Online Resource 1
LG <- fread(file.path(RES, "claims_ledger.csv"))
dn <- c(Spec = "Specification", Add1 = "Addendum 1", Add2 = "Addendum 2", Add3 = "Addendum 3", Add4 = "Addendum 4")
esc <- function(s) gsub("%", "\\\\%", gsub("_", "\\\\_", s))
rows <- LG[, sprintf("%s & %s & %s & %s \\\\", ifelse(c(TRUE, declaration[-1] != declaration[-.N]), dn[declaration], ""),
                     esc(text), esc(number), ifelse(verdict == "held", "held", "\\textbf{failed}"))]
writeLines(c("\\begin{tabular}{@{}lp{6.2cm}p{5.6cm}l@{}}", "\\toprule", "Declaration & Claim & Result & Verdict \\\\", "\\midrule",
             rows, "\\bottomrule", "\\end{tabular}"), file.path(ROOT, "sp", "esm", "tables", "tab_S13.tex"))

## ------------------------------------------------------------------------------------------- Table 4 (plasmode)
PL0 <- fread(file.path(RES, "plasmode_summary.csv")); S7 <- fread(file.path(RES, "study7_summary.csv"))
N00 <- fread(file.path(RES, "plasmode", "logit_k00.csv.gz"))
adj0 <- function(tr, s) { R <- fread(file.path(RES, "plasmode", sprintf("%s_k00.csv.gz", tr)))
  mean(R[[s]] > quantile(N00[[s]], .95, type = 1, na.rm = TRUE), na.rm = TRUE) }
## (a) one run (Study 8): the correct model with k patients' tbsa x4, capped at 100; (b) k = 0 from the Addendum 1 run,
## k >= 1 from Study 7 (uncapped tbsa x4 under the alternative)
B8 <- S8[design == "burn"][order(k)]
rowsN <- B8[, sprintf("%d & %s & %s & %s & %s & %s & %s & %s \\\\", k, f3(EF), f3(HL), f3(PH), f3(Stukel), f3(StukelW), f3(Cubic), f3(GiViTI))]
rowsA <- unlist(lapply(c("loglog", "cloglog"), function(tr) c(
  sprintf("%s-type & 0 & %s & %s & %s & %s & Add.~1 \\\\", tr, f3(adj0(tr, "EF")), f3(adj0(tr, "HL")),
          f3(PL0[truth == tr & k == 0]$EF), f3(PL0[truth == tr & k == 0]$HL)),
  S7[truth == tr, sprintf(" & %d & %s & %s & %s & %s & Study~7 \\\\", k, f3(adj_EF), f3(adj_HL), f3(EF), f3(HL))])))
writeLines(c("\\begin{tabular}{@{}lrrrrrrr@{}}", "\\toprule",
             "\\multicolumn{8}{@{}l}{(a) Correct model; $k$ patients with tbsa multiplied by 4, capped at 100: raw rate at 5\\%} \\\\",
             "$k$ & EF & HL & P--H & Stukel & Stukel-W & Cubic & GiViTI \\\\", "\\midrule", rowsN, "\\midrule",
             "\\multicolumn{8}{@{}l}{(b) Asymmetric truths, $k$ patients with tbsa multiplied by 4: power of EF and HL} \\\\",
             "Truth & $k$ & \\multicolumn{2}{c}{size-adjusted} & \\multicolumn{2}{c}{raw} & Run & \\\\",
             " & & EF & HL & EF & HL & & \\\\", "\\midrule", rowsA, "\\bottomrule", "\\end{tabular}"),
           file.path(TAB, "tab_plasmode.tex"))
## computing time (Study 10), Online Resource 1
TT <- fread(file.path(RES, "study10_timing.csv"))
fsec <- function(x) if (x < 0.1) sprintf("%.4f", x) else if (x < 10) sprintf("%.2f", x) else sprintf("%.0f", x)
fx <- function(x) if (x < 10) sprintf("%.1f", x) else format(round(x), big.mark = ",")
TT[, test := sub("le Cessie-van Houwelingen \\(O\\(n\\^3\\) form\\)", "le Cessie--van Houwelingen, original $O(n^3)$ form", test)]
rows <- TT[, sprintf("%s & %d & %s & %s \\\\", test, n, fsec(sec), fx(times_EF)), by = seq_len(nrow(TT))]$V1
writeLines(c("\\begin{tabular}{@{}lrrr@{}}", "\\toprule", "Test & $n$ & seconds & times $\\EF$ \\\\", "\\midrule",
             rows, "\\bottomrule", "\\end{tabular}"), file.path(ROOT, "sp", "esm", "tables", "tab_S8.tex"))
cat("figures and tables written\n")
