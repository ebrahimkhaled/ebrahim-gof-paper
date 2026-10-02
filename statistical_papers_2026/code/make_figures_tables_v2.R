## make_figures_tables_v2.R -- figures and tables of the revised Statistical Papers manuscript, from the result files.
## Figures: cairo PDF, Arial, 8-10 pt at final size (text width 131 mm), no titles inside the artwork.
## Fig1 map of the second-order gain with the simulated cells marked; Fig2 robustness screen (stage 1);
## Fig3 second-order prediction against simulation; Fig4 many small groups (power against m); Fig5 beetle read-out.
ROOT <- "."
suppressMessages(library(data.table))
FIG <- file.path(ROOT, "sp", "figures"); TAB <- file.path(ROOT, "sp", "tables"); RES <- file.path(ROOT, "results")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE); dir.create(TAB, showWarnings = FALSE, recursive = TRUE)
fig <- function(f, w_mm, h_mm) {
  grDevices::cairo_pdf(file.path(FIG, f), width = w_mm / 25.4, height = h_mm / 25.4, family = "Arial", pointsize = 8)
  par(cex = 1, mgp = c(1.9, 0.5, 0), tcl = -0.25, las = 1)
}
## colours by family: risk-ordered partition (greens), covariate clusters (blues), record-level directed (browns)
col <- c(EF = "#01665E", HL = "#35978F", PH = "#80CDC1", Tsiatis = "#2166AC", Xie = "#67A9CF",
         Stukel = "#8C510A", Cubic = "#BF812D", GiViTI = "#DFC27D")
lty <- c(EF = 1, HL = 1, PH = 1, Tsiatis = 2, Xie = 2, Stukel = 3, Cubic = 3, GiViTI = 3)
pch <- c(EF = 16, HL = 1, PH = 2, Tsiatis = 15, Xie = 0, Stukel = 17, Cubic = 6, GiViTI = 4)
f2 <- function(v) { s <- sprintf("%+.2f", v); s[s %in% c("+0.00", "-0.00")] <- "0.00"; s }
f3 <- function(v) sprintf("%.3f", v)
mm <- function(s) gsub("(^|[ (])-([0-9])", "\\1$-$\\2", s)          # typographic minus in table cells

## ------------------------------------------------------------------------------------------- Fig. 1 (map)
M  <- fread(file.path(RES, "N2_sign_map.csv"))
S6 <- fread(file.path(RES, "study6_summary.csv"))
pal <- colorRampPalette(c("#8C510A", "#D8B365", "#F5F5F5", "#5AB4AC", "#01665E"))(41); zlim <- c(-0.35, 0.35)
fig("Fig1.pdf", 131, 112)
layout(matrix(c(1:6, 7, 7, 7), 3, 3, byrow = TRUE), heights = c(1, 1, 0.26))
par(mar = c(2.5, 3.0, 1.5, 0.3), oma = c(0, 0.6, 0, 0))
for (sd_ in c("cloglog", "loglog")) for (sp_ in c(0.5, 1, 2)) {
  S <- M[M$side == sd_ & M$s == sp_]
  Z <- rbind(as.matrix(dcast(S, lam ~ mu, value.var = "gain_n1000")[, -1]), 0)   # lambda = 1: the logit itself, gain 0
  mus <- sort(unique(S$mu)); lams <- c(sort(unique(S$lam)), 1)
  image(seq_along(mus), seq_along(lams), t(Z), col = pal, zlim = zlim, axes = FALSE, xlab = "", ylab = "")
  axis(1, seq_along(mus), mus); axis(2, seq_along(lams), c(head(lams, -1), "1"), cex.axis = 0.95)
  for (i in seq_along(mus)) for (j in seq_along(lams))
    text(i, j, f2(Z[j, i]), cex = 1.0, col = if (abs(Z[j, i]) > 0.2) "white" else "black")
  ## simulated cells: outlined, simulated gain (size-adjusted, else raw) printed under the prediction
  for (r in seq_len(nrow(S6))) {
    x <- S6[r]; if (x$side != sd_ || x$s != sp_) next
    i <- match(x$mu, mus); j <- match(x$lam, lams); if (is.na(i) || is.na(j)) next
    rect(i - 0.5, j - 0.5, i + 0.5, j + 0.5, border = "black", lwd = 1.4)
    sim <- if (is.na(x$sim_adj)) x$sim_raw else x$sim_adj
    text(i, j - 0.32, sprintf("(%s)", f2(sim)), cex = 0.9, col = if (abs(Z[j, i]) > 0.2) "white" else "black")
  }
  box()
  mtext(sprintf("(%s) %s, spread %s", letters[(sd_ == "loglog") * 3 + match(sp_, c(0.5, 1, 2))], sd_, sp_),
        side = 3, line = 0.25, cex = 0.95, adj = 0)
  if (sp_ == 0.5) mtext(expression(lambda), side = 2, line = 2.0, cex = 1, las = 1)
  mtext(expression(mu), side = 1, line = 1.5, cex = 1)
}
par(mar = c(1.9, 9, 0.9, 9))
image(seq(zlim[1], zlim[2], length.out = 41), 1, matrix(seq(zlim[1], zlim[2], length.out = 41)), col = pal,
      axes = FALSE, xlab = "", ylab = "")
axis(1, at = seq(-0.3, 0.3, 0.1), labels = f2(seq(-0.3, 0.3, 0.1)))
mtext("gain in power of EF over HL", side = 3, line = 0.1, cex = 0.95)
dev.off()

## ------------------------------------------------------------------------------------------- Fig. 2 (stage 1)
S4 <- fread(file.path(RES, "study4_summary.csv"))
tests <- c("EF", "HL", "PH", "Tsiatis", "Xie", "Stukel", "Cubic", "GiViTI")
lab <- c(EF = "EF", HL = "HL", PH = "Pigeon-Heyse", Tsiatis = "Tsiatis", Xie = "Xie", Stukel = "Stukel",
         Cubic = "cubic LR", GiViTI = "GiViTI")
fig("Fig2.pdf", 131, 62)
par(mfrow = c(1, 2), mar = c(3.0, 3.1, 1.3, 0.4))
for (mu in c(4, -4)) {
  S <- rbind(S4[design == "D1" & k == 0], S4[design == "D1" & mult == mu])[order(k)]
  x <- seq_len(nrow(S))
  plot(NA, xlim = range(x), ylim = c(0, 1), xaxt = "n", xlab = "Corrupted records in 1000",
       ylab = "False-alarm rate at 5%")
  abline(h = 0.10, lty = 3, col = "grey55"); abline(h = 0.05, lty = 2, col = "grey55")
  for (t in rev(tests)) lines(x, S[[t]], type = "b", col = col[t], lty = lty[t], pch = pch[t], cex = 0.8, lwd = 1.2)
  axis(1, x, S$k)
  mtext(if (mu == 4) "(a) covariate multiplied by 4" else "(b) covariate multiplied by -4", side = 3, line = 0.2,
        adj = 0, cex = 0.95)
  if (mu == 4) legend("topleft", lab[tests], col = col[tests], lty = lty[tests], pch = pch[tests], bty = "n",
                      cex = 0.85, ncol = 2, lwd = 1.2)
}
dev.off()

## ------------------------------------------------------------------------------------------- Fig. 3 (second order)
V <- fread(file.path(RES, "verify_N1_second_order.csv"))
se <- sqrt(2 * 0.25 / 2000)
fig("Fig3.pdf", 84, 86)
par(mar = c(3.6, 3.6, 0.5, 0.5), mgp = c(2.3, 0.5, 0))
pd <- c(D1 = 16, D3 = 17, D4 = 1); cd <- c(D1 = col[["EF"]], D3 = col[["Tsiatis"]], D4 = "grey40")
plot(V$second_order, V$simulated, xlim = c(-0.06, 0.20), ylim = c(-0.06, 0.12), type = "n",
     xlab = "Predicted gain (Theorem 3)", ylab = "Simulated gain, size-adjusted")
abline(h = 0, v = 0, col = "grey85"); abline(0, 1, lty = 2, col = "grey40")
segments(V$second_order, V$simulated - 1.96 * se, V$second_order, V$simulated + 1.96 * se, col = adjustcolor("grey50", .45))
points(V$second_order, V$simulated, pch = pd[V$design], col = cd[V$design], cex = 0.85)
legend("bottomright", c("balanced (D1)", "skewed covariate (D3)", "strong discrimination (D4)"), pch = pd, col = cd,
       bty = "n", cex = 0.85)
dev.off()

## ------------------------------------------------------------------------------------------- Fig. 4 (many small groups)
P <- fread(file.path(ROOT, "theory", "sparse_sim", "out", "power_summary.csv"))
P <- P[n == 4000]
fig("Fig4.pdf", 131, 58)
par(mfrow = c(1, 3), mar = c(3.0, 3.1, 1.3, 0.4))
for (sc in c("cloglog", "loglog", "quadneg")) {
  S <- P[scen == sc][order(m)]
  plot(NA, xlim = range(S$m), ylim = c(0, 1), log = "x", xlab = "Records per group, m", ylab = "Power",
       xaxt = "n")
  axis(1, S$m)
  lines(S$m, S$pred_EF, col = col[["EF"]], lwd = 1.3); points(S$m, S$pow_zEF, col = col[["EF"]], pch = 16, cex = 1.25)
  lines(S$m, S$pred_HL, col = col[["Stukel"]], lwd = 1.3, lty = 2); points(S$m, S$pow_zHL, col = col[["Stukel"]], pch = 1, cex = 1.25)
  abline(h = 0.05, lty = 3, col = "grey55")
  A <- sign(S$A[1])
  mtext(sprintf("(%s) %s, A %s 0", letters[match(sc, c("cloglog", "loglog", "quadneg"))],
                c(cloglog = "cloglog", loglog = "log-log", quadneg = "omitted square")[sc], if (A < 0) "<" else ">"),
        side = 3, line = 0.2, adj = 0, cex = 0.95)
  if (sc == "cloglog") legend("topleft", c("EF", "HL"), col = col[c("EF", "Stukel")], pch = c(16, 1), lty = c(1, 2),
                               bty = "n", cex = 0.9)
}
dev.off()

## ------------------------------------------------------------------------------------------- Fig. 5 (beetle)
B <- readRDS(file.path(RES, "beetle.rds"))
fig("Fig5.pdf", 84, 64)
par(mar = c(3.4, 3.1, 0.5, 0.5))
bp <- barplot(rbind(B$z$r, B$z$b * B$z$r), beside = TRUE, col = c("grey78", col[["EF"]]), border = NA,
              names.arg = sprintf("%.2f", B$z$pbar), ylim = c(-2, 2.2), xlab = "", ylab = "Contribution")
mtext("Mean fitted risk of the dose group", side = 1, line = 2.1)
abline(h = 0, col = "grey40")
legend("topright", c(expression(r[g]), expression(c[g] == b[g] * r[g])), fill = c("grey78", col[["EF"]]), border = NA,
       bty = "n", cex = 0.9)
dev.off()

## ------------------------------------------------------------------------------------------- Table 1 (stage 1)
rows <- c()
for (d in c("D1", "D5")) {
  S <- rbind(S4[design == d & k == 0], S4[design == d & k > 0])[order(-mult, k)]
  for (r in seq_len(nrow(S))) {
    x <- S[r]
    err <- if (x$k == 0) "none" else if (x$mult > 0) "$\\times4$" else "$\\times(-4)$"
    cells <- sapply(tests, function(t) { v <- x[[t]]; if (v > 0.10) paste0("\\textbf{", f3(v), "}") else f3(v) })
    rows <- c(rows, sprintf("%s & %s & %d & %s \\\\", if (r == 1) c(D1 = "D1 (two covariates)", D5 = "D5 (six covariates)")[d] else "",
                            err, x$k, paste(cells, collapse = " & ")))
  }
  rows <- c(rows, "\\addlinespace")
}
writeLines(c("\\begin{tabular}{@{}llrrrrrrrrr@{}}", "\\toprule",
             " & & & \\multicolumn{3}{c}{risk-ordered groups} & \\multicolumn{2}{c}{covariate clusters} & \\multicolumn{3}{c}{record-level directed} \\\\",
             "\\cmidrule(lr){4-6}\\cmidrule(lr){7-8}\\cmidrule(l){9-11}",
             "Design & Error & $k$ & EF & HL & P--H & Tsiatis & Xie & Stukel & Cubic & GiViTI \\\\", "\\midrule",
             head(rows, -1), "\\bottomrule", "\\end{tabular}"), file.path(TAB, "tab_screen.tex"))

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
S5 <- fread(file.path(RES, "study5_summary.csv"))[truth != "logit" & n == 1000]
tl <- c(cloglog = "cloglog", loglog = "log-log", probit = "probit", quad = "omitted square", inter = "omitted interaction")
S5[, truth := factor(truth, names(tl))]; setorder(S5, design, truth)
rows <- S5[, {
  best <- c("EF", "HL", "PH")[which.max(c(adj_EF, adj_HL, adj_PH))]
  bf <- function(t, v) if (t == best && abs(adj_EF - adj_HL) > 0.01) paste0("\\textbf{", f3(v), "}") else f3(v)
  sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\", if (truth == "cloglog") design else "", tl[as.character(truth)],
          mm(sprintf("%.2f", A)), bf("EF", adj_EF), bf("HL", adj_HL), bf("PH", adj_PH), f3(raw_EF), f3(raw_HL), f3(adj_Tsiatis),
          f3(adj_Stukel)) }, by = .(design, truth)]$V1
rows <- unlist(lapply(split(rows, rep(1:5, each = 5)), function(r) c(r, "\\addlinespace")))
writeLines(c("\\begin{tabular}{@{}llrrrrrrrr@{}}", "\\toprule",
             " & & & \\multicolumn{3}{c}{size-adjusted} & \\multicolumn{2}{c}{raw} & \\multicolumn{2}{c}{reference (adjusted)} \\\\",
             "\\cmidrule(lr){4-6}\\cmidrule(lr){7-8}\\cmidrule(l){9-10}",
             "Design & Truth & $A$ & EF & HL & P--H & EF & HL & Tsiatis & Stukel \\\\", "\\midrule",
             head(rows, -1), "\\bottomrule", "\\end{tabular}"), file.path(TAB, "tab_power.tex"))

## ------------------------------------------------------------------------------------------- Table 4 (plasmode)
PL0 <- fread(file.path(RES, "plasmode_summary.csv")); S7 <- fread(file.path(RES, "study7_summary.csv"))
N00 <- fread(file.path(RES, "plasmode", "logit_k00.csv.gz"))
adj0 <- function(tr, s) { R <- fread(file.path(RES, "plasmode", sprintf("%s_k00.csv.gz", tr)))
  mean(R[[s]] > quantile(N00[[s]], .95, type = 1, na.rm = TRUE), na.rm = TRUE) }
rowsN <- c(sprintf("correct model & 0 & %s & %s & %s & %s & --- & --- \\\\",
                   f3(PL0[truth == "logit" & k == 0]$EF), f3(PL0[truth == "logit" & k == 0]$HL),
                   "---", f3(PL0[truth == "logit" & k == 0]$Stukel)),
           S7[truth == "logit", sprintf(" & %d & %s & %s & %s & %s & %s & %s \\\\", k, f3(EF), f3(HL), f3(PH), f3(Stukel), f3(Cubic), f3(GiViTI))])
rowsA <- unlist(lapply(c("loglog", "cloglog"), function(tr) c(
  sprintf("%s-type truth & 0 & %s & %s & \\multicolumn{4}{l}{raw: EF %s, HL %s} \\\\", tr, f3(adj0(tr, "EF")), f3(adj0(tr, "HL")),
          f3(PL0[truth == tr & k == 0]$EF), f3(PL0[truth == tr & k == 0]$HL)),
  S7[truth == tr, sprintf(" & %d & %s & %s & \\multicolumn{4}{l}{raw: EF %s, HL %s} \\\\", k, f3(adj_EF), f3(adj_HL), f3(EF), f3(HL))])))
writeLines(c("\\begin{tabular}{@{}lrrrrrrr@{}}", "\\toprule",
             "\\multicolumn{8}{@{}l}{(a) Correct model: raw false-alarm rate at 5\\%} \\\\",
             "Truth & $k$ & EF & HL & P--H & Stukel & Cubic & GiViTI \\\\", "\\midrule", rowsN, "\\midrule",
             "\\multicolumn{8}{@{}l}{(b) Asymmetric truths: size-adjusted power of EF and HL} \\\\",
             "Truth & $k$ & EF & HL & \\multicolumn{4}{l}{} \\\\", "\\midrule", rowsA, "\\bottomrule", "\\end{tabular}"),
           file.path(TAB, "tab_plasmode.tex"))
cat("figures and tables written\n")
