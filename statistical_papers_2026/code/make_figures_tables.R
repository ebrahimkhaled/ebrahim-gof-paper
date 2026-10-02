## make_figures_tables.R -- every figure and table of the Statistical Papers manuscript, from the result files.
## Figures: base R, Arial, 8-12 pt at final size, no titles inside the artwork, vector PDF named Fig1..Fig4.
## Tables: LaTeX fragments written into sp/tables/.
ROOT <- "."
source(file.path(ROOT, "code", "ef_exact.R"))
suppressMessages(library(data.table))
FIG <- file.path(ROOT, "sp", "figures"); TAB <- file.path(ROOT, "sp", "tables")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE); dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
RES <- file.path(ROOT, "results")
col_EF <- "#D55E00"; col_HL <- "#0072B2"; col_ST <- "#555555"         # colour-blind safe (Okabe-Ito)
pdf_open <- function(f, w_mm, h_mm) {
  grDevices::cairo_pdf(file.path(FIG, f), width = w_mm / 25.4, height = h_mm / 25.4, family = "Arial", pointsize = 9)
  par(cex = 1, mgp = c(2, 0.6, 0), tcl = -0.3, las = 1)
}
## ---------------------------------------------------------------------------------------------- Fig. 1
## where EF beats HL: second-order gain at n = 1000 over link asymmetry (lambda) and location (mu), three spreads
M <- fread(file.path(RES, "N2_sign_map.csv"))
## brown (loss) -> white -> deep green (gain): the BrBG diverging scheme, readable with red-green colour blindness
pal <- colorRampPalette(c("#8C510A", "#D8B365", "#F5F5F5", "#5AB4AC", "#01665E"))(41)
zlim <- c(-0.35, 0.35)
pdf_open("Fig1.pdf", 131, 105)
layout(matrix(c(1:6, 7, 7, 7), 3, 3, byrow = TRUE), heights = c(1, 1, 0.28))
par(mar = c(2.4, 2.6, 1.5, 0.3), oma = c(0.4, 0, 0, 0))
for (sd_ in c("cloglog", "loglog")) for (sp_ in c(0.5, 1, 2)) {
  S <- M[M$side == sd_ & M$s == sp_]
  Z <- as.matrix(dcast(S, lam ~ mu, value.var = "gain_n1000")[, -1])
  mus <- sort(unique(S$mu)); lams <- sort(unique(S$lam))
  image(seq_along(mus), seq_along(lams), t(Z), col = pal, zlim = zlim, axes = FALSE, xlab = "", ylab = "")
  axis(1, seq_along(mus), mus, cex.axis = 0.9); axis(2, seq_along(lams), lams, cex.axis = 0.9)
  for (i in seq_along(mus)) for (j in seq_along(lams))
    text(i, j, sprintf("%+.2f", Z[j, i]), cex = 0.85, col = if (abs(Z[j, i]) > 0.2) "white" else "black")
  box()
  mtext(sprintf("%s, spread %s", sd_, sp_), side = 3, line = 0.3, cex = 0.9, adj = 0)
  if (sp_ == 0.5) mtext(expression(lambda), side = 2, line = 2.1, cex = 0.85, las = 0)
  mtext(expression(mu), side = 1, line = 1.6, cex = 0.85)
}
par(mar = c(1.6, 12, 0.6, 12))
image(seq(zlim[1], zlim[2], length.out = 41), 1, matrix(seq(zlim[1], zlim[2], length.out = 41)), col = pal,
      axes = FALSE, xlab = "", ylab = "")
axis(1, at = seq(-0.3, 0.3, 0.1), labels = sprintf("%+.1f", seq(-0.3, 0.3, 0.1)), cex.axis = 0.8)
mtext("power of EF minus power of HL", side = 3, line = 0.1, cex = 0.75)
dev.off()

## ---------------------------------------------------------------------------------------------- Fig. 2
## Theorem 3 against simulation: predicted second-order gain vs simulated size-adjusted gain, 60 cells
V <- fread(file.path(RES, "verify_N1_second_order.csv"))
se <- sqrt(2 * 0.25 / 2000)                       # conservative MC SE of a difference of two powers
pdf_open("Fig2.pdf", 84, 84)
par(mar = c(3.2, 4.0, 0.6, 0.6), mgp = c(2.6, 0.6, 0))
pch_d <- c(D1 = 16, D3 = 17, D4 = 1)
lim <- c(-0.06, 0.20)
plot(V$second_order, V$simulated, xlim = lim, ylim = c(-0.06, 0.11), pch = pch_d[V$design],
     col = ifelse(V$design == "D4", col_ST, col_EF), xlab = "Predicted gain (Theorem 3)",
     ylab = "Simulated gain, size-adjusted", cex = 0.8)
abline(0, 1, lty = 2, col = "grey40"); abline(h = 0, v = 0, col = "grey80")
arrows(V$second_order, V$simulated - 1.96 * se, V$second_order, V$simulated + 1.96 * se, length = 0,
       col = adjustcolor("grey50", 0.5))
legend("topleft", c("balanced risks (D1)", "skewed covariate (D3)", "strong discrimination (D4)"),
       pch = pch_d, col = c(col_EF, col_EF, col_ST), bty = "n", cex = 0.8)
dev.off()

## ---------------------------------------------------------------------------------------------- Fig. 3
## corrupted records: false-alarm rate against the number of corrupted records, both error types
st3 <- list.files(file.path(RES, "study3"), full.names = TRUE)
R3 <- rbindlist(lapply(st3, function(f) {
  x <- fread(f); m <- regmatches(basename(f), regexec("_k(\\d+)_m(neg)?(\\d+)", basename(f)))[[1]]
  data.table(k = as.integer(m[2]), mult = if (m[3] == "neg") -as.numeric(m[4]) else as.numeric(m[4]),
             HL = mean(x$p_HL_chisq < .05), EF = mean(x$p_EF_chisq < .05), Stukel = mean(x$p_Stukel < .05, na.rm = TRUE))
}))
base <- R3[k == 0]
pdf_open("Fig3.pdf", 131, 62)
par(mfrow = c(1, 2), mar = c(3.2, 3.3, 1.4, 0.4))
for (mu in c(4, -4)) {
  S <- rbind(base[, .(k, HL, EF, Stukel)], R3[mult == mu, .(k, HL, EF, Stukel)])[order(k)]
  x <- seq_along(S$k)
  plot(x, S$Stukel, type = "b", pch = 15, col = col_ST, ylim = c(0, 1), xaxt = "n",
       xlab = "Corrupted records in 1000", ylab = "False-alarm rate at 5%", cex = 0.8)
  lines(x, S$HL, type = "b", pch = 16, col = col_HL, cex = 0.8)
  lines(x, S$EF, type = "b", pch = 17, col = col_EF, cex = 0.8)
  axis(1, x, S$k); abline(h = c(0.05, 0.10), lty = c(2, 3), col = "grey50")
  mtext(if (mu == 4) "(a) covariate times 4" else "(b) covariate times -4", side = 3, line = 0.2,
        adj = 0, cex = 0.8)
  if (mu == 4) legend("topleft", c("Stukel", "HL", "EF"), pch = c(15, 16, 17), col = c(col_ST, col_HL, col_EF),
                      lty = 1, bty = "n", cex = 0.8)
}
dev.off()

## ---------------------------------------------------------------------------------------------- Fig. 4
## the beetle data, one group per dose: standardized residuals r_g and directional contributions c_g
data("fbeetle", package = "VGAM"); d <- fbeetle
i <- rep(seq_len(nrow(d)), d$n)
y <- unlist(lapply(seq_len(nrow(d)), function(k) rep(c(1, 0), c(d$dead[k], d$n[k] - d$dead[k]))))
x <- d$logdose[i]; fl <- glm(y ~ x, binomial); z <- ef_parts(y, fitted(fl), model.matrix(fl), 8)
pdf_open("Fig4.pdf", 84, 66)
par(mar = c(3.2, 3.4, 0.6, 0.6))
bp <- barplot(rbind(z$r, z$b * z$r), beside = TRUE, col = c("grey75", col_EF), border = NA,
              names.arg = sprintf("%.2f", z$pbar), ylim = c(-2, 2.6), cex.names = 0.9,
              xlab = "Mean fitted risk of the dose group", ylab = "Contribution")
abline(h = 0, col = "grey40")
legend("topright", c(expression(r[g]), expression(c[g] == b[g] * r[g])), fill = c("grey75", col_EF), border = NA,
       bty = "n", cex = 0.8)
dev.off()

## ---------------------------------------------------------------------------------------------- tables
f3 <- function(v) sprintf("%.3f", v)
## Table 1: size at 5% (Study 1) with the smallest expected count
T1 <- fread(file.path(RES, "sparse_threshold.csv"))
S1 <- rbindlist(lapply(list.files(file.path(RES, "study1"), full.names = TRUE), function(f) {
  x <- fread(f); m <- regmatches(basename(f), regexec("^(D\\d)_logit_n(\\d+)", basename(f)))[[1]]
  data.table(design = m[2], n = as.integer(m[3]), Stukel = mean(x$p_Stukel < .05, na.rm = TRUE))
}))
T1 <- merge(T1, S1, by = c("design", "n"))[order(design, n)]
lab <- c(D1 = "balanced", D2 = "rare events", D3 = "skewed covariate", D4 = "strong discrimination", D5 = "six covariates")
rows <- T1[, sprintf("%s & %d & %.2f & %s & %s & %s \\\\", ifelse(n == 100, lab[design], ""), n, min_expected,
                     f3(size_HL), ifelse(size_EF < 0.0354 | size_EF > 0.0646, paste0("\\textbf{", f3(size_EF), "}"), f3(size_EF)),
                     f3(Stukel))]
rows <- unlist(lapply(split(rows, rep(1:5, each = 5)), function(r) c(r, "\\addlinespace")))
writeLines(c("\\begin{tabular}{@{}lrrrrr@{}}", "\\toprule",
             "Design & $n$ & $e_{\\min}$ & HL & EF & Stukel \\\\", "\\midrule", head(rows, -1), "\\bottomrule", "\\end{tabular}"),
           file.path(TAB, "tab_size.tex"))

## Table 2: power at n = 1000 (size-adjusted) with the alignment and the second-order prediction
S2 <- fread(file.path(RES, "study2_summary.csv")); V2 <- fread(file.path(RES, "verify_N1_second_order.csv"))
T2 <- merge(S2[n == 1000], V2[n == 1000, .(design, truth, second_order)], by = c("design", "truth"))
tl <- c(cloglog = "cloglog", loglog = "log-log", probit = "probit", quad = "square", inter = "interaction")
T2[, truth := factor(truth, names(tl))]; T2 <- T2[order(design, truth)]
rows2 <- T2[, sprintf("%s & %s & %+.2f & %+.3f & %+.3f & %s & %s & %s \\\\", ifelse(truth == "cloglog", c(D1 = "D1 balanced", D3 = "D3 skewed", D4 = "D4 strong")[design], ""),
                      tl[as.character(truth)], A, second_order, adj_EFchi - adj_HLchi, f3(adj_HLchi), f3(adj_EFchi), f3(adj_Stk))]
rows2 <- gsub("-0\\.000", "0.000", rows2, fixed = FALSE)          # a rounded zero carries no sign
rows2 <- gsub("([ &])-([0-9])", "\\1$-$\\2", rows2)                  # a typographic minus
rows2 <- unlist(lapply(split(rows2, rep(1:3, each = 5)), function(r) c(r, "\\addlinespace")))
writeLines(c("\\begin{tabular}{@{}llrrrrrr@{}}", "\\toprule",
             "Design & Truth & $A$ & Predicted & Simulated & HL & EF & Stukel \\\\",
             " & & & \\multicolumn{2}{c}{gain of EF} & \\multicolumn{3}{c}{size-adjusted power} \\\\",
             "\\midrule", head(rows2, -1), "\\bottomrule", "\\end{tabular}"), file.path(TAB, "tab_power.tex"))

## Table 3: the plasmode study
PL <- fread(file.path(RES, "plasmode_summary.csv")); N0 <- fread(file.path(RES, "plasmode", "logit_k00.csv.gz"))
adj <- function(tr, s) { R <- fread(file.path(RES, "plasmode", sprintf("%s_k00.csv.gz", tr)))
  mean(R[[s]] > quantile(N0[[s]], .95, type = 1, na.rm = TRUE), na.rm = TRUE) }
rows3 <- c(
  sprintf("Correct model & raw & %s & %s & %s \\\\", f3(PL[truth == "logit" & k == 0]$HL), f3(PL[truth == "logit" & k == 0]$EF),
          f3(PL[truth == "logit" & k == 0]$Stukel)),
  sapply(c(1, 2, 5, 10), function(kk) sprintf("\\quad %d corrupted record%s & raw & %s & %s & %s \\\\", kk, ifelse(kk > 1, "s", ""),
          f3(PL[truth == "logit" & k == kk]$HL), f3(PL[truth == "logit" & k == kk]$EF), f3(PL[truth == "logit" & k == kk]$Stukel))),
  "\\addlinespace",
  sprintf("loglog-type truth ($A=%.1f$) & size-adjusted & %s & %s & %s \\\\", -6.307, f3(adj("loglog", "HL")), f3(adj("loglog", "EF")), f3(adj("loglog", "Stukel"))),
  sprintf("cloglog-type truth ($A=+%.1f$) & size-adjusted & %s & %s & %s \\\\", 8.910, f3(adj("cloglog", "HL")), f3(adj("cloglog", "EF")), f3(adj("cloglog", "Stukel"))))
writeLines(c("\\begin{tabular}{@{}llrrr@{}}", "\\toprule", "Setting & Rate & HL & EF & Stukel \\\\", "\\midrule", rows3,
             "\\bottomrule", "\\end{tabular}"), file.path(TAB, "tab_plasmode.tex"))
cat("figures and tables written\n")
