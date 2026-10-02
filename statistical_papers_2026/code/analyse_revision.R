## analyse_revision.R -- Studies 4-7 (Addendum 2, hash c9de6124): every summary the revised paper uses, and the
## verdict on each claim written before the run. Missing p-values (a test that declined) count as no rejection.
ROOT <- "."
suppressMessages(library(data.table))
source(file.path(ROOT, "code", "ef_exact.R"))
RES <- file.path(ROOT, "results")
tests <- c("EF", "HL", "PH", "Tsiatis", "Xie", "Stukel", "Cubic", "GiViTI")
rej  <- function(p) mean(!is.na(p) & p < .05)
adjr <- function(p1, p0) { cr <- quantile(p0[!is.na(p0)], .05, type = 1); mean(!is.na(p1) & p1 <= cr) }
rd <- function(st, id) fread(file.path(RES, paste0("study", st), paste0(id, ".csv.gz")))

## ------------------------------------------------------------------ Study 4: the robustness screen
S4 <- rbindlist(lapply(list.files(file.path(RES, "study4"), "gz$"), function(f) {
  x <- rd(4, sub(".csv.gz", "", f)); m <- regmatches(f, regexec("^(D\\d)_logit_n01000_k(\\d+)_m(neg)?(\\d+)", f))[[1]]
  data.table(design = m[2], k = as.integer(m[3]), mult = if (m[4] == "neg") -4 else if (m[3] == "00") 1 else 4,
             t(sapply(tests, function(t) rej(x[[t]]))))
}))
setorder(S4, design, mult, k); fwrite(S4, file.path(RES, "study4_summary.csv"))
cat("== Study 4 (raw false-alarm rates)\n"); print(S4, digits = 3)
D1 <- S4[design == "D1"]
cat("\nclaims: (a) Stukel, Cubic > 0.10 at k=1 x4 in D1:", D1[k == 1 & mult == 4, Stukel > .10 & Cubic > .10],
    "\n        (b) GiViTI > 0.10 by k=2 x4 in D1:", D1[k %in% 1:2 & mult == 4, any(GiViTI > .10)],
    "\n        (c) EF, HL, PH <= 0.10 up to k=10 x4, both designs:", S4[mult == 4, all(EF <= .10 & HL <= .10 & PH <= .10)],
    "\n        (d) Tsiatis and Xie > 0.10 by k=10 x4 in D1:", D1[k == 10 & mult == 4, Tsiatis > .10 & Xie > .10], "\n")

## ------------------------------------------------------------------ alignment A for every design x truth
src <- readLines(file.path(ROOT, "code", "run_ef_studies.R"))
eval(parse(text = src[grep("^design_data <- function", src):(grep("^## one replicate", src) - 1)]))
Afile <- file.path(RES, "alignment_all_designs.csv")
if (!file.exists(Afile)) {
  AA <- rbindlist(lapply(c("D1", "D2", "D3", "D4", "D5"), function(design) rbindlist(lapply(c("cloglog", "loglog", "probit", "quad", "inter"), function(truth) {
    set.seed(777); d <- design_data(design, 1e6); pt <- true_p(truth, d); X <- cbind(1, d$X)
    ps <- suppressWarnings(glm.fit(X, pt, family = quasibinomial()))$fitted.values
    g <- ef_groups(ps, 10); ng <- tabulate(g, 10)
    dg <- as.numeric(rowsum(pt - ps, g, reorder = TRUE)) / ng; pb <- as.numeric(rowsum(ps, g, reorder = TRUE)) / ng
    data.table(design, truth, A = sum((1 - 2 * pb) * dg / (pb * (1 - pb))), event_rate = mean(pt))
  }))))
  fwrite(AA, Afile)
}
AA <- fread(Afile)

## ------------------------------------------------------------------ Study 5: power within the partition family
S5 <- rbindlist(lapply(c("D1", "D2", "D3", "D4", "D5"), function(design) rbindlist(lapply(c(500L, 1000L, 2000L), function(n) {
  N0 <- rd(5, sprintf("%s_logit_n%05d_k00_m1", design, n))
  rbindlist(lapply(c("logit", "cloglog", "loglog", "probit", "quad", "inter"), function(truth) {
    x <- rd(5, sprintf("%s_%s_n%05d_k00_m1", design, truth, n))
    row <- data.table(design, truth, n, meanC = mean(x$C))
    for (t in c(tests, "EF5", "HL5", "EF20", "HL20")) {
      set(row, j = paste0("raw_", t), value = rej(x[[t]]))
      set(row, j = paste0("adj_", t), value = adjr(x[[t]], N0[[t]]))
    }
    row }))
}))))
S5 <- merge(S5, AA, by = c("design", "truth"), all.x = TRUE)
S5[, `:=`(gain_adj = adj_EF - adj_HL, gain_raw = raw_EF - raw_HL, gain5 = adj_EF5 - adj_HL5, gain20 = adj_EF20 - adj_HL20)]
fwrite(S5, file.path(RES, "study5_summary.csv"))
A5 <- S5[truth != "logit"]
cat("\n== Study 5 (size-adjusted power, n = 1000)\n")
print(A5[n == 1000, .(design, truth, A = round(A, 2), EF = adj_EF, HL = adj_HL, PH = adj_PH, gain = round(gain_adj, 3),
                      raw_EF, raw_HL, Ts = adj_Tsiatis, Xie = adj_Xie, Stk = adj_Stukel)], digits = 3)
## claim (e): among EF, HL, PH, EF best where A < 0, HL or PH best where A > 0 (cells with a clear difference)
A5[, best := c("EF", "HL", "PH")[max.col(as.matrix(.SD), ties.method = "first")], .SDcols = c("adj_EF", "adj_HL", "adj_PH")]
cl <- A5[abs(gain_adj) > 0.01]
cat("\nclaim (e): cells with |gain| > 0.01:", nrow(cl), "; EF best where A<0:", cl[A < 0, sum(best == "EF")], "of", cl[A < 0, .N],
    "; HL/PH best where A>0:", cl[A > 0, sum(best != "EF")], "of", cl[A > 0, .N], "\n")
## claim (f): sign of mean(C) equals sign of A where |A| > 0.2
cf <- A5[abs(A) > 0.2]
cat("claim (f): sign(mean C) = sign(A) in", cf[, sum(sign(meanC) == sign(A))], "of", nrow(cf), "cells with |A| > 0.2 (",
    round(100 * cf[, mean(sign(meanC) == sign(A))], 1), "% ; claim needs >= 90%)\n")
cat("sign of size-adjusted gain = -sign(A) where |gain| > 0.01:", cl[, sum(sign(gain_adj) == -sign(A))], "of", nrow(cl), "\n")
cat("G effect (mean size-adjusted gain EF - HL over alternatives): G=5", round(mean(A5$gain5), 4), " G=10", round(mean(A5$gain_adj), 4),
    " G=20", round(mean(A5$gain20), 4), "\n")
cat("PH vs HL, size-adjusted: mean difference", round(mean(A5$adj_PH - A5$adj_HL), 4), "; raw:", round(mean(A5$raw_PH - A5$raw_HL), 4), "\n")
cat("sizes (logit rows): EF", paste(range(S5[truth == "logit"]$raw_EF), collapse = "-"), " HL", paste(range(S5[truth == "logit"]$raw_HL), collapse = "-"),
    " PH", paste(range(S5[truth == "logit"]$raw_PH), collapse = "-"), "\n")

## ------------------------------------------------------------------ Study 6: the map cells
M <- fread(file.path(RES, "N2_sign_map.csv"))
null_of <- function(mu, s) { f <- file.path(RES, "study6", sprintf("map_cloglog_l1.00_mu%+d_s%.1f.csv.gz", as.integer(mu), s)); if (file.exists(f)) fread(f) else NULL }
S6 <- rbindlist(lapply(list.files(file.path(RES, "study6"), "gz$"), function(f) {
  m <- regmatches(f, regexec("^map_(\\w+)_l([0-9.]+)_mu([+-]\\d+)_s([0-9.]+)", f))[[1]]
  side <- m[2]; lam <- as.numeric(m[3]); mu <- as.numeric(m[4]); s <- as.numeric(m[5])
  x <- fread(file.path(RES, "study6", f)); N0 <- null_of(mu, s)
  pred <- if (lam < 1) M[M$side == side & M$lam == lam & M$mu == mu & M$s == s]$gain_n1000 else 0
  data.table(side, lam, mu, s, pred = pred, raw_EF = rej(x$EF), raw_HL = rej(x$HL),
             adj_EF = if (is.null(N0)) NA_real_ else adjr(x$EF, N0$EF), adj_HL = if (is.null(N0)) NA_real_ else adjr(x$HL, N0$HL))
}))
S6[, `:=`(sim_raw = raw_EF - raw_HL, sim_adj = adj_EF - adj_HL)]
fwrite(S6, file.path(RES, "study6_summary.csv"))
cat("\n== Study 6 (map cells, n = 1000)\n"); print(S6, digits = 3)

## ------------------------------------------------------------------ Study 7: plasmode with corrupted records
S7 <- rbindlist(lapply(c("logit", "loglog", "cloglog"), function(truth) rbindlist(lapply(c(1L, 5L, 10L), function(k) {
  x <- rd(7, sprintf("burn_%s_k%02d", truth, k)); N0 <- rd(7, sprintf("burn_logit_k%02d", k))
  data.table(truth, k, t(sapply(c("EF", "HL", "PH", "Stukel", "Cubic", "GiViTI"), function(t) rej(x[[t]]))),
             adj_EF = adjr(x$EF, N0$EF), adj_HL = adjr(x$HL, N0$HL), PH_NA = mean(is.na(x$PH)))
}))))
fwrite(S7, file.path(RES, "study7_summary.csv"))
cat("\n== Study 7 (burn plasmode with k corrupted records)\n"); print(S7, digits = 3)
