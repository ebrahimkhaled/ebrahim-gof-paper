## analyse_claims.R -- the verdict on every claim written in the declarations before its study ran, recomputed from the
## per-replicate result files: STUDY_SPEC (72603ab2), Addendum 1 (42f8c0cd), Addendum 2 (c9de6124) and Addendum 3
## (8594a3e6). Writes results/claims_ledger.csv (one row per claim: declaration, claim, number, verdict) and the summaries
## of Studies 8 and 9 that the revised paper uses. Missing p-values (a test that declined) count as no rejection.
ROOT <- "."
suppressMessages(library(data.table))
source(file.path(ROOT, "code", "ef_exact.R"))
RES <- file.path(ROOT, "results")
rej  <- function(p) mean(!is.na(p) & p < .05)
adjr <- function(p1, p0) { cr <- quantile(p0[!is.na(p0)], .05, type = 1); mean(!is.na(p1) & p1 <= cr) }
rd   <- function(dir, id) fread(file.path(RES, dir, paste0(id, ".csv.gz")))
L <- list()
claim <- function(decl, id, text, number, held) L[[length(L) + 1]] <<- data.table(declaration = decl, claim = id, text, number,
                                                                                verdict = if (isTRUE(held)) "held" else "failed")

## ------------------------------------------------------------------ STUDY_SPEC: Studies 1-3
S1 <- rbindlist(lapply(list.files(file.path(RES, "study1"), "gz$"), function(f) {
  x <- fread(file.path(RES, "study1", f)); m <- regmatches(f, regexec("^(D\\d)_logit_n(\\d+)", f))[[1]]
  data.table(design = m[2], n = as.integer(m[3]), EFx = rej(x$p_EF_exact), EFc = rej(x$p_EF_chisq), HLc = rej(x$p_HL_chisq))
}))
bad <- S1[n >= 200 & abs(EFx - .05) > .015]
claim("Spec", "S1a", "EF with its finite-sample law holds 5% within 0.015 in every cell with n >= 200",
      if (nrow(bad)) paste(sprintf("%s n=%d: %.4f", bad$design, bad$n, bad$EFx), collapse = "; ") else "all within", nrow(bad) == 0)
dep <- S1[abs(EFc - .05) > .015 | abs(HLc - .05) > .015]
claim("Spec", "S1b", "EF and HL with chi-square(G-2) depart from 5% by more than 0.015 in at least one cell",
      if (nrow(dep)) paste(sprintf("%s n=%d: EF %.4f HL %.4f", dep$design, dep$n, dep$EFc, dep$HLc), collapse = "; ") else "none", nrow(dep) > 0)
S2 <- fread(file.path(RES, "study2_summary.csv"))
S2a <- S2[abs(A) > 2 * seC]
claim("Spec", "S2a", "sign of mean C equals sign of A in every cell where |A| > 2 Monte Carlo SE of mean C",
      sprintf("%d of %d", S2a[, sum(sign(meanC) == sign(A))], nrow(S2a)), S2a[, all(sign(meanC) == sign(A))])
for (tr in c("cloglog", "loglog")) {
  x <- S2[truth == tr]
  claim("Spec", paste0("S2b-", tr), sprintf("EF >= HL (size-adjusted) against the %s link", tr),
        sprintf("%d of %d cells", x[, sum(adj_EFchi >= adj_HLchi)], nrow(x)), x[, all(adj_EFchi >= adj_HLchi)])
}
x <- S2[truth == "probit"]
claim("Spec", "S2b-probit", "EF about equal to HL against the probit link (|difference| <= 0.01)",
      sprintf("%d of %d cells; largest |difference| %.3f", x[, sum(abs(adj_EFchi - adj_HLchi) <= .01)], nrow(x), x[, max(abs(adj_EFchi - adj_HLchi))]),
      x[, all(abs(adj_EFchi - adj_HLchi) <= .01)])
x <- S2[truth == "quad"]
claim("Spec", "S2b-quad", "EF <= HL against the omitted square",
      sprintf("%d of %d cells (A < 0 in all of them)", x[, sum(adj_EFchi <= adj_HLchi)], nrow(x)), x[, all(adj_EFchi <= adj_HLchi)])
S3 <- rbindlist(lapply(list.files(file.path(RES, "study3"), "gz$"), function(f) {
  x <- fread(file.path(RES, "study3", f)); m <- regmatches(f, regexec("_k(\\d+)_m(neg)?(\\d)", f))[[1]]
  data.table(k = as.integer(m[2]), mult = if (m[3] == "neg") -4 else 4, EFx = rej(x$p_EF_exact), Stk = rej(x$p_Stukel))
}))
claim("Spec", "S3a", "EF (finite-sample law) stays at or below 0.10 up to k = 10 exaggerated records",
      sprintf("largest %.4f", S3[mult == 4, max(EFx)]), S3[mult == 4, all(EFx <= .10)])
claim("Spec", "S3b", "Stukel exceeds 0.10 at k = 1",
      sprintf("x4: %.4f; x(-4): %.4f", S3[k == 1 & mult == 4, Stk], S3[k == 1 & mult == -4, Stk]), S3[k == 1, all(Stk > .10)])

## ------------------------------------------------------------------ Addendum 1: plasmode
P <- function(id) fread(file.path(RES, "plasmode", paste0(id, ".csv.gz")))
p0 <- P("logit_k00")
for (tr in c("loglog", "cloglog")) {
  x <- P(paste0(tr, "_k00")); e <- adjr(x$p_EF_chisq, p0$p_EF_chisq); h <- adjr(x$p_HL_chisq, p0$p_HL_chisq)
  claim("Add1", paste0("P-", tr), sprintf("%s-type truth: EF %s powerful than HL (size-adjusted)", tr, if (tr == "loglog") "more" else "less"),
        sprintf("EF %.3f, HL %.3f", e, h), if (tr == "loglog") e > h else e < h)
}
cor_ <- rbindlist(lapply(c(1, 2, 5, 10), function(k) { x <- P(sprintf("logit_k%02d", k))
  data.table(k, EF = rej(x$p_EF_chisq), HL = rej(x$p_HL_chisq), Stk = rej(x$p_Stukel)) }))
claim("Add1", "P-corr-a", "EF and HL stay at or below 0.10 with k <= 10 corrupted patients",
      sprintf("largest EF %.4f, HL %.4f", max(cor_$EF), max(cor_$HL)), all(cor_$EF <= .10 & cor_$HL <= .10))
claim("Add1", "P-corr-b", "Stukel exceeds 0.10 at k = 1", sprintf("%.4f", cor_[k == 1, Stk]), cor_[k == 1, Stk > .10])

## ------------------------------------------------------------------ Addendum 2: Studies 4-6
S4 <- fread(file.path(RES, "study4_summary.csv")); D1 <- S4[design == "D1"]
claim("Add2", "a", "Stukel and cubic LR exceed 0.10 at k = 1 (x4) in D1",
      sprintf("Stukel %.3f, cubic %.3f", D1[k == 1 & mult == 4, Stukel], D1[k == 1 & mult == 4, Cubic]), D1[k == 1 & mult == 4, Stukel > .10 & Cubic > .10])
claim("Add2", "b", "GiViTI exceeds 0.10 by k = 2 (x4) in D1",
      sprintf("k=1 %.3f, k=2 %.3f", D1[k == 1 & mult == 4, GiViTI], D1[k == 2 & mult == 4, GiViTI]), D1[k %in% 1:2 & mult == 4, any(GiViTI > .10)])
claim("Add2", "c", "EF, HL and Pigeon-Heyse at or below 0.10 up to k = 10 (x4), both designs",
      sprintf("largest %.4f", S4[mult == 4, max(EF, HL, PH)]), S4[mult == 4, all(EF <= .10 & HL <= .10 & PH <= .10)])
claim("Add2", "d", "Tsiatis and Xie exceed 0.10 by k = 10 (x4) in D1",
      sprintf("Tsiatis %.3f, Xie %.3f", D1[k == 10 & mult == 4, Tsiatis], D1[k == 10 & mult == 4, Xie]), D1[k == 10 & mult == 4, Tsiatis > .10 & Xie > .10])
S5 <- fread(file.path(RES, "study5_summary.csv"))[truth != "logit"]
S5[, best := c("EF", "HL", "PH")[max.col(cbind(adj_EF, adj_HL, adj_PH), ties.method = "first")]]
vis <- S5[round(abs(gain_adj), 3) > .010]   # rounded: two cells differ by exactly 0.010
claim("Add2", "e", "among EF, HL, Pigeon-Heyse: EF best where A < 0, HL or Pigeon-Heyse best where A > 0 (cells with |EF - HL| > 0.01)",
      sprintf("A < 0: %d of %d; A > 0: %d of %d", vis[A < 0, sum(best == "EF")], vis[A < 0, .N], vis[A > 0, sum(best != "EF")], vis[A > 0, .N]),
      vis[A < 0, all(best == "EF")] && vis[A > 0, all(best != "EF")])
f5 <- S5[abs(A) > .2]
claim("Add2", "f", "sign of mean C equals sign of A in at least 90% of cells with |A| > 0.2",
      sprintf("%d of %d", f5[, sum(sign(meanC) == sign(A))], nrow(f5)), f5[, mean(sign(meanC) == sign(A))] >= .9)
S6 <- fread(file.path(RES, "study6_summary.csv"))
decl6 <- data.table(side = c("cloglog", "cloglog", "loglog", "loglog", "cloglog", "cloglog", "loglog"),
                    lam = c(0, 0, 0, 0, 0.25, 0.5, 0.25), mu = c(2, -2, -2, 2, 1, 0, -1), s = c(2, 2, 2, 2, 1, 1, 1))
x6 <- merge(decl6, S6, by = c("side", "lam", "mu", "s"))[abs(predicted) > .02]
x6[, sim := fifelse(is.na(sim_adj), sim_raw, sim_adj)]
claim("Add2", "map", "sign of the simulated gain equals the sign of the map's prediction where |prediction| > 0.02",
      sprintf("%d of %d declared cells", x6[, sum(sign(sim) == sign(predicted))], nrow(x6)), x6[, all(sign(sim) == sign(predicted))])

## ------------------------------------------------------------------ Addendum 3: Studies 8 and 9
tests8 <- c("EF", "HL", "PH", "Tsiatis", "Xie", "Stukel", "StukelW", "Cubic", "GiViTI")
if (length(list.files(file.path(RES, "study8"), "gz$")) == 40) {
  S8 <- rbindlist(lapply(list.files(file.path(RES, "study8"), "gz$"), function(f) {
    x <- fread(file.path(RES, "study8", f)); id <- sub(".csv.gz", "", f)
    m <- regmatches(id, regexec("^(D\\d|burn)_(none|mult|yflip|cap)(neg)?(\\d)?_k(\\d+)$", id))[[1]]
    err <- m[3]; mult <- if (err == "mult") (if (m[4] == "neg") -1 else 1) * as.numeric(m[5]) else NA_real_
    data.table(design = m[2], err, mult, k = as.integer(m[6]), B = nrow(x), t(sapply(tests8, function(t) rej(x[[t]]))))
  }))
  S8[, error := fcase(err == "none", "none", err == "mult" & mult == 4, "x4", err == "mult" & mult == -4, "x(-4)",
                      err == "mult" & mult == -1, "sign", err == "yflip", "y flipped", err == "cap", "tbsa x4, capped")]
  setorder(S8, design, error, k); fwrite(S8, file.path(RES, "study8_summary.csv"))
  g <- S8[design == "D1" & error == "sign" & k <= 10]
  claim("Add3", "g", "in-range sign errors (D1): EF, HL, Pigeon-Heyse <= 0.10 up to k = 10, and one of Stukel, cubic LR, GiViTI > 0.10 by k = 10",
        sprintf("grouped largest %.3f; at k=10 Stukel %.3f, cubic %.3f, GiViTI %.3f", g[, max(EF, HL, PH)], g[k == 10, Stukel], g[k == 10, Cubic], g[k == 10, GiViTI]),
        g[, all(EF <= .10 & HL <= .10 & PH <= .10)] && g[k == 10, max(Stukel, Cubic, GiViTI) > .10])
  h <- S8[design == "D1" & error == "y flipped" & k <= 20]
  claim("Add3", "h", "response misclassification (D1): EF, HL, Pigeon-Heyse <= 0.10 up to k = 20",
        sprintf("largest %.3f", h[, max(EF, HL, PH)]), h[, all(EF <= .10 & HL <= .10 & PH <= .10)])
  i_ <- S8[design == "D1" & error == "x(-4)" & k == 1]
  claim("Add3", "i", "Stukel-W below Stukel at k = 1, x(-4), D1", sprintf("Stukel-W %.3f, Stukel %.3f", i_$StukelW, i_$Stukel), i_$StukelW < i_$Stukel)
  j_ <- S8[design == "burn"]
  claim("Add3", "j", "plasmode with capped tbsa: EF <= 0.10 up to k = 10", sprintf("largest %.3f", j_[, max(EF)]), j_[, all(EF <= .10)])
}
if (length(list.files(file.path(RES, "study9"), "gz$")) == 90) {
  SH <- fread(file.path(RES, "study9_shifts.csv"))
  Afile <- file.path(RES, "alignment_calibrated.csv")
  if (!file.exists(Afile)) {
    src <- readLines(file.path(ROOT, "code", "run_ef_studies.R"))
    eval(parse(text = src[grep("^design_data <- function", src):(grep("^## one replicate", src) - 1)]))
    AA <- rbindlist(lapply(seq_len(nrow(SH[truth != "logit"])), function(r) {
      s <- SH[truth != "logit"][r]; set.seed(777); d <- design_data(s$design, 1e6); d$eta <- d$eta + s$shift
      pt <- true_p(s$truth, d); X <- cbind(1, d$X)
      ps <- suppressWarnings(glm.fit(X, pt, family = quasibinomial()))$fitted.values
      g <- ef_groups(ps, 10); ng <- tabulate(g, 10)
      dg <- as.numeric(rowsum(pt - ps, g, reorder = TRUE)) / ng; pb <- as.numeric(rowsum(ps, g, reorder = TRUE)) / ng
      data.table(design = s$design, truth = s$truth, A = sum((1 - 2 * pb) * dg / (pb * (1 - pb))), event_rate = mean(pt))
    }))
    fwrite(AA, Afile)
  }
  AA <- fread(Afile)
  tests9 <- c("EF", "HL", "PH", "Tsiatis", "Xie", "Stukel", "StukelW", "Cubic")
  S9 <- rbindlist(lapply(c("D1", "D2", "D3", "D4", "D5"), function(design) rbindlist(lapply(c(500L, 1000L, 2000L), function(n) {
    N0 <- rd("study9", sprintf("%s_logit_n%05d", design, n))
    rbindlist(lapply(c("logit", "cloglog", "loglog", "probit", "quad", "inter"), function(truth) {
      x <- rd("study9", sprintf("%s_%s_n%05d", design, truth, n))
      row <- data.table(design, truth, n, B = nrow(x), meanC = mean(x$C))
      for (t in tests9) { set(row, j = paste0("raw_", t), value = rej(x[[t]])); set(row, j = paste0("adj_", t), value = adjr(x[[t]], N0[[t]])) }
      ## Monte Carlo SE of the paired size-adjusted EF - HL difference
      cE <- quantile(N0$EF, .05, type = 1); cH <- quantile(N0$HL, .05, type = 1)
      dEH <- (x$EF <= cE) - (x$HL <= cH); set(row, j = "se_gain", value = sd(dEH) / sqrt(length(dEH)))
      row }))
  }))))
  S9 <- merge(S9, AA, by = c("design", "truth"), all.x = TRUE)
  S9[, gain_adj := adj_EF - adj_HL]
  fwrite(S9, file.path(RES, "study9_summary.csv"))
  A9 <- S9[truth != "logit"]
  A9[, best := c("EF", "HL", "PH")[max.col(cbind(adj_EF, adj_HL, adj_PH), ties.method = "first")]]
  A9[, all_low := pmax(adj_EF, adj_HL, adj_PH) <= .05]
  v9 <- A9[round(abs(gain_adj), 3) > .010]
  claim("Add3", "e'", "event rate held fixed: EF best of EF, HL, Pigeon-Heyse in every cell with A < 0 and |EF - HL| > 0.01",
        sprintf("%d of %d (cells where all three have power <= 0.05: %d)", v9[A < 0, sum(best == "EF")], v9[A < 0, .N], v9[A < 0, sum(all_low)]),
        v9[A < 0, all(best == "EF")])
  claim("Add3", "f'", "sign of the EF - HL difference opposite to the sign of A in at least 90% of cells with |EF - HL| > 0.01",
        sprintf("%d of %d", v9[, sum(sign(gain_adj) == -sign(A))], nrow(v9)), v9[, mean(sign(gain_adj) == -sign(A))] >= .9)
}
LL <- rbindlist(L); fwrite(LL, file.path(RES, "claims_ledger.csv"))
print(LL[, .(declaration, claim, verdict, number)], right = FALSE)
