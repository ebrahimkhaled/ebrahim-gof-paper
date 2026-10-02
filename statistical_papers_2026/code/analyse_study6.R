## analyse_study6.R -- the simulated map cells against the map's second-order prediction (Addendum 2, Study 6).
ROOT <- "."
suppressMessages(library(data.table))
RES <- file.path(ROOT, "results")
rej  <- function(p) mean(!is.na(p) & p < .05)
adjr <- function(p1, p0) { cr <- quantile(p0[!is.na(p0)], .05, type = 1); mean(!is.na(p1) & p1 <= cr) }
MAP <- fread(file.path(RES, "N2_sign_map.csv"))
files <- list.files(file.path(RES, "study6"), "gz$")
S6 <- rbindlist(lapply(files, function(f) {
  m <- regmatches(f, regexec("^map_([a-z]+)_l([0-9.]+)_mu([+-][0-9]+)_s([0-9]+\\.[0-9])\\.csv\\.gz$", f))[[1]]
  sd_ <- m[2]; lm_ <- as.numeric(m[3]); mu_ <- as.numeric(m[4]); sp_ <- as.numeric(m[5])
  x  <- fread(file.path(RES, "study6", f))
  nf <- file.path(RES, "study6", sprintf("map_cloglog_l1.00_mu%+d_s%.1f.csv.gz", as.integer(mu_), sp_))
  N0 <- if (file.exists(nf)) fread(nf) else NULL
  pr <- if (lm_ < 1) MAP$gain_n1000[MAP$side == sd_ & abs(MAP$lam - lm_) < 1e-9 & MAP$mu == mu_ & MAP$s == sp_] else 0
  data.table(side = sd_, lam = lm_, mu = mu_, s = sp_, predicted = pr, raw_EF = rej(x$EF), raw_HL = rej(x$HL),
             adj_EF = if (is.null(N0)) NA_real_ else adjr(x$EF, N0$EF), adj_HL = if (is.null(N0)) NA_real_ else adjr(x$HL, N0$HL))
}))
S6[, `:=`(sim_raw = raw_EF - raw_HL, sim_adj = adj_EF - adj_HL)]
setorder(S6, -lam, side, mu)
fwrite(S6, file.path(RES, "study6_summary.csv"))
print(S6, digits = 3)
ok <- S6[lam < 1 & abs(predicted) > 0.02]
cat("\nsign of simulated gain (size-adjusted where available, else raw) = sign of prediction:",
    ok[, sum(sign(ifelse(is.na(sim_adj), sim_raw, sim_adj)) == sign(predicted))], "of", nrow(ok), "\n")
