## run_addendum4.R -- Studies 11-13 of declarations/ADDENDUM4_large_sample.md.
##   Rscript run_addendum4.R 11|12|13 [B]
## One file per cell in results/study<k>/, one row per replicate, finished cells skipped on restart; each row carries
## fingerprints of its data. A second argument B runs a smoke test into results/smoke<k>/.
ROOT <- "."
suppressMessages({ library(parallel); library(data.table) })
args  <- commandArgs(TRUE)
study <- as.integer(args[1]); smoke <- length(args) > 1
B     <- if (smoke) as.integer(args[2]) else 2000L
OUT   <- file.path(ROOT, "results", paste0(if (smoke) "smoke" else "study", study))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
src <- readLines(file.path(ROOT, "code", "run_ef_studies.R"))
gen_src <- src[grep("^design_data <- function", src):(grep("^## one replicate", src) - 1)]
SH <- fread(file.path(ROOT, "results", "study9_shifts.csv"))

## the tests of Studies 11-13 on one fitted data set
large_tests <- function(y, X, giv = FALSE) {
  d  <- data.frame(y = y, X[, -1, drop = FALSE])
  fq <- bt_fit(list(d = d, f = as.formula(paste("y ~", paste(colnames(X)[-1], collapse = " + ")))))
  if (is.null(fq) || !isTRUE(fq$fit$converged)) return(NULL)
  p <- fq$p_raw; Xf <- fq$X; n <- length(y)
  o <- sample.int(n)
  ef10 <- ef_test(y[o], p[o], Xf[o, , drop = FALSE], 10)
  sP <- ls_stats(y, p, Xf, paul_G(n, y)); s25 <- ls_stats(y, p, Xf, floor(n / 25))
  stk <- bt_stukel(fq)
  qe <- quantile(fq$eta, c(.025, .975), names = FALSE); fw <- fq; fw$eta <- pmin(pmax(fq$eta, qe[1]), qe[2])
  jw <- bt_stukel_joint_stat(fw)
  pick <- c("p_HL_chisq", "p_EF_chisq", "p_HL_norm", "p_EF_norm", "p_L", "p_LW", "p_comb")
  c(G_paul = sP[["G"]], HL_stat = unname(ef10["HL"]),
    EF10 = unname(ef10["p_EF_chisq"]), HL10 = unname(ef10["p_HL_chisq"]),
    setNames(sP[pick], paste0(pick, "_P")), setNames(s25[pick], paste0(pick, "_25")),
    Stukel = unname(stk["Stk.joint"]),
    StukelW = if (is.finite(jw$chi)) pchisq(jw$chi, jw$df, lower.tail = FALSE) else NA_real_)
}

one <- function(rep, cell) {
  set.seed(cell$seed + rep)
  if (cell$design == "W") {                       # wide risk range (Online Resource 1, Table S3): eta = 1.2 x1 + 0.5 x2
    x1 <- runif(cell$n, -3, 3); x2 <- rbinom(cell$n, 1, .5); d <- list(X = cbind(x1, x2), eta = 1.2 * x1 + 0.5 * x2)
  } else d <- design_data(cell$design, cell$n)
  if (cell$truth == "shared") {                   # Study 11: records with similar risk share a deviation
    b <- ceiling(rank(d$eta, ties.method = "random") / 25)
    u <- rnorm(max(b), 0, cell$tau)
    y <- rbinom(cell$n, 1, plogis(d$eta + u[b]))
  } else {
    d$eta <- d$eta + cell$shift
    y <- rbinom(cell$n, 1, true_p(cell$truth, d))
  }
  X <- cbind(`(Intercept)` = 1, d$X); colnames(X)[-1] <- paste0("x", seq_len(ncol(X) - 1))
  if (cell$k > 0) { bad <- sample.int(cell$n, cell$k); X[bad, 2] <- cell$mult * X[bad, 2] }
  if (min(sum(y), length(y) - sum(y)) < 10) return(NULL)
  out <- tryCatch(large_tests(y, X), error = function(e) NULL)
  if (is.null(out)) return(NULL)
  c(rep = rep, fp_events = sum(y), fp_x = round(sum(X[, 2]), 6), out)
}

cells <- switch(as.character(study),
  `11` = CJ(design = c("D1", "W"), truth = c("logit", "shared"), tau = c(0, 0.3, 0.6), n = 4000L)[
           (truth == "logit" & tau == 0) | (truth == "shared" & tau > 0)],
  `12` = CJ(design = c("D1", "D2", "D4"), truth = c("logit", "cloglog", "loglog", "quad"), n = c(2000L, 5000L, 10000L, 20000L), tau = 0),
  `13` = rbind(data.table(design = "D1", truth = "logit", n = 5000L, tau = 0, k = 0L, mult = 1),
               CJ(design = "D1", truth = "logit", n = 5000L, tau = 0, k = c(5L, 25L), mult = c(4, -4))))
if (is.null(cells$k)) cells[, `:=`(k = 0L, mult = 1)]
cells <- merge(cells, SH[, .(design, truth, shift)], by = c("design", "truth"), all.x = TRUE, sort = FALSE)
cells[is.na(shift), shift := 0]
cells[, id := sprintf("%s_%s_n%05d_t%.1f_k%02d_m%s", design, truth, n, tau, k, gsub("-", "neg", mult))]
cells[, seed := 8e8 + study * 1e7 + .I * 1e4]

cl <- makePSOCKcluster(20)
clusterExport(cl, c("one", "large_tests", "gen_src", "ROOT"))
invisible(clusterEvalQ(cl, { eval(parse(text = gen_src)); source(file.path(ROOT, "code", "rivals.R"))
                             source(file.path(ROOT, "code", "large_sample.R")); NULL }))
for (j in seq_len(nrow(cells))) {
  cell <- as.list(cells[j]); f <- file.path(OUT, paste0(cell$id, ".csv.gz"))
  if (file.exists(f)) next
  t0 <- Sys.time()
  R <- rbindlist(lapply(clusterApplyLB(cl, seq_len(B), one, cell = cell),
                        function(v) if (is.null(v)) NULL else as.data.table(as.list(v))))
  dup <- sum(duplicated(R[, .(fp_events, fp_x, HL_stat)]))
  fwrite(R, f)
  r <- function(t) mean(R[[t]] < .05, na.rm = TRUE)
  cat(sprintf("%-34s B=%d dup=%d G=%d | HLc %.3f EFc %.3f HLn %.3f EFn %.3f L %.3f LW %.3f comb %.3f | Stk %.3f StkW %.3f [%.0fs]\n",
              cell$id, nrow(R), dup, as.integer(median(R$G_paul)), r("p_HL_chisq_P"), r("p_EF_chisq_P"), r("p_HL_norm_P"),
              r("p_EF_norm_P"), r("p_L_P"), r("p_LW_P"), r("p_comb_P"), r("Stukel"), r("StukelW"),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
stopCluster(cl)
cat("study", study, "done\n")
