## run_ef_studies.R -- Studies 1-3 of declarations/STUDY_SPEC_EF_ESJ.md (hash 72603ab2).
##   Rscript run_ef_studies.R 1|2|3
## Each cell is written to results/study<k>/<cell>.csv.gz as soon as it finishes, one row per replicate,
## and a finished cell is skipped on restart.
ROOT <- "."
suppressMessages({ library(parallel); library(data.table) })
study <- as.integer(commandArgs(TRUE)[1])
OUT <- file.path(ROOT, "results", paste0("study", study)); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

## covariates and linear predictor of each design
design_data <- function(design, n) {
  switch(design,
    D1 = { x1 <- runif(n, -3, 3); x2 <- rbinom(n, 1, .5); list(X = cbind(x1, x2), eta = 0.6 * x1 + 0.5 * x2) },
    D2 = { x1 <- runif(n, -3, 3); x2 <- rbinom(n, 1, .5); list(X = cbind(x1, x2), eta = -2.8 + 0.6 * x1 + 0.5 * x2) },
    D3 = { x1 <- rexp(n) - 1;    x2 <- rbinom(n, 1, .5); list(X = cbind(x1, x2), eta = -0.3 + 0.9 * x1 + 0.5 * x2) },
    D4 = { x1 <- runif(n, -3, 3); x2 <- rbinom(n, 1, .5); list(X = cbind(x1, x2), eta = 1.6 * x1 + 0.8 * x2) },
    D5 = { Z <- cbind(matrix(rnorm(4 * n), n), matrix(rbinom(2 * n, 1, .3), n)); list(X = Z, eta = 0.4 * rowSums(Z)) })
}
true_p <- function(truth, d) {
  e <- d$eta; x1 <- d$X[, 1]; x2 <- d$X[, 2]
  switch(truth,
    logit   = plogis(e),
    cloglog = 1 - exp(-exp(e)),
    loglog  = exp(-exp(-e)),
    probit  = pnorm(e),
    quad    = plogis(e + 0.15 * x1^2),
    inter   = plogis(e + 0.5 * x1 * x2))
}

## one replicate: generate, optionally corrupt, fit the linear logit model, run every test
one_rep <- function(rep, cell) {
  set.seed(cell$seed + rep)
  d <- design_data(cell$design, cell$n)
  y <- rbinom(cell$n, 1, true_p(cell$truth, d))
  X <- cbind(1, d$X)
  if (cell$k > 0) {                                   # Study 3: records recorded with a wrong covariate
    bad <- sample.int(cell$n, cell$k)
    X[bad, 2] <- cell$mult * X[bad, 2]
  }
  if (min(sum(y), cell$n - sum(y)) < 2) return(NULL)
  f <- suppressWarnings(glm.fit(X, y, family = binomial()))
  if (!isTRUE(f$converged)) return(NULL)
  p <- f$fitted.values; eta <- f$linear.predictors
  c(rep = rep, ef_test(y, p, X), stukel_joint(y, p, X, eta), events = sum(y))
}

cells <- switch(study,
  `1` = CJ(design = c("D1", "D2", "D3", "D4", "D5"), truth = "logit", n = c(100L, 200L, 500L, 1000L, 5000L)),
  `2` = rbind(CJ(design = c("D1", "D3", "D4"), truth = c("cloglog", "loglog", "probit", "quad", "inter"),
                 n = c(200L, 500L, 1000L, 2000L)),
              CJ(design = c("D1", "D3", "D4"), truth = "logit", n = 2000L)),
  `3` = rbind(data.table(design = "D1", truth = "logit", n = 1000L, k = 0L, mult = 1),
              CJ(design = "D1", truth = "logit", n = 1000L, k = c(1L, 2L, 5L, 10L), mult = c(4, -4))))
if (is.null(cells$k)) { cells[, k := 0L]; cells[, mult := 1] }
cells[, id := sprintf("%s_%s_n%05d_k%02d_m%s", design, truth, n, k, gsub("-", "neg", mult))]
cells[, seed := 20261002L * 0 + study * 1e8 + .I * 1e5]

cl <- makePSOCKcluster(20)
clusterExport(cl, c("design_data", "true_p", "one_rep", "ROOT"))
invisible(clusterEvalQ(cl, source(file.path(ROOT, "code", "ef_exact.R"))))
B <- 2000L
for (j in seq_len(nrow(cells))) {
  cell <- as.list(cells[j]); f <- file.path(OUT, paste0(cell$id, ".csv.gz"))
  if (file.exists(f)) next
  t0 <- Sys.time()
  R <- rbindlist(lapply(clusterApplyLB(cl, seq_len(B), one_rep, cell = cell),
                        function(v) if (is.null(v)) NULL else as.data.table(as.list(v))))
  fwrite(R, f)
  cat(sprintf("%-30s B=%d  size/power: HLchi %.3f HLx %.3f EFchi %.3f EFx %.3f Stk %.3f  [%.0fs]\n", cell$id, nrow(R),
              mean(R$p_HL_chisq < .05), mean(R$p_HL_exact < .05), mean(R$p_EF_chisq < .05),
              mean(R$p_EF_exact < .05), mean(R$p_Stukel < .05, na.rm = TRUE),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
stopCluster(cl)
cat("study", study, "done\n")
