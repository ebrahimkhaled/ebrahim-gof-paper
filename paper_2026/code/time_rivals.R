## time_rivals.R -- Study 10 of declarations/ADDENDUM3_inrange_calibrated.md: elapsed time, one core, of EF (with HL),
## Stukel, the le Cessie-van Houwelingen test in its original O(n^3) form (the form of the public code
## smwrStats::leCessie.test: n x n hat matrix and n x n products) and BAGofT (BAGofT package, default settings) on one
## data set of design D1. Median of repeated runs when a run is short.
ROOT <- "."
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) { RhpcBLASctl::blas_set_num_threads(1); RhpcBLASctl::omp_set_num_threads(1) }
suppressMessages(library(data.table))
source(file.path(ROOT, "code", "rivals.R"))
src <- readLines(file.path(ROOT, "code", "run_ef_studies.R"))
eval(parse(text = src[grep("^design_data <- function", src):(grep("^## one replicate", src) - 1)]))

## le Cessie-van Houwelingen, O(n^3) form: H = V X (X'VX)^{-1} X', M = (I-H)' R (I-H)
lecessie_n3 <- function(fit) {
  y <- fit$y; p <- fitted(fit); r <- y - p; N <- length(y)
  covs <- model.frame(fit)[, -1, drop = FALSE]
  dl <- lapply(covs, function(x) if (is.numeric(x)) as.numeric(0.5 * stats::dist(scale(x))^2)
               else { xx <- as.numeric(as.factor(x)); nc <- length(unique(xx))
                      as.numeric((stats::dist(xx, method = "manhattan") != 0) * nc / (nc - 1)) })
  D <- matrix(0, N, N); D[lower.tri(D)] <- sqrt(rowSums(as.data.frame(dl))); D <- D + t(D)
  R <- pmax(1 - D / mean(D), 0)
  Q <- sum(as.numeric(r %*% R) * r)
  X <- model.matrix(fit); mu2 <- p * (1 - p)
  H <- (mu2 * X) %*% solve(crossprod(X, mu2 * X)) %*% t(X)
  IH <- diag(N) - H
  M <- t(IH) %*% R %*% IH; M <- (M + t(M)) / 2
  V <- diag(mu2)
  EQ <- sum(diag(M) * mu2); mu4 <- mu2 * (1 - 3 * mu2)
  VarQ <- sum(diag(M)^2 * (mu4 - 3 * mu2^2)) + 2 * sum(diag(M %*% V %*% M %*% V))
  stats::pchisq(Q * 2 * EQ / VarQ, 2 * EQ^2 / VarQ, lower.tail = FALSE)
}

tm <- function(f, reps) median(sapply(seq_len(reps), function(i) system.time(f())[["elapsed"]]))
out <- list()
for (n in c(1000L, 5000L)) {
  set.seed(20261002 + n)
  d <- design_data("D1", n); dat <- data.frame(y = rbinom(n, 1, plogis(d$eta)), d$X)
  fit <- glm(y ~ x1 + x2, binomial, data = dat); p <- fitted(fit); X <- model.matrix(fit); y <- dat$y
  fq <- bt_fit(list(d = dat, f = y ~ x1 + x2))
  ## fast tests: time 200 calls at once and divide, so the timer resolution does not matter
  out[[length(out) + 1]] <- data.table(n, test = "EF and HL", sec = tm(function() for (i in 1:200) ef_test(y, p, X, 10), 11) / 200)
  out[[length(out) + 1]] <- data.table(n, test = "Stukel", sec = tm(function() for (i in 1:200) bt_stukel(fq), 11) / 200)
  out[[length(out) + 1]] <- data.table(n, test = "le Cessie-van Houwelingen (O(n^3) form)",
                                       sec = tm(function() lecessie_n3(fit), if (n == 1000) 3 else 1))
  if (n == 1000)
    out[[length(out) + 1]] <- data.table(n, test = "BAGofT",
      sec = tm(function() BAGofT::BAGofT(testModel = BAGofT::testGlmBi(formula = y ~ x1 + x2, link = "logit"), data = dat), 1))
  print(rbindlist(out)); flush.console()
  fwrite(rbindlist(out), file.path(ROOT, "results", "study10_timing.csv"))
}
R <- rbindlist(out)
R[, times_EF := sec / sec[test == "EF and HL"], by = n]
fwrite(R, file.path(ROOT, "results", "study10_timing.csv"))
print(R)
