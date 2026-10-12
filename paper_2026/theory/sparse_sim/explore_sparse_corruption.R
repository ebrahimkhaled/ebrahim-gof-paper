## explore_sparse_corruption.R -- EXPLORATORY (after Addendum 2): in the many-small-groups regime, is the discarded
## part C (standardized: Z_L, to first order the Osius-Rojek score) fragile to corrupted records, while the pair
## statistic EF is not? Correct logistic model, eta = 0.6 x1 + 0.5 x2, n = 4000, groups of m records by fitted risk,
## k records with x1 multiplied by 4 or -4 before the fit. Normal references of Theorem S1: Z_EF and Z_HL one-sided
## upper 5%; Z_L two-sided 5%. 1000 replicates a cell, each written to disk.
here <- "./theory/sparse_sim"
suppressMessages({ library(parallel); library(data.table) })
OUT <- file.path(here, "out", "raw_corrupt"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
one <- function(rep, m, k, mult) {
  set.seed(7.7e8 + 1e5 * m + 1e3 * k + 10 * (mult > 0) + rep)
  X <- sg_design(4000); y <- rbinom(4000, 1, sg_truth(X, "null"))
  if (k > 0) { bad <- sample.int(4000, k); X[bad, 2] <- mult * X[bad, 2] }
  p <- glm.fit(X, y, family = binomial())$fitted.values
  s <- sg_stats(y, p, X, m)
  ## sigma_L^2 = sum c^2 tau / V^2 - t' M^-1 t, t = sum c U / V  (Theorem S1)
  o <- order(p); g <- integer(4000); g[o] <- (seq_len(4000) - 1L) %/% m + 1L
  G <- max(g); ng <- tabulate(g, G); if (ng[G] < m) { g[g == G] <- G - 1L; G <- G - 1L; ng <- tabulate(g, G) }
  v <- p * (1 - p); tau <- rowsum(v, g, reorder = TRUE)[, 1]; pb <- rowsum(p, g, reorder = TRUE)[, 1] / ng
  V <- ng * pb * (1 - pb); cc <- 1 - 2 * pb; U <- rowsum(v * X, g, reorder = TRUE)
  tt <- colSums(cc * U / V); sL <- sqrt(sum(cc^2 * tau / V^2) - drop(tt %*% solve(crossprod(X, v * X), tt)))
  c(rep = rep, zEF = unname(s["zEF"]), zHL = unname(s["zHL"]), zL = unname(s["C"]) / sL)
}
cl <- makePSOCKcluster(16)
clusterExport(cl, c("one", "here")); invisible(clusterEvalQ(cl, source(file.path(here, "sg_core.R"))))
res <- list()
for (m in c(5L, 25L)) for (k in c(0L, 1L, 5L, 10L)) for (mult in if (k == 0) 1 else c(4, -4)) {
  f <- file.path(OUT, sprintf("m%02d_k%02d_%s.csv", m, k, if (mult < 0) "neg" else "pos"))
  R <- if (file.exists(f)) fread(f) else { R <- rbindlist(lapply(parLapply(cl, 1:1000, one, m = m, k = k, mult = mult),
                                                                 function(v) as.data.table(as.list(v)))); fwrite(R, f); R }
  res[[length(res) + 1]] <- data.table(m, k, mult, EF = mean(R$zEF > qnorm(.95)), HL = mean(R$zHL > qnorm(.95)),
                                       ZL_two_sided = mean(abs(R$zL) > qnorm(.975)))
}
stopCluster(cl)
S <- rbindlist(res); fwrite(S, file.path(here, "out", "corrupt_summary.csv")); print(S, digits = 3)
