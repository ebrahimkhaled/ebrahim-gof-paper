## check_mean_convergence.R -- does E(HL) approach the first-order limit sum(lam) as n grows?
## If the gap shrinks with n, the 6% excess at n = 500 is a finite-sample effect; if it does not, the
## covariance in ef_exact.R misses a term (for example the randomness of the group boundaries).
source("./code/ef_exact.R")
suppressMessages(library(parallel))
one <- function(rep, n) {
  set.seed(1e6 + 1e4 * n %% 9973 + rep)
  x1 <- runif(n, -3, 3); x2 <- rbinom(n, 1, .5)
  y  <- rbinom(n, 1, plogis(0.6 * x1 + 0.5 * x2))
  X  <- cbind(1, x1, x2)
  f  <- glm.fit(X, y, family = binomial())
  ef_test(y, f$fitted.values, X)[c("HL", "EF", "C", "sum_lam", "p_HL_exact", "p_EF_exact", "p_HL_chisq")]
}
cl <- makePSOCKcluster(20)
clusterEvalQ(cl, source("./code/ef_exact.R"))
for (n in c(500, 2000, 10000, 50000)) {
  R <- do.call(rbind, parLapply(cl, 1:2000, one, n = n))
  cat(sprintf("n=%6d  E(HL)=%.3f  E(EF)=%.3f  sum(lam)=%.3f  ratio=%.3f  var(HL)=%.2f vs 2*sum(lam^2)~?  size HLx=%.4f EFx=%.4f HLchi=%.4f\n",
              n, mean(R[, "HL"]), mean(R[, "EF"]), mean(R[, "sum_lam"]), mean(R[, "HL"]) / mean(R[, "sum_lam"]),
              var(R[, "HL"]), mean(R[, "p_HL_exact"] < .05), mean(R[, "p_EF_exact"] < .05), mean(R[, "p_HL_chisq"] < .05)))
}
stopCluster(cl)
