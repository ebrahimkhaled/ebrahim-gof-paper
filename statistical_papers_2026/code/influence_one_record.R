## influence_one_record.R -- Proposition 3 illustrated: one record's covariate multiplied by k, everything else fixed.
## Reports the change in HL, C, EF and Stukel's statistic (a) at the clean fit's coefficients and (b) after
## refitting, against the bound of Proposition 3 at fixed coefficients,
##   |dEF| <= sum_k [ (2|r_k| + d_k) d_k + |b_k| d_k ] + (pbar terms),  d_k = 2 / sqrt(V_k),
## evaluated here in its simple envelope  G (2 max|r| + 2/sqrt(Vmin) + 1) * 2/sqrt(Vmin).
ROOT <- "."
source(file.path(ROOT, "code", "ef_exact.R"))
set.seed(20261002)
n <- 1000; G <- 10
x1 <- runif(n, -3, 3); x2 <- rbinom(n, 1, .5); y <- rbinom(n, 1, plogis(0.6 * x1 + 0.5 * x2))
X0 <- cbind(1, x1, x2); f0 <- glm.fit(X0, y, family = binomial()); b0 <- f0$coefficients
## the record to corrupt: a non-event at x1 near 1.5; multiplying x1 pushes it to ever higher predicted risk,
## so its outcome becomes ever more surprising -- the hostile case for a calibration test
i <- which(y == 0)[which.min(abs(x1[y == 0] - 1.5))]
stat <- function(X, beta) {
  eta <- as.numeric(X %*% beta); p <- plogis(eta)
  t <- ef_test(y, p, X); s <- stukel_joint(y, p, X, eta)
  c(HL = unname(t["HL"]), C = unname(t["C"]), EF = unname(t["EF"]), Stukel = unname(s["Stukel"]))
}
base <- stat(X0, b0)
z0 <- ef_parts(y, f0$fitted.values, X0, G)
env <- G * (2 * max(abs(z0$r)) + 2 / sqrt(min(z0$V)) + 1) * 2 / sqrt(min(z0$V))
out <- do.call(rbind, lapply(c(2, 4, 8, 16, 100, 1000), function(k) {
  X <- X0; X[i, 2] <- k * X0[i, 2]
  fixed <- stat(X, b0) - base
  f1 <- suppressWarnings(glm.fit(X, y, family = binomial())); refit <- stat(X, f1$coefficients) - base
  data.frame(k = k, x1 = round(X[i, 2], 2),
             dHL_fixed = fixed["HL"], dEF_fixed = fixed["EF"], dStukel_fixed = fixed["Stukel"],
             dHL_refit = refit["HL"], dEF_refit = refit["EF"], dStukel_refit = refit["Stukel"], row.names = NULL)
}))
cat("record", i, " x1 =", round(x1[i], 3), " y =", y[i], "\nbase:", round(base, 3),
    "\nenvelope of Proposition 3 at fixed coefficients:", round(env, 2), "\n\n")
print(round(out, 3))
write.csv(out, file.path(ROOT, "results", "influence_one_record.csv"), row.names = FALSE)
