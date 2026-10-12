## check_identity.R -- numerical check of the decomposition HL = G + C + Q, EF = G + Q, with
## Q = sum_g V_g^{-1} sum_{i != j in g} (y_i - pbar_g)(y_j - pbar_g), for any partition and any fitted risks.
source("./code/ef_exact.R")
set.seed(3)
for (G in c(2, 10, 50, 250)) {
  n <- 1000; x <- runif(n, -3, 3); y <- rbinom(n, 1, plogis(0.6 * x))
  X <- cbind(1, x); p <- glm.fit(X, y, family = binomial())$fitted.values
  t <- ef_test(y, p, X, G); g <- ef_groups(p, G)
  Q <- sum(sapply(split(seq_len(n), g), function(I) {
    pb <- mean(p[I]); V <- length(I) * pb * (1 - pb); d <- y[I] - pb
    (sum(d)^2 - sum(d^2)) / V }))
  cat(sprintf("G = %3d: HL - (G + C + Q) = %.2e   EF - (G + Q) = %.2e\n", G, t["HL"] - (G + t["C"] + Q), t["EF"] - (G + Q)))
}
