## pop.R -- population-level (N = 2e6) drift terms Delta and A of Theorem 2 and the predicted power
source("sg_core.R")
set.seed(20261002)
N  <- 2e6
X  <- sg_design(N)
ms <- c(2, 5, 10, 25, 50); ns <- c(1000, 4000, 16000)
p0 <- plogis(drop(X %*% c(0, 0.6, 0.5)))
R  <- sg_R(p0, X)
cat(sprintf("R (null design) = %.4f\n", R))
out <- list()
for (sc in c("cloglog", "loglog", "quadpos", "quadneg")) {
  pt <- sg_truth(X, sc)
  fit <- suppressWarnings(glm.fit(X, pt, family = quasibinomial()))
  ps <- fit$fitted.values; h <- pt - ps
  Rs <- sg_R(ps, X)
  for (m in ms) {
    ord <- order(ps); g <- integer(N); g[ord] <- (seq_len(N) - 1L) %/% m + 1L
    G  <- max(g)
    dl <- rowsum(h, g)[, 1]; h2 <- rowsum(h^2, g)[, 1]; ch <- rowsum((1 - 2 * ps) * h, g)[, 1]
    pb <- rowsum(ps, g)[, 1] / m; V <- m * pb * (1 - pb)
    Dl <- sum((dl^2 - h2) / V); A <- sum(ch / V)
    for (n in ns) {
      Gn <- n / m; f <- n / N
      sE <- sqrt(2 * Gn * (1 - 1 / m)); sH <- sqrt(2 * Gn * (1 - 1 / m) + Gn * Rs / m)
      tE <- Dl * f / sE; tH <- (Dl + A) * f / sH
      out[[length(out) + 1]] <- data.frame(scen = sc, n = n, m = m, Delta = Dl * f, A = A * f,
        sdEF = sE, sdHL = sH, thetaEF = tE, thetaHL = tH,
        powEF = pnorm(tE - qnorm(.95)), powHL = pnorm(tH - qnorm(.95)),
        powHL2 = pnorm(tH - qnorm(.975)) + pnorm(-tH - qnorm(.975)), R = Rs)
    }
  }
  cat(sc, "beta* =", round(fit$coefficients, 3), " R =", round(Rs, 4),
      " A/n (m=10) per record:", signif(out[[length(out)]]$A / 16000, 3), "\n")
}
res <- do.call(rbind, out)
write.csv(res, "out/pop_predictions.csv", row.names = FALSE)
print(res, digits = 3)
