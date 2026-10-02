## theorem2_predictions.R -- what Theorem 2 predicts for every Study-2 cell, from a pseudo-true fit on 10^6 records.
## For a fixed departure the linear logit fit converges to a pseudo-true p*. With groups at the quantiles of p*,
##   d_g = mean_g(p_true - p*),  delta_g = n_g d_g,  V_g = n_g pbar_g (1 - pbar_g),
## so A(delta) = sum (1 - 2 pbar_g) delta_g / V_g = sum (1 - 2 pbar_g) d_g / (pbar_g (1 - pbar_g)) does not depend on n:
## E(C) -> A for every n. The power of EF and HL follows from the Theorem-2 law with mean m = V^{-1/2} delta.
ROOT <- "."
source(file.path(ROOT, "code", "ef_exact.R"))
suppressMessages(library(data.table))
src <- readLines(file.path(ROOT, "code", "run_ef_studies.R"))
eval(parse(text = src[grep("^design_data <- function", src):(grep("^## one replicate", src) - 1)]))

G <- 10; N <- 1e6
## upper alpha point of the null law sum lam chi2_1 - kappa, by bisection on qf_upper
crit <- function(lam, delta, kappa, alpha = .05) uniroot(function(q) qf_upper(q, lam, delta, kappa) - alpha,
                                                       c(-kappa + 1e-6, 200), tol = 1e-6)$root
pred <- list()
for (design in c("D1", "D3", "D4")) for (truth in c("cloglog", "loglog", "probit", "quad", "inter")) {
  set.seed(777)
  d <- design_data(design, N); pt <- true_p(truth, d); X <- cbind(1, d$X)
  y <- rbinom(N, 1, pt)
  f <- glm.fit(X, y, family = binomial()); ps <- f$fitted.values
  g <- ef_groups(ps, G)
  dg <- as.numeric(rowsum(pt - ps, g, reorder = TRUE)) / tabulate(g, G)
  pb <- as.numeric(rowsum(ps, g, reorder = TRUE)) / tabulate(g, G)
  A  <- sum((1 - 2 * pb) * dg / (pb * (1 - pb)))
  ## covariance per record (scale-free), then scaled to each n
  z1 <- ef_parts(y, ps, X, G)                     # Sigma does not depend on n (ratios of sums)
  es <- eigen(z1$Sigma, symmetric = TRUE); keep <- es$values > 1e-10 * max(es$values)
  lam <- es$values[keep]; P <- es$vectors
  for (n in c(200, 500, 1000, 2000)) {
    ng <- n / G
    V  <- ng * pb * (1 - pb)
    m  <- ng * dg / sqrt(V)
    b  <- (1 - 2 * pb) / sqrt(V)
    beta <- as.numeric(crossprod(P, b)); eta <- as.numeric(crossprod(P, m))
    ## null laws and their 5% points
    k0 <- sum(beta[keep]^2) / 4; d0 <- beta[keep]^2 / (4 * lam)
    qEF <- crit(lam, d0, k0); qHL <- crit(lam, rep(0, length(lam)), 0)
    ## alternative laws (Theorem 2)
    kA <- sum(beta[keep]^2) / 4 - sum(eta[!keep]^2 - beta[!keep] * eta[!keep])
    dA <- (2 * eta[keep] - beta[keep])^2 / (4 * lam)
    pEF <- qf_upper(qEF, lam, dA, kA)
    pHL <- qf_upper(qHL, lam, eta[keep]^2 / lam, -sum(eta[!keep]^2))
    pred[[length(pred) + 1]] <- data.table(design, truth, n, A = A, pred_power_EF = pEF, pred_power_HL = pHL,
                                           pred_gain = pEF - pHL)
  }
  cat(design, truth, "A =", round(A, 3), "\n")
}
P <- rbindlist(pred)
fwrite(P, file.path(ROOT, "results", "theorem2_predictions.csv"))
print(P, digits = 3)
