## ef_exact.R -- the grouped Farrington (EF) test with its exact limiting null distribution.
##
## EF = HL - C on G equal-size risk groups, HL = sum r_g^2, C = sum b_g r_g, with
##   r_g = (o_g - e_g) / sqrt(V_g),  V_g = n_g pbar_g (1 - pbar_g),  b_g = (1 - 2 pbar_g) / sqrt(V_g).
## To first order the residual vector r is N(m, Sigma) with
##   Sigma = V^{-1/2} (Lambda - U M^{-1} U') V^{-1/2},
## Lambda = diag(sum_{i in g} v_i), U = rows sum_{i in g} v_i x_i', M = X'WX, v_i = p_i (1 - p_i).
## Writing Sigma = P diag(lam) P' and beta = P'b, under the null
##   EF = sum_j lam_j (xi_j - beta_j / (2 sqrt(lam_j)))^2 - kappa,   kappa = sum_j beta_j^2 / 4,
## so P(EF > q) = P( sum_j lam_j chi2_1(delta_j) > q + kappa ),  delta_j = beta_j^2 / (4 lam_j),
## computed by Imhof's method (Theorem 1 of the paper).
suppressMessages(requireNamespace("CompQuadForm"))

## equal-size groups by rank of the fitted risk, sizes differing by at most one, ties in record order
ef_groups <- function(p, G) {
  o <- order(p, seq_along(p))
  g <- integer(length(p))
  g[o] <- cut(seq_along(p), breaks = G, labels = FALSE)
  g
}

## the pieces every statistic below needs: grouped residuals, their first-order covariance, the weights
ef_parts <- function(y, p, X, G = 10, groups = NULL) {
  ## groups: an optional vector of group labels 1..G (e.g. one group per covariate pattern); default: rank groups
  g   <- if (is.null(groups)) ef_groups(p, G) else as.integer(factor(groups))
  G   <- max(g)
  v   <- p * (1 - p)
  ng  <- tabulate(g, G)
  o   <- as.numeric(rowsum(y, g, reorder = TRUE))
  e   <- as.numeric(rowsum(p, g, reorder = TRUE))
  pb  <- e / ng
  V   <- ng * pb * (1 - pb)
  Lam <- as.numeric(rowsum(v, g, reorder = TRUE))
  U   <- rowsum(v * X, g, reorder = TRUE)
  M   <- crossprod(X, v * X)
  S   <- (diag(Lam, G) - U %*% solve(M, t(U))) / sqrt(outer(V, V))
  r   <- (o - e) / sqrt(V)
  list(r = r, b = (1 - 2 * pb) / sqrt(V), Sigma = (S + t(S)) / 2, pbar = pb, V = V, o = o, e = e)
}

## upper tail of  sum lam_j chi2_1(delta_j) - kappa  at q (Imhof; Davies as a fallback)
qf_upper <- function(q, lam, delta = rep(0, length(lam)), kappa = 0) {
  x <- q + kappa
  if (x <= 0) return(1)
  p <- CompQuadForm::imhof(x, lambda = lam, h = rep(1, length(lam)), delta = delta,
                           epsabs = 1e-8, epsrel = 1e-8, limit = 1e4)$Qq
  if (!is.finite(p) || p < 0 || p > 1)
    p <- CompQuadForm::davies(x, lambda = lam, h = rep(1, length(lam)), delta = delta, acc = 1e-8, lim = 1e5)$Qq
  min(max(p, 0), 1)
}

## EF and HL, each with the classical chi-square(G - 2) reference and with its exact limiting law
ef_test <- function(y, p, X, G = 10, groups = NULL) {
  z    <- ef_parts(y, p, X, G, groups)
  G    <- length(z$r)
  HL   <- sum(z$r^2)
  C    <- sum(z$b * z$r)
  EF   <- HL - C
  es   <- eigen(z$Sigma, symmetric = TRUE)
  keep <- es$values > 1e-10 * max(es$values)
  lam  <- es$values[keep]
  beta <- as.numeric(crossprod(es$vectors[, keep, drop = FALSE], z$b))
  kap  <- sum(beta^2) / 4
  dlt  <- beta^2 / (4 * lam)
  c(HL = HL, C = C, EF = EF, kappa = kap,
    p_HL_chisq = stats::pchisq(HL, G - 2, lower.tail = FALSE),
    p_HL_exact = qf_upper(HL, lam),
    p_EF_chisq = stats::pchisq(EF, G - 2, lower.tail = FALSE),
    p_EF_exact = qf_upper(EF, lam, dlt, kap),
    sum_lam = sum(lam))
}

## Second-order null mean of the grouped residuals. The maximum-likelihood fit is biased away from zero
## by b = M^{-1} sum_i x_i h_i (p_i - 1/2) (Firth 1993; h_i = v_i x_i' M^{-1} x_i), and the fitted risk has a
## curvature term, so E(p_hat_i - p_i) = v_i x_i' b + (1/2) v_i (1 - 2 p_i) x_i' M^{-1} x_i to order 1/n.
## Summed over a group this is the O(1) mean of o_g - e_g under a correct model; the sum over all groups is
## exactly zero (the intercept's score equation), so the bias only moves residuals between groups.
ef_null_mean <- function(p, X, g, G) {
  v  <- p * (1 - p)
  M  <- crossprod(X, v * X)
  Q  <- rowSums((X %*% solve(M)) * X)             # x_i' M^{-1} x_i
  h  <- v * Q
  bb <- solve(M, crossprod(X, h * (p - 0.5)))
  d  <- v * as.numeric(X %*% bb) + 0.5 * v * (1 - 2 * p) * Q
  -as.numeric(rowsum(d, g, reorder = TRUE))       # E(o_g - e_g)
}

## EF with the refined reference: the Theorem-2 law evaluated at the second-order null mean m0
ef_test_refined <- function(y, p, X, G = 10) {
  z    <- ef_parts(y, p, X, G)
  g    <- ef_groups(p, G)
  m0   <- ef_null_mean(p, X, g, G) / sqrt(z$V)
  EF   <- sum(z$r^2) - sum(z$b * z$r)
  es   <- eigen(z$Sigma, symmetric = TRUE)
  keep <- es$values > 1e-10 * max(es$values)
  lam  <- es$values[keep]
  P    <- es$vectors
  beta <- as.numeric(crossprod(P, z$b)); eta <- as.numeric(crossprod(P, m0))
  ## EF = sum_{lam>0} lam_j chi2_1((2 eta_j - beta_j)^2 / (4 lam_j)) - sum beta_j^2/4 + sum_{lam=0}(eta_j^2 - beta_j eta_j)
  kap  <- sum(beta[keep]^2) / 4 - sum(eta[!keep]^2 - beta[!keep] * eta[!keep])
  dlt  <- (2 * eta[keep] - beta[keep])^2 / (4 * lam)
  ## HL refined the same way (b = 0)
  c(EF = EF, C = sum(z$b * z$r), C0 = sum(z$b * m0),
    p_EF_refined = qf_upper(EF, lam, dlt, kap),
    p_HL_refined = qf_upper(sum(z$r^2), lam, eta[keep]^2 / lam, -sum(eta[!keep]^2)))
}

## Stukel's (1988) joint score test for the two shape directions, chi-square(2)
stukel_joint <- function(y, p, X, eta) {
  Z <- cbind(0.5 * eta^2 * (p >= 0.5), -0.5 * eta^2 * (p < 0.5))
  Z <- Z[, colSums(Z != 0) > 0, drop = FALSE]
  v <- p * (1 - p)
  u <- crossprod(Z, y - p)
  A <- crossprod(Z, v * Z) - crossprod(Z, v * X) %*% solve(crossprod(X, v * X), crossprod(X, v * Z))
  s <- tryCatch(as.numeric(crossprod(u, solve(A, u))), error = function(e) NA_real_)
  c(Stukel = s, p_Stukel = stats::pchisq(s, ncol(Z), lower.tail = FALSE))
}
