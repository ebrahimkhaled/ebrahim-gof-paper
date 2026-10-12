## large_sample.R -- the grouped statistics with many groups (Theorems 3-4 of the paper) for any number of groups G:
## EF and HL with their plug-in normal references, the bounded direction component L_W, and the combined test.
##
##   ls_groups(p, G)          equal-size groups (sizes differ by at most one) by rank of the fitted risk, ties at random
##   paul_G(n, y)             the number of groups of Paul, Pennell and Lemeshow (2013), for 1000 < n <= 25000
##   ls_stats(y, p, X, G)     statistics and p-values; X includes the intercept
##
## The normal-reference moments are those of theory/sparse_sim/sg_core.R (checked against simulation in Online
## Resource 1, Tables S1-S4), written here for a group vector instead of a fixed group size.
## L_W = sum_g zw_g (o_g - e_g) / n_g with zw_g = zeta(pbar_g), zeta(pi) = (1 - 2 pi) / {pi (1 - pi)}, winsorized at
## the 2.5% and 97.5% quantiles of zeta over the records, so that no group carries an unbounded weight. Its null
## variance is the sigma_L^2 of Theorem 3 with these weights: sum_g zw_g^2 tau_g / n_g^2 - t' M^-1 t,
## t = sum_g zw_g U_g / n_g. Without winsorizing it is C / sigma_L, to first order the Osius-Rojek statistic.
## The combined test is Fisher's, -2 {log p_EF + log p_LW} against chi-square(4): the two components are asymptotically
## independent standard normal variables (Theorem 3(iii)).

ls_groups <- function(p, G) {
  n <- length(p); o <- order(p, runif(n))
  g <- integer(n); g[o] <- ceiling(seq_len(n) * G / n)
  g
}

paul_G <- function(n, y) {
  n1 <- sum(y)
  max(10, min(n1 / 2, (n - n1) / 2, 2 + 8 * (n / 1000)^2)) |> floor()
}

ls_stats <- function(y, p, X, G, g = NULL) {
  if (is.null(g)) g <- ls_groups(p, G)
  G   <- max(g)
  ng  <- tabulate(g, G)
  v   <- p * (1 - p)
  o   <- rowsum(y, g, reorder = TRUE)[, 1]
  e   <- rowsum(p, g, reorder = TRUE)[, 1]
  pb  <- e / ng
  V   <- ng * pb * (1 - pb)
  cc  <- 1 - 2 * pb
  HL  <- sum((o - e)^2 / V)
  C   <- sum(cc * (o - e) / V)
  EF  <- HL - C
  tau <- rowsum(v, g, reorder = TRUE)[, 1]
  k3  <- rowsum(v * (1 - 2 * p), g, reorder = TRUE)[, 1]
  k4  <- rowsum(v * (1 - 6 * v), g, reorder = TRUE)[, 1]
  U   <- rowsum(v * X, g, reorder = TRUE)
  K   <- rowsum(v * (1 - 2 * p) * X, g, reorder = TRUE)
  Mi  <- solve(crossprod(X, v * X))
  dtr <- sum(diag(crossprod(U, U / V) %*% Mi))
  mu  <- sum(tau / V) - dtr
  k1  <- colSums((tau / V) * cc * U / V)
  k2  <- colSums(K / V)
  vHL <- sum((2 * tau^2 + k4) / V^2) - 2 * drop(k1 %*% Mi %*% k2) + drop(k1 %*% Mi %*% k1)
  dd  <- k1 - colSums(cc * U / V)
  k3v <- colSums((K - cc * U) / V)
  vEF <- sum((2 * tau^2 + k4 - 2 * cc * k3 + cc^2 * tau) / V^2) - 2 * drop(dd %*% Mi %*% k3v) + drop(dd %*% Mi %*% dd)
  muEF <- sum(tau / V) - dtr * (1 - 1 / mean(ng))
  zEF <- (EF - muEF) / sqrt(vEF); zHL <- (HL - mu) / sqrt(vHL)
  ## direction components: raw (Osius-Rojek type) and winsorized
  zeta <- function(q) (1 - 2 * q) / (q * (1 - q))
  dir_z <- function(zw) {
    L  <- sum(zw * (o - e) / ng)
    tt <- colSums(zw * U / ng)
    s2 <- sum(zw^2 * tau / ng^2) - drop(tt %*% Mi %*% tt)
    if (s2 > 0) L / sqrt(s2) else NA_real_
  }
  qz  <- quantile(zeta(p), c(.025, .975), names = FALSE)
  zL  <- dir_z(zeta(pb))
  zLW <- dir_z(pmin(pmax(zeta(pb), qz[1]), qz[2]))
  pEFn <- pnorm(zEF, lower.tail = FALSE); pLW <- 2 * pnorm(-abs(zLW))
  c(G = G, HL = HL, C = C, EF = EF,
    p_HL_chisq = pchisq(HL, G - 2, lower.tail = FALSE), p_EF_chisq = pchisq(EF, G - 2, lower.tail = FALSE),
    zHL = zHL, zEF = zEF, p_HL_norm = pnorm(zHL, lower.tail = FALSE), p_EF_norm = pEFn,
    zL = zL, zLW = zLW, p_L = 2 * pnorm(-abs(zL)), p_LW = pLW,
    p_comb = pchisq(-2 * (log(pEFn) + log(pLW)), 4, lower.tail = FALSE))
}

## the score statistic for a common within-group correlation (Proposition of Section 4): with the HL table's
## working model (each record of group g has risk pbar_g), U_g = sum_{i != j in g} z_i z_j, z_i = (y_i - pbar_g) /
## sqrt(pbar_g (1 - pbar_g)); the score for rho is U / 2 with U = sum_g U_g; numerically EF - G = sum_g U_g / n_g
ls_score_rho <- function(y, p, g) {
  ng <- tabulate(g); pb <- (rowsum(p, g, reorder = TRUE)[, 1]) / ng
  z  <- (y - pb[g]) / sqrt(pb[g] * (1 - pb[g]))
  s1 <- rowsum(z, g, reorder = TRUE)[, 1]; s2 <- rowsum(z^2, g, reorder = TRUE)[, 1]
  Ug <- s1^2 - s2
  c(U = sum(Ug), EF_minus_G = sum(Ug / ng))
}
