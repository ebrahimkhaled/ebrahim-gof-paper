## sg_core.R -- grouped HL and grouped Farrington (EF) statistics in the sparse-group regime
## (G = n/m groups of m records each), with the plug-in null moments of Theorem 1.
##
## For groups g of size n_g formed on the fitted risk p:
##   tau_g = sum v_i, k3_g = sum v_i(1-2p_i), k4_g = sum v_i(1-6v_i), U_g = sum v_i x_i, K_g = sum v_i(1-2p_i) x_i
##   V_g = n_g pbar(1-pbar), c_g = 1-2pbar
## mean:   mu    = sum tau/V - tr(U' V^-1 U M^-1)   (HL);  muEF = sum tau/V - (1-1/m) tr(...)  (EF)
## HL var: varHL = sum (2tau^2+k4)/V^2 - 2 k1'M^-1 k2 + k1'M^-1 k1,
##         k1 = sum (tau/V) c U/V,  k2 = sum K/V
## EF var: varEF = sum (2tau^2+k4-2c k3+c^2 tau)/V^2 - 2 d'M^-1 k3v + d'M^-1 d,
##         d = k1 - sum c U/V,  k3v = sum (K - cU)/V
## simple: varEF ~ 2 sum (1 - 1/n_g);  varHL ~ varEF + G R/m (R = weighted RSS of (1-2p)/v on X)

sg_stats <- function(y, p, X, m, ord = NULL) {
  n <- length(y)
  if (is.null(ord)) ord <- order(p)
  g <- integer(n); g[ord] <- (seq_len(n) - 1L) %/% m + 1L
  G <- max(g)
  ## fold a short last group into its neighbour
  ng <- tabulate(g, G)
  if (ng[G] < m) { g[g == G] <- G - 1L; G <- G - 1L; ng <- tabulate(g, G) }
  v  <- p * (1 - p)
  o  <- rowsum(y, g, reorder = TRUE)[, 1]
  e  <- rowsum(p, g, reorder = TRUE)[, 1]
  pb <- e / ng
  V  <- ng * pb * (1 - pb)
  cc <- 1 - 2 * pb
  HL <- sum((o - e)^2 / V)
  C  <- sum(cc * (o - e) / V)
  EF <- HL - C
  tau <- rowsum(v, g, reorder = TRUE)[, 1]
  k3  <- rowsum(v * (1 - 2 * p), g, reorder = TRUE)[, 1]
  k4  <- rowsum(v * (1 - 6 * v), g, reorder = TRUE)[, 1]
  U   <- rowsum(v * X, g, reorder = TRUE)
  K   <- rowsum(v * (1 - 2 * p) * X, g, reorder = TRUE)
  M   <- crossprod(X, v * X)
  Mi  <- solve(M)
  dtr <- sum(diag(crossprod(U, U / V) %*% Mi))
  mu  <- sum(tau / V) - dtr
  k1  <- colSums((tau / V) * cc * U / V)
  k2  <- colSums(K / V)
  vHL <- sum((2 * tau^2 + k4) / V^2) - 2 * drop(k1 %*% Mi %*% k2) + drop(k1 %*% Mi %*% k1)
  dd  <- k1 - colSums(cc * U / V)
  k3v <- colSums((K - cc * U) / V)
  vEF <- sum((2 * tau^2 + k4 - 2 * cc * k3 + cc^2 * tau) / V^2) - 2 * drop(dd %*% Mi %*% k3v) +
         drop(dd %*% Mi %*% dd)
  muEF <- sum(tau / V) - dtr * (1 - 1 / mean(ng))
  c(G = G, HL = HL, C = C, EF = EF, mu = mu, muEF = muEF, d = dtr, sdHL = sqrt(vHL), sdEF = sqrt(vEF),
    sdEF_simple = sqrt(2 * sum(1 - 1 / ng)),
    zHL = (HL - mu) / sqrt(vHL), zEF = (EF - muEF) / sqrt(vEF),
    pHL_chisq = pchisq(HL, G - 2, lower.tail = FALSE))
}

## R: weighted residual mean square of zeta = (1-2p)/v on X with weights v
sg_R <- function(p, X) {
  v <- p * (1 - p); z <- (1 - 2 * p) / v
  f <- lm.wfit(X, z, v)
  sum(v * f$residuals^2) / length(p)
}

## data-generating mechanisms: eta = 0.6 x1 + 0.5 x2
sg_design <- function(n) {
  x1 <- runif(n, -3, 3); x2 <- rbinom(n, 1, 0.5)
  cbind(1, x1, x2)
}
sg_truth <- function(X, scen, gam = 0.12) {
  eta <- 0.6 * X[, 2] + 0.5 * X[, 3]
  switch(scen,
    null    = plogis(eta),
    cloglog = 1 - exp(-exp(eta)),
    loglog  = exp(-exp(-eta)),
    quadpos = plogis(eta + gam * (X[, 2]^2 - 3)),
    quadneg = plogis(eta - gam * (X[, 2]^2 - 3)))
}
