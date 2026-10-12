## probit_alignment.R -- the alignment A of a probit truth misfitted by the logit (Section 3.3 of the paper):
## eta ~ N(0, s^2), y ~ Bernoulli(Phi(eta)), pseudo-true logistic fit of y on eta from 10^6 records, ten groups by
## fitted risk, A = sum_g (1 - 2 pbar_g) d_g / {pbar_g (1 - pbar_g)}, d_g the mean bias of the fitted risks in group g.
set.seed(20261012)
N <- 1e6
for (s in c(1, 2, 3)) {
  eta <- rnorm(N, 0, s)
  p_true <- pnorm(eta)
  # pseudo-true fit: the logistic score equations solved at E(y) = p_true (no response noise)
  fit <- suppressWarnings(glm(p_true ~ eta, family = quasibinomial))
  p_fit <- fitted(fit)
  g <- cut(rank(p_fit, ties.method = "random"), breaks = 10, labels = FALSE)
  pbar <- tapply(p_fit, g, mean)
  d <- tapply(p_true - p_fit, g, mean)
  A <- sum((1 - 2 * pbar) * d / (pbar * (1 - pbar)))
  cat(sprintf("s = %d: A = %.2f\n", s, A))
}
