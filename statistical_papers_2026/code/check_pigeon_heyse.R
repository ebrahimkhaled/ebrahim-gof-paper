## check_pigeon_heyse.R -- is Pigeon-Heyse just HL with a within-group variance correction, and why is it conservative?
## J2 = sum (o_g - e_g)^2 / Lambda_g with Lambda_g = sum_{i in g} p_i (1 - p_i), referred to chi-square(G - 1).
## Same first-order theory as Theorem 1: r~ = Lambda^{-1/2} s ~ N(0, Sigma_PH), Sigma_PH = I - Lambda^{-1/2} U M^{-1} U' Lambda^{-1/2},
## so J2 -> sum mu_j chi2_1 with mu_j the eigenvalues of Sigma_PH. Compare the eigenvalue sum with G - 1 and G - 2.
ROOT <- "."
source(file.path(ROOT, "code", "rivals.R"))
src <- readLines(file.path(ROOT, "code", "run_ef_studies.R"))
eval(parse(text = src[grep("^design_data <- function", src):(grep("^## one replicate", src) - 1)]))
set.seed(11)
for (design in c("D1", "D3", "D4", "D5")) {
  d <- design_data(design, 1000); y <- rbinom(1000, 1, plogis(d$eta)); X <- cbind(1, d$X)
  p <- glm.fit(X, y, family = binomial())$fitted.values
  g <- ef_groups(p, 10); v <- p * (1 - p)
  s <- as.numeric(rowsum(y - p, g, reorder = TRUE)); Lam <- as.numeric(rowsum(v, g, reorder = TRUE))
  ng <- tabulate(g, 10); pb <- as.numeric(rowsum(p, g, reorder = TRUE)) / ng; V <- ng * pb * (1 - pb)
  U <- rowsum(v * X, g, reorder = TRUE); M <- crossprod(X, v * X)
  J2 <- sum(s^2 / Lam); HL <- sum(s^2 / V)
  ## the battery's Pigeon-Heyse (quantile groups) for comparison
  fq <- bt_fit(list(d = data.frame(y = y, d$X), f = as.formula(paste("y ~", paste(colnames(d$X), collapse = "+")))))
  pPH <- bt_ph(fq$y, fq$ph, 10)
  SPH <- diag(10) - (U %*% solve(M, t(U))) / sqrt(outer(Lam, Lam))
  SHL <- (diag(Lam) - U %*% solve(M, t(U))) / sqrt(outer(V, V))
  mu <- eigen(SPH, symmetric = TRUE, only.values = TRUE)$values
  la <- eigen(SHL, symmetric = TRUE, only.values = TRUE)$values
  cat(sprintf("%s: J2 %.2f (battery p %.3f, own p %.3f)  HL %.2f | Lambda/V range %.3f-%.3f\n", design, J2, pPH,
              pchisq(J2, 9, lower.tail = FALSE), HL, min(Lam / V), max(Lam / V)))
  cat("   PH eigenvalues:", round(mu, 3), " sum", round(sum(mu), 2), "(reference uses G-1 = 9)\n")
  cat("   HL eigenvalues:", round(la, 3), " sum", round(sum(la), 2), "(reference uses G-2 = 8)\n")
  ## size of J2 under chi2(9), chi2(8) and its own weighted limit, by quick simulation of the Gaussian limit
  z <- matrix(rnorm(2e5 * 10), 2e5) %*% (eigen(SPH, symmetric = TRUE)$vectors %*% diag(sqrt(pmax(mu, 0))))
  Q <- rowSums(z^2)
  cat(sprintf("   limit law of J2: P(> chi2_9 5%% point) = %.4f, P(> chi2_8 5%% point) = %.4f\n",
              mean(Q > qchisq(.95, 9)), mean(Q > qchisq(.95, 8))))
}
