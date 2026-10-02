## check_large_sample.R -- quick checks of large_sample.R before any declared run:
## (1) EF - G equals sum_g U_g / n_g exactly; (2) the score U is uncorrelated with the beta-score (third moments);
## (3) a small null run: sizes of the normal references, L_W and the combined test at n = 5000 with Paul's G.
ROOT <- "."
source(file.path(ROOT, "code", "large_sample.R"))
set.seed(1)
n <- 3000; x1 <- runif(n, -3, 3); x2 <- rbinom(n, 1, .5); X <- cbind(1, x1, x2)
y <- rbinom(n, 1, plogis(0.6 * x1 + 0.5 * x2)); p <- glm.fit(X, y, family = binomial())$fitted.values
g <- ls_groups(p, 120); s <- ls_stats(y, p, X, g = g); sc <- ls_score_rho(y, p, g)
cat(sprintf("EF - G = %.12f | sum U_g/n_g = %.12f\n", s[["EF"]] - s[["G"]], sc[["EF_minus_G"]]))
## orthogonality: cov(U, x' (y - pi)) over replicates at the true risk
pr <- plogis(0.6 * x1 + 0.5 * x2); gg <- ls_groups(pr, 120)
R <- t(replicate(4000, { yy <- rbinom(n, 1, pr); c(ls_score_rho(yy, pr, gg)[["U"]], colSums(X * (yy - pr))) }))
cat("correlations of U with the beta-score components:", round(cor(R)[1, -1], 3), "\n")
## small null run with Paul's G
sim <- replicate(400, {
  yy <- rbinom(5000, 1, plogis(0.6 * runif(5000, -3, 3)))
})
library(parallel)
cl <- makePSOCKcluster(16); clusterExport(cl, "ROOT"); invisible(clusterEvalQ(cl, source(file.path(ROOT, "code", "large_sample.R"))))
S <- do.call(rbind, parLapply(cl, 1:1000, function(r) {
  set.seed(5e6 + r); n <- 5000; x1 <- runif(n, -3, 3); x2 <- rbinom(n, 1, .5); X <- cbind(1, x1, x2)
  y <- rbinom(n, 1, plogis(0.6 * x1 + 0.5 * x2)); p <- glm.fit(X, y, family = binomial())$fitted.values
  ls_stats(y, p, X, paul_G(n, y))
}))
stopCluster(cl)
cat("G =", unique(S[, "G"]), "\n")
print(round(colMeans(S[, c("p_HL_chisq", "p_EF_chisq", "p_HL_norm", "p_EF_norm", "p_L", "p_LW", "p_comb")] < .05), 3))
