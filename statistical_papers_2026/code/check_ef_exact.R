## check_ef_exact.R -- quick sanity checks of ef_exact.R before any study is run:
## (1) the statistic agrees with ebrahim.gof's ef.gof on the same data and groups;
## (2) under the null the simulated mean of EF matches the theoretical mean sum(lam) - 0;
## (3) the exact p-value is uniform-looking on a small run.
source("./code/ef_exact.R")
set.seed(20261002)
gen <- function(n) {
  x1 <- runif(n, -3, 3); x2 <- rbinom(n, 1, .5)
  y  <- rbinom(n, 1, plogis(0.6 * x1 + 0.5 * x2))
  data.frame(y, x1, x2)
}
d   <- gen(1000)
fit <- glm(y ~ x1 + x2, binomial, data = d)
X   <- model.matrix(fit); p <- fitted(fit)
out <- ef_test(d$y, p, X)
print(round(out, 4))
pk <- tryCatch(ebrahim.gof::ef.gof(d$y, p, G = 10), error = function(e) conditionMessage(e))
print(pk)

R <- t(replicate(400, {
  d <- gen(500); f <- glm(y ~ x1 + x2, binomial, data = d)
  ef_test(d$y, fitted(f), model.matrix(f))
}))
cat("\nmean EF", round(mean(R[, "EF"]), 3), " mean sum(lam)", round(mean(R[, "sum_lam"]), 3),
    " mean HL", round(mean(R[, "HL"]), 3), "\n")
cat("size at 5%: HL chisq", mean(R[, "p_HL_chisq"] < .05), " HL exact", mean(R[, "p_HL_exact"] < .05),
    " EF chisq", mean(R[, "p_EF_chisq"] < .05), " EF exact", mean(R[, "p_EF_exact"] < .05), "\n")
cat("mean kappa", signif(mean(R[, "kappa"]), 3), " mean C", round(mean(R[, "C"]), 4), "\n")
