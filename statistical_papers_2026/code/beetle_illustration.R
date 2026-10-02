## beetle_illustration.R -- Bliss (1935) flour-beetle mortality (VGAM::fbeetle), the textbook asymmetric-link data.
## The data have eight doses, so the natural grouping is one group per dose (each covariate pattern kept whole; the
## numbers of beetles per dose differ). Reported for the logistic fit: HL, C and EF with one group per dose, the
## read-out c_g = b_g r_g, and the alignment A obtained by taking the cloglog fit as the truth; then the same for the
## cloglog fit. A rank grouping with G = 10 and random tie-breaking is shown as a check that the order of tied
## records does not drive the result (the fixed "events first" order manufactures misfit, see the paper, Section 2.1).
ROOT <- "."
source(file.path(ROOT, "code", "ef_exact.R"))
data("fbeetle", package = "VGAM")
d <- fbeetle
i <- rep(seq_len(nrow(d)), d$n)
y <- unlist(lapply(seq_len(nrow(d)), function(k) rep(c(1, 0), c(d$dead[k], d$n[k] - d$dead[k]))))
x <- d$logdose[i]; dose <- i
fl <- glm(y ~ x, binomial("logit")); fc <- glm(y ~ x, binomial("cloglog"))
X <- model.matrix(fl)
cat("n =", length(y), " deaths =", sum(y), " mean fitted risk (logit) =", round(mean(fitted(fl)), 3),
    "\nAIC logit", round(AIC(fl), 2), " cloglog", round(AIC(fc), 2), "\n\n")
t <- ef_test(y, fitted(fl), X, groups = dose); z <- ef_parts(y, fitted(fl), X, groups = dose)
ng <- tabulate(dose)
dg <- as.numeric(rowsum(fitted(fc) - fitted(fl), dose, reorder = TRUE)) / ng
A  <- sum((1 - 2 * z$pbar) * dg / (z$pbar * (1 - z$pbar)))
cat(sprintf("one group per dose (G = 8): HL %.2f (p %.3f)  C %.2f  EF %.2f (p %.3f)  predicted A %.2f\n",
            t["HL"], t["p_HL_chisq"], t["C"], t["EF"], t["p_EF_chisq"], A))
print(round(data.frame(dose = round(d$logdose, 4), n = ng, pbar = z$pbar, r = z$r, b = z$b, c_g = z$b * z$r), 3))
tc <- ef_test(y, fitted(fc), model.matrix(fc), groups = dose)
cat(sprintf("\ncloglog fit, one group per dose: HL %.2f (p %.3f)  C %.2f  EF %.2f (p %.3f)\n",
            tc["HL"], tc["p_HL_chisq"], tc["C"], tc["EF"], tc["p_EF_chisq"]))
## rank grouping, G = 10, random tie-breaking: the distribution over 200 random orders
set.seed(20261002)
R <- t(replicate(200, { o <- sample.int(length(y)); ef_test(y[o], fitted(fl)[o], X[o, ], 10)[c("p_HL_chisq", "p_EF_chisq")] }))
cat(sprintf("\nrank groups G = 10, random tie order (200 orders): median p HL %.3f, EF %.3f; range HL %.3f-%.3f\n",
            median(R[, 1]), median(R[, 2]), min(R[, 1]), max(R[, 1])))
saveRDS(list(z = z, t = t, A = A, tc = tc, ng = ng, logdose = d$logdose), file.path(ROOT, "results", "beetle.rds"))
