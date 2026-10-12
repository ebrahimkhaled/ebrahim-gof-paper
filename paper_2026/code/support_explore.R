## support_explore.R -- EXPLORATORY (beyond Addendum 5): on the SUPPORT data, how the p-values of HL and EF for M0
## depend on the number of groups, and the alignment A and pair drift of the M0 misfit when M1's fitted risks are taken
## as the truth (the bias d_g = mean of p1 - p0 in each group of M0's fitted risks).
ROOT <- "."
CSV  <- Sys.getenv("SUPPORT_CSV", "support2.csv")
suppressMessages({ library(data.table); library(splines) })
source(file.path(ROOT, "code", "rivals.R")); source(file.path(ROOT, "code", "large_sample.R"))
s <- read.csv(CSV, na.strings = c("NaN", "NA", ""))
d <- data.frame(y = as.integer(s$hospdead), age = s$age, meanbp = s$meanbp, hrt = s$hrt, resp = s$resp, temp = s$temp,
                crea = s$crea, sod = s$sod, wblc = s$wblc, num.co = s$num.co, sex = factor(s$sex))
d <- na.omit(d); set.seed(20260931L); d <- d[sample(nrow(d)), ]; rownames(d) <- NULL
fit0 <- glm(y ~ age + meanbp + hrt + resp + temp + crea + sod + wblc + num.co + sex, binomial, data = d)
fit1 <- glm(y ~ age + ns(meanbp, 3) + hrt + resp + temp + crea + sod + wblc + num.co + sex, binomial, data = d)
p0 <- fitted(fit0); p1 <- fitted(fit1); X <- model.matrix(fit0); y <- d$y; n <- length(y)
out <- rbindlist(lapply(c(10, 20, 50, 100, 200, 355, 631), function(G) {
  set.seed(4242); g <- ls_groups(p0, G); s_ <- ls_stats(y, p0, X, g = g)
  ng <- tabulate(g, G); pb <- rowsum(p0, g, reorder = TRUE)[, 1] / ng; V <- ng * pb * (1 - pb)
  h <- p1 - p0; delta <- rowsum(h, g, reorder = TRUE)[, 1]
  A <- sum((1 - 2 * pb) * delta / V)
  Delta <- sum((delta^2 - rowsum(h^2, g, reorder = TRUE)[, 1]) / V)
  data.table(G = G, m = round(n / G, 1), p_HL_chisq = s_[["p_HL_chisq"]], p_EF_chisq = s_[["p_EF_chisq"]],
             p_HL_norm = s_[["p_HL_norm"]], p_EF_norm = s_[["p_EF_norm"]], A = A, Delta = Delta,
             noncentrality_HL = sum(delta^2 / V))
}))
fwrite(out, file.path(ROOT, "results", "support_explore.csv")); print(out, digits = 3)
