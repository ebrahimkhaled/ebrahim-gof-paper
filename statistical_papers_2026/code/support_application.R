## support_application.R -- Addendum 5: the SUPPORT data, models M0 and M1 (declarations/ADDENDUM5_support.md).
## The data file support2.csv is public (Vanderbilt Biostatistics, also Hmisc::getHdata(support2)); set SUPPORT_CSV
## to its location. Writes results/support_application.csv and results/support_bootstrap.csv.gz.
ROOT <- "."
CSV  <- Sys.getenv("SUPPORT_CSV", "support2.csv")
suppressMessages({ library(data.table); library(splines); library(parallel) })
source(file.path(ROOT, "code", "rivals.R")); source(file.path(ROOT, "code", "large_sample.R"))

s <- read.csv(CSV, na.strings = c("NaN", "NA", ""))
d <- data.frame(y = as.integer(s$hospdead), age = s$age, meanbp = s$meanbp, hrt = s$hrt, resp = s$resp, temp = s$temp,
                crea = s$crea, sod = s$sod, wblc = s$wblc, num.co = s$num.co, sex = factor(s$sex))
d <- na.omit(d); set.seed(20260931L); d <- d[sample(nrow(d)), ]; rownames(d) <- NULL
f0 <- y ~ age + meanbp + hrt + resp + temp + crea + sod + wblc + num.co + sex
f1 <- y ~ age + ns(meanbp, 3) + hrt + resp + temp + crea + sod + wblc + num.co + sex

R_const <- function(p, X) { v <- p * (1 - p); z <- (1 - 2 * p) / v; f <- lm.wfit(X, z, v); sum(v * f$residuals^2) / length(p) }

analyse <- function(form, label) {
  fit <- glm(form, binomial, data = d); X <- model.matrix(fit); p <- fitted(fit); y <- d$y; n <- nrow(d)
  G <- paul_G(n, y); m <- n / G; R <- R_const(p, X)
  fHL <- 1 - 1 / m + R / (2 * m); fEF <- 1 - 1 / m
  size_chi <- function(fac) pnorm(qnorm(.95) / sqrt(fac), lower.tail = FALSE)    # normal approximation of chi2(G-2)
  set.seed(4242); s_ <- ls_stats(y, p, X, G)
  o <- sample.int(n); e10 <- ef_test(y[o], p[o], X[o, , drop = FALSE], 10)
  fq <- bt_fit(list(d = transform(d, y = y), f = form)); stk <- bt_stukel(fq)
  qe <- quantile(fq$eta, c(.025, .975), names = FALSE); fw <- fq; fw$eta <- pmin(pmax(fq$eta, qe[1]), qe[2])
  jw <- bt_stukel_joint_stat(fw)
  data.table(model = label, n = n, events = sum(y), G = G, m = round(m, 1), R = R, factor_EF = fEF, factor_HL = fHL,
             implied_size_chi_EF = size_chi(fEF), implied_size_chi_HL = size_chi(fHL),
             HL = s_[["HL"]], EF = s_[["EF"]], p_HL_chisq = s_[["p_HL_chisq"]], p_HL_norm = s_[["p_HL_norm"]],
             p_EF_chisq = s_[["p_EF_chisq"]], p_EF_norm = s_[["p_EF_norm"]],
             p_HL10 = unname(e10["p_HL_chisq"]), p_EF10 = unname(e10["p_EF_chisq"]),
             p_Stukel = unname(stk["Stk.joint"]), p_StukelW = if (is.finite(jw$chi)) pchisq(jw$chi, jw$df, lower.tail = FALSE) else NA)
}
A <- rbind(analyse(f0, "M0"), analyse(f1, "M1"))

## parametric bootstrap under the fitted M0 (B = 999)
fit0 <- glm(f0, binomial, data = d); p0 <- fitted(fit0)
cl <- makePSOCKcluster(20)
clusterExport(cl, c("d", "f0", "p0", "ROOT"))
invisible(clusterEvalQ(cl, { library(data.table); library(splines); source(file.path(ROOT, "code", "large_sample.R")); NULL }))
BT <- rbindlist(parLapply(cl, 1:999, function(b) {
  set.seed(7e7 + b); db <- d; db$y <- rbinom(nrow(d), 1, p0)
  fb <- suppressWarnings(glm(f0, binomial, data = db)); X <- model.matrix(fb); pb <- fitted(fb)
  s <- ls_stats(db$y, pb, X, paul_G(nrow(db), db$y))
  data.table(b = b, HL = s[["HL"]], EF = s[["EF"]], G = s[["G"]], p_HL_chisq = s[["p_HL_chisq"]], p_HL_norm = s[["p_HL_norm"]],
             p_EF_chisq = s[["p_EF_chisq"]], p_EF_norm = s[["p_EF_norm"]])
}))
stopCluster(cl)
fwrite(BT, file.path(ROOT, "results", "support_bootstrap.csv.gz"))
A[model == "M0", `:=`(p_HL_boot = (1 + sum(BT$HL >= HL)) / 1000, p_EF_boot = (1 + sum(BT$EF >= EF)) / 1000,
                      boot_size_HL_chisq = mean(BT$p_HL_chisq < .05), boot_size_HL_norm = mean(BT$p_HL_norm < .05),
                      boot_size_EF_chisq = mean(BT$p_EF_chisq < .05), boot_size_EF_norm = mean(BT$p_EF_norm < .05))]
fwrite(A, file.path(ROOT, "results", "support_application.csv"))
print(t(A))
