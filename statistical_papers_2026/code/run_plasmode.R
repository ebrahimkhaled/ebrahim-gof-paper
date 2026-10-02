## run_plasmode.R -- declarations/ADDENDUM1_plasmode.md (hash 42f8c0cd).
## Real covariates of the 1000 burn-injury patients; outcomes simulated from a logit, loglog-type or cloglog-type
## truth with the observed 15% mean risk; the textbook linear logit model refitted to every simulated data set.
ROOT <- "."
source(file.path(ROOT, "code", "ef_exact.R"))
suppressMessages({ library(parallel); library(data.table); library(aplore3) })
OUT <- file.path(ROOT, "results", "plasmode"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

d0  <- transform(burn1000, y = as.integer(death == "Dead"))
f0  <- glm(y ~ age + tbsa + race + inh_inj + flame, binomial, data = d0)
X0  <- model.matrix(f0); eta <- as.numeric(X0 %*% coef(f0)) - coef(f0)[1]   # risk score without the intercept
rate <- mean(d0$y)
link <- list(logit   = function(e) plogis(e),
             loglog  = function(e) exp(-exp(-e)),
             cloglog = function(e) 1 - exp(-exp(e)))
## intercept giving the observed mean risk under each truth
icpt <- sapply(link, function(F) uniroot(function(a) mean(F(a + eta)) - rate, c(-20, 20))$root)
ptrue <- lapply(names(link), function(k) link[[k]](icpt[k] + eta)); names(ptrue) <- names(link)

## the prediction: A from the pseudo-true logit fit to the true risks on these 1000 rows
predA <- sapply(names(link), function(k) {
  fp <- suppressWarnings(glm.fit(X0, ptrue[[k]], family = quasibinomial()))
  ps <- fp$fitted.values; g <- ef_groups(ps, 10); ng <- tabulate(g, 10)
  dg <- as.numeric(rowsum(ptrue[[k]] - ps, g, reorder = TRUE)) / ng
  pb <- as.numeric(rowsum(ps, g, reorder = TRUE)) / ng
  sum((1 - 2 * pb) * dg / (pb * (1 - pb)))
})
cat("observed death rate", round(rate, 3), "\npredicted A:", round(predA, 3), "\n")

one <- function(rep, truth, k) {
  set.seed(5e8 + 1e5 * match(truth, c("logit", "loglog", "cloglog")) + 1e4 * k + rep)
  o <- sample.int(nrow(X0))                               # random order = random tie-breaking in the grouping
  X <- X0[o, ]; y <- rbinom(nrow(X), 1, ptrue[[truth]][o])
  if (k > 0) { bad <- sample.int(nrow(X), k); X[bad, "tbsa"] <- 4 * X[bad, "tbsa"] }
  f <- suppressWarnings(glm.fit(X, y, family = binomial()))
  if (!isTRUE(f$converged)) return(NULL)
  c(rep = rep, ef_test(y, f$fitted.values, X), stukel_joint(y, f$fitted.values, X, f$linear.predictors))
}
cl <- makePSOCKcluster(20)
clusterExport(cl, c("one", "X0", "ptrue", "ROOT"))
invisible(clusterEvalQ(cl, source(file.path(ROOT, "code", "ef_exact.R"))))
cells <- rbind(data.table(truth = c("logit", "loglog", "cloglog"), k = 0L),
               data.table(truth = "logit", k = c(1L, 2L, 5L, 10L)))
res <- list()
for (j in seq_len(nrow(cells))) {
  tr <- cells$truth[j]; k <- cells$k[j]
  f <- file.path(OUT, sprintf("%s_k%02d.csv.gz", tr, k))
  R <- if (file.exists(f)) fread(f) else {
    R <- rbindlist(lapply(clusterApplyLB(cl, 1:2000, one, truth = tr, k = k),
                          function(v) if (is.null(v)) NULL else as.data.table(as.list(v))))
    fwrite(R, f); R }
  res[[j]] <- R[, .(truth = tr, k = k, B = .N, HL = mean(p_HL_chisq < .05), EF = mean(p_EF_chisq < .05),
                    Stukel = mean(p_Stukel < .05, na.rm = TRUE), meanC = mean(C))]
}
stopCluster(cl)
S <- rbindlist(res)
## size-adjusted power for the two alternatives, against the logit-truth null
N0 <- fread(file.path(OUT, "logit_k00.csv.gz"))
for (tr in c("loglog", "cloglog")) {
  R <- fread(file.path(OUT, sprintf("%s_k00.csv.gz", tr)))
  adj <- sapply(c(HL = "HL", EF = "EF"), function(s) mean(R[[s]] > quantile(N0[[s]], .95, type = 1)))
  adjS <- mean(R$Stukel > quantile(N0$Stukel, .95, type = 1, na.rm = TRUE), na.rm = TRUE)
  cat(sprintf("%-8s predicted A %.3f | size-adjusted power: HL %.3f  EF %.3f  (EF - HL %+.3f)  Stukel %.3f\n",
              tr, predA[tr], adj["HL"], adj["EF"], adj["EF"] - adj["HL"], adjS))
}
print(S, digits = 3)
fwrite(S, file.path(ROOT, "results", "plasmode_summary.csv"))
