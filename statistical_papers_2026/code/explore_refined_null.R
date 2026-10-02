## explore_refined_null.R -- EXPLORATORY (not part of the hashed specification): does the second-order
## null mean explain EF's conservativeness in sparse designs, and does the refined reference fix it?
## Cells: the Study-1 failures (D2 and D4 at n = 100, 200, 500) and two clean controls (D1 n = 200, D5 n = 200).
## Fresh seeds (9e8 + ...), so this does not reuse Study 1's data.
ROOT <- "."
source(file.path(ROOT, "code", "ef_exact.R"))
suppressMessages({ library(parallel); library(data.table) })
src <- readLines(file.path(ROOT, "code", "run_ef_studies.R"))
eval(parse(text = src[grep("^design_data <- function", src):(grep("^true_p <- function", src) - 1)]))
one <- function(rep, design, n) {
  set.seed(9e8 + rep + 1e4 * n + 1e6 * match(design, c("D1", "D2", "D3", "D4", "D5")))
  d <- design_data(design, n); y <- rbinom(n, 1, plogis(d$eta)); X <- cbind(1, d$X)
  if (min(sum(y), n - sum(y)) < 2) return(NULL)
  f <- suppressWarnings(glm.fit(X, y, family = binomial())); if (!isTRUE(f$converged)) return(NULL)
  c(ef_test(y, f$fitted.values, X)[c("p_EF_chisq", "p_EF_exact", "p_HL_chisq", "p_HL_exact")],
    ef_test_refined(y, f$fitted.values, X))
}
cl <- makePSOCKcluster(20)
clusterExport(cl, c("design_data", "one", "ROOT"))
invisible(clusterEvalQ(cl, source(file.path(ROOT, "code", "ef_exact.R"))))
cells <- rbind(CJ(design = c("D2", "D4"), n = c(100L, 200L, 500L)), data.table(design = c("D1", "D5"), n = 200L))
out <- rbindlist(lapply(seq_len(nrow(cells)), function(j) {
  R <- do.call(rbind, clusterApplyLB(cl, 1:4000, one, design = cells$design[j], n = cells$n[j]))
  data.table(cells[j], B = nrow(R), meanC = mean(R[, "C"]), meanC0 = mean(R[, "C0"]),
             EF_chisq = mean(R[, "p_EF_chisq"] < .05), EF_exact = mean(R[, "p_EF_exact"] < .05),
             EF_refined = mean(R[, "p_EF_refined"] < .05),
             HL_chisq = mean(R[, "p_HL_chisq"] < .05), HL_refined = mean(R[, "p_HL_refined"] < .05))
}))
stopCluster(cl)
options(width = 180); print(out, digits = 3)
fwrite(out, file.path(ROOT, "results", "explore_refined_null.csv"))
