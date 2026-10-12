## analyse_study2.R -- Study 2: power (raw and size-adjusted) and the Theorem-2 check, E(C) = A(delta).
## Size adjustment uses the matched null of the same design: Study-1 cells for n <= 1000 and the n = 2000
## logit cells run inside Study 2. MCSE of a mean of C is sd(C)/sqrt(B).
ROOT <- "."
suppressMessages(library(data.table))
rd <- function(st, design, truth, n) fread(file.path(ROOT, "results", paste0("study", st),
                                                     sprintf("%s_%s_n%05d_k00_m1.csv.gz", design, truth, n)))
P <- fread(file.path(ROOT, "results", "theorem2_predictions.csv"))
tests <- c(HLchi = "p_HL_chisq", HLx = "p_HL_exact", EFchi = "p_EF_chisq", EFx = "p_EF_exact", Stk = "p_Stukel")
out <- list()
for (design in c("D1", "D3", "D4")) for (truth in c("cloglog", "loglog", "probit", "quad", "inter")) for (n in c(200, 500, 1000, 2000)) {
  A0 <- if (n <= 1000) rd(1, design, "logit", n) else rd(2, design, "logit", n)
  A1 <- rd(2, design, truth, n)
  row <- list(design = design, truth = truth, n = n, meanC = mean(A1$C), seC = sd(A1$C) / sqrt(nrow(A1)),
              A = P$A[P$design == design & P$truth == truth & P$n == n])
  for (t in names(tests)) {
    pv0 <- A0[[tests[t]]]; pv1 <- A1[[tests[t]]]
    crit <- quantile(pv0, .05, na.rm = TRUE, type = 1)       # size-adjusted: reject when p below the null 5% point
    row[[paste0("raw_", t)]] <- mean(pv1 < .05, na.rm = TRUE)
    row[[paste0("adj_", t)]] <- mean(pv1 <= crit, na.rm = TRUE)
  }
  out[[length(out) + 1]] <- as.data.table(row)
}
S <- rbindlist(out)
S[, A_gap_in_se := (meanC - A) / seC]
S[, gain_adj := adj_EFchi - adj_HLchi]
fwrite(S, file.path(ROOT, "results", "study2_summary.csv"))
options(width = 220)
print(S[, .(design, truth, n, A = round(A, 3), meanC = round(meanC, 3), gap_se = round(A_gap_in_se, 1),
            HL = adj_HLchi, EF = adj_EFchi, EFx = adj_EFx, gain = round(gain_adj, 3), Stk = adj_Stk)])
cat("\nsign agreement of mean(C) and A where |A| > 2 MCSE:",
    S[abs(A) > 2 * seC, sum(sign(meanC) == sign(A))], "of", S[abs(A) > 2 * seC, .N], "\n")
