## explore_N3_signfree.R -- EXPLORATORY (N3, not in the hashed spec): the sign-free statistic T = HL + |C|,
## computed from the stored Study-1/Study-2 replicates (same data sets as HL and EF, so all comparisons are paired).
## Size-adjusted power at 5%: each statistic's critical value is the 95% point of its own matched-null values.
ROOT <- "."
suppressMessages(library(data.table))
rd <- function(st, design, truth, n) fread(file.path(ROOT, "results", paste0("study", st),
                                                     sprintf("%s_%s_n%05d_k00_m1.csv.gz", design, truth, n)))
out <- list()
for (design in c("D1", "D3", "D4")) for (n in c(200, 500, 1000, 2000)) {
  N0 <- if (n <= 1000) rd(1, design, "logit", n) else rd(2, design, "logit", n)
  s0 <- list(HL = N0$HL, EF = N0$EF, T = N0$HL + abs(N0$C))
  crit <- lapply(s0, quantile, probs = .95, type = 1)
  for (truth in c("cloglog", "loglog", "probit", "quad", "inter")) {
    A1 <- rd(2, design, truth, n)
    s1 <- list(HL = A1$HL, EF = A1$EF, T = A1$HL + abs(A1$C))
    pw <- mapply(function(s, cr) mean(s > cr), s1, crit)
    out[[length(out) + 1]] <- data.table(design, truth, n, HL = pw["HL"], EF = pw["EF"], T = pw["T"],
                                         EF_gain = pw["EF"] - pw["HL"], T_gain = pw["T"] - pw["HL"])
  }
}
S <- rbindlist(out)
fwrite(S, file.path(ROOT, "results", "explore_N3_signfree.csv"))
options(width = 200); print(S, digits = 3)
cat("\nmean gain over all 60 cells: EF", round(mean(S$EF_gain), 4), " T", round(mean(S$T_gain), 4),
    "\nworst cell: EF", round(min(S$EF_gain), 3), " T", round(min(S$T_gain), 3),
    "\ncells where T beats EF by > 0.01:", sum(S$T - S$EF > .01), " EF beats T by > 0.01:", sum(S$EF - S$T > .01), "\n")
