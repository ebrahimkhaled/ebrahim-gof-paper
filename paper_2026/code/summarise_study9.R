## summarise_study9.R -- the numbers the text quotes from Study 9 (event rate held fixed) and the timing of Study 10
ROOT <- "."
suppressMessages(library(data.table))
S <- fread(file.path(ROOT, "results", "study9_summary.csv")); A <- S[truth != "logit"]; A[, g := round(adj_EF - adj_HL, 3)]
v <- A[abs(g) > .010]
cat("cells with |EF-HL| > 0.01:", nrow(v), " A<0:", v[A < 0, .N], " A>0:", v[A > 0, .N], "\n")
cat("sign misses:\n"); print(v[sign(g) != -sign(A), .(design, truth, n, A = round(A, 2), adj_EF, adj_HL, g)])
cat("A > 0 cells:\n"); print(A[A > 0, .(design, truth, n, A = round(A, 2), adj_EF, adj_HL, adj_PH, g)])
cat("largest gain", A[, max(g)], "at\n"); print(A[g == max(g), .(design, truth, n)])
cat("largest loss", A[, min(g)], "at\n"); print(A[g == min(g), .(design, truth, n)])
cat("mean gain", round(A[, mean(g)], 4), "\n")
L <- S[truth == "logit"]
for (t in c("EF", "HL", "PH", "StukelW", "Stukel", "Tsiatis")) cat("null raw", t, range(L[[paste0("raw_", t)]]), "\n")
cat("PH - HL: adjusted", round(A[, mean(adj_PH - adj_HL)], 4), " raw", round(A[, mean(raw_PH - raw_HL)], 4), "\n")
cat("Stukel-W vs EF at n=1000, wrong links:\n"); print(A[n == 1000 & truth %in% c("cloglog", "loglog", "probit"), .(design, truth, EF = adj_EF, StkW = adj_StukelW, Stk = adj_Stukel)])
cat("se of paired gain: range", range(A$se_gain), "\n")
print(fread(file.path(ROOT, "results", "alignment_calibrated.csv"))[, .(design, truth, A = round(A, 2), ev = round(event_rate, 3))])
