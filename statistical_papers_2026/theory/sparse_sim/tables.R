## tables.R -- markdown tables for SPARSE_GROUP_THEOREM.md from the summary CSVs
suppressMessages(library(data.table))
nul <- fread("out/null_summary.csv"); pw <- fread("out/power_summary.csv")
cb  <- fread("out/combo_summary.csv"); wd <- fread("out/wide_summary.csv")
f2 <- function(x) formatC(x, format = "f", digits = 2); f3 <- function(x) formatC(x, format = "f", digits = 3)
md <- function(dt) {
  h <- paste0("| ", paste(names(dt), collapse = " | "), " |")
  s <- paste0("|", paste(rep("---:", ncol(dt)), collapse = "|"), "|")
  b <- apply(dt, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  paste(c(h, s, b), collapse = "\n")
}
t1 <- nul[oracle == FALSE, .(n, m, G = n / m,
  `sd EF emp/pred` = paste0(f2(sd_EF), " / ", f2(pred_sd_EF)),
  `sd HL emp/pred` = paste0(f2(sd_HL), " / ", f2(pred_sd_HL)),
  `mean EF-G emp/pred` = paste0(f2(mean_EF_minus_G), " / ", f2(pred_mu_EF)),
  `mean C` = f3(mean_C),
  `Z_EF` = f3(size_zEF), `Z_HL` = f3(size_zHL), `abs Z_HL` = f3(size_zHL2),
  `HL chi2(G-2)` = f3(size_HLchisq), `EF chi2(G-2)` = f3(size_EFchisq), `skew Z_EF` = f2(skew_zEF))]
t1o <- nul[oracle == TRUE, .(n, m, `sd EF emp/pred` = paste0(f2(sd_EF), " / ", f2(pred_sd_EF)),
  `sd HL emp/pred` = paste0(f2(sd_HL), " / ", f2(pred_sd_HL)), `Z_EF` = f3(size_zEF), `Z_HL` = f3(size_zHL))]
t2 <- wd[, .(n, m, `R hat` = f2(Rhat), `var(HL)/var(EF) emp` = f3(var_ratio), `plug-in pred` = f3(pred_ratio),
  `1+R/(2(m-1))` = f3(simple_ratio), `Z_EF` = f3(size_zEF), `Z_HL` = f3(size_zHL), `HL chi2` = f3(size_HLchisq))]
p <- merge(pw, cb[, .(scen, n, m, rej_zL2, rej_fisher)], by = c("scen", "n", "m"))
t3 <- p[, .(alt = scen, n, m, Delta = f2(Delta), A = f2(A),
  `EF emp/pred` = paste0(f3(pow_zEF), " / ", f3(pred_EF)),
  `HL emp/pred` = paste0(f3(pow_zHL), " / ", f3(pred_HL)),
  `abs HL emp/pred` = paste0(f3(pow_zHL2), " / ", f3(pred_HL2)),
  `HL chi2` = f3(pow_HLchisq), `Z_L 2-sided` = f3(rej_zL2), `Fisher(Z_EF,Z_L)` = f3(rej_fisher))]
t4 <- cb[scen == "null", .(n, m, `sd C emp/pred` = paste0(f3(sd_C), " / ", f3(pred_sd_C)),
  `cor(C,EF)` = f3(corr_C_EF), `Z_L size` = f3(rej_zL2), `Fisher size` = f3(rej_fisher))]
writeLines(c("### T1", md(t1), "", "### T1o", md(t1o), "", "### T2", md(t2), "", "### T3", md(t3), "",
             "### T4", md(t4)), "out/tables.md")
cat("written out/tables.md\n")
