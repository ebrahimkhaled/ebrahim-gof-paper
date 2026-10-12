## analyse.R -- summarise out/raw/*.csv into the tables of SPARSE_GROUP_THEOREM.md
suppressMessages(library(data.table))
f  <- list.files("out/raw", pattern = "\\.csv$", full.names = TRUE)
d  <- rbindlist(lapply(f, fread))
d[, oracle := as.logical(oracle)]
d[, pEF_chisq := pchisq(EF, G - 2, lower.tail = FALSE)]
z  <- qnorm(.95); z2 <- qnorm(.975)

## ---- null: moments and size ----
nul <- d[scen == "null", .(R = .N,
  mean_HL_minus_G = mean(HL - G), pred_mu_HL = mean(mu - G),
  mean_EF_minus_G = mean(EF - G), pred_mu_EF = mean(muEF - G),
  mean_C = mean(C), d = mean(d),
  sd_HL = sd(HL), pred_sd_HL = sqrt(mean(sdHL^2)),
  sd_EF = sd(EF), pred_sd_EF = sqrt(mean(sdEF^2)), simple_sd_EF = sqrt(mean(sdEF_simple^2)),
  var_ratio = var(HL) / var(EF),
  size_zEF = mean(zEF > z), size_zHL = mean(zHL > z), size_zHL2 = mean(abs(zHL) > z2),
  size_HLchisq = mean(pHL_chisq < .05), size_EFchisq = mean(pEF_chisq < .05),
  skew_zEF = mean(((zEF - mean(zEF)) / sd(zEF))^3), skew_zHL = mean(((zHL - mean(zHL)) / sd(zHL))^3)),
  by = .(oracle, n, m)][order(oracle, n, m)]
fwrite(nul, "out/null_summary.csv")

## ---- alternatives: power ----
pw <- d[scen != "null" & oracle == FALSE, .(R = .N,
  pow_zEF = mean(zEF > z), pow_zHL = mean(zHL > z), pow_zHL2 = mean(abs(zHL) > z2),
  pow_HLchisq = mean(pHL_chisq < .05),
  shift_EF = mean(EF - muEF), shift_HL = mean(HL - mu), mean_C = mean(C)), by = .(scen, n, m)]
pr <- fread("out/pop_predictions.csv")
pw <- merge(pw, pr[, .(scen, n, m, Delta, A, thetaEF, thetaHL, pred_EF = powEF, pred_HL = powHL, pred_HL2 = powHL2)],
            by = c("scen", "n", "m"))[order(scen, n, m)]
fwrite(pw, "out/power_summary.csv")
options(width = 250)
print(nul, digits = 3); print(pw, digits = 3)
