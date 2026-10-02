## combo.R -- the C component alone and a Fisher combination of the two asymptotically independent parts
suppressMessages(library(data.table))
d <- rbindlist(lapply(list.files("out/raw", "[.]csv$", full.names = TRUE), fread))[oracle == FALSE]
d[, sL := sqrt(pmax(sdHL^2 - sdEF^2, 1e-12))]
d[, zL := C / sL]
d[, pQ := pnorm(zEF, lower.tail = FALSE)]
d[, pL := 2 * pnorm(-abs(zL))]
d[, pF := pchisq(-2 * (log(pQ) + log(pL)), 4, lower.tail = FALSE)]
s <- d[, .(sd_C = sd(C), pred_sd_C = sqrt(mean(sL^2)), corr_C_EF = cor(C, EF),
           rej_zL2 = mean(pL < .05), rej_fisher = mean(pF < .05), rej_zEF = mean(pQ < .05),
           rej_zHL2 = mean(abs(zHL) > qnorm(.975))), by = .(scen, n, m)][order(scen, n, m)]
fwrite(s, "out/combo_summary.csv"); options(width = 200); print(s, digits = 3)
