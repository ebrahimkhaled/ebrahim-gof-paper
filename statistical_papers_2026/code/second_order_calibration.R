## second_order_calibration.R -- calibration of the second-order prediction over ALL 60 settings, without conditioning
## on the simulated outcome (the mock referee's point M3): sign agreement where the PREDICTED gain exceeds 0.01, the
## disagreements, and mean absolute error and rank correlation by design.
suppressMessages(library(data.table))
V <- fread("./results/verify_N1_second_order.csv")
P <- V[abs(second_order) > 0.01]
cat("settings with |predicted| > 0.01:", nrow(P), "; sign agrees:", sum(sign(P$simulated) == sign(P$second_order)), "\n")
print(P[sign(simulated) != sign(second_order), .(design, truth, n, predicted = round(second_order, 3), simulated)])
for (d in c("D1", "D3", "D4")) { x <- V[design == d]
  cat(d, ": MAE", round(mean(abs(x$simulated - x$second_order)), 4), " Spearman", round(cor(x$simulated, x$second_order, method = "spearman"), 2), "\n") }
cat("all 60: MAE", round(mean(abs(V$simulated - V$second_order)), 4), " Spearman", round(cor(V$simulated, V$second_order, method = "spearman"), 2), "\n")
cat("negative predictions:", V[second_order < -0.01, .N], "; of them simulated negative:", V[second_order < -0.01 & simulated < 0, .N], "\n")
