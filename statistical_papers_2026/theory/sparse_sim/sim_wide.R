## sim_wide.R -- null check of the HL/EF variance gap on a design with a wide risk range
## (eta = 1.2 x1 + 0.5 x2, x1 ~ U(-3,3)), where R = weighted RSS of (1-2p)/v on X is large.
## Chunked, crash-safe and resumable like sim.R.
library(parallel)
wd <- normalizePath(".")
dir.create("out/raw_wide", showWarnings = FALSE, recursive = TRUE)
ms <- c(2, 5, 10, 25, 50); ns <- c(4000); chunk <- 50; nrep <- 2000
tasks <- expand.grid(n = ns, k = seq_len(nrep / chunk))
tasks$file <- sprintf("out/raw_wide/wide_%d_%03d.csv", tasks$n, tasks$k)
tasks <- tasks[!file.exists(tasks$file), ]
run_chunk <- function(i) {
  tk <- tasks[i, ]; set.seed(9e6 + tk$k)
  rows <- list()
  for (r in seq_len(chunk)) {
    X <- sg_design(tk$n); p0 <- plogis(drop(X %*% c(0, 1.2, 0.5)))
    y <- rbinom(tk$n, 1, p0)
    p <- glm.fit(X, y, family = binomial())$fitted.values
    for (m in ms) rows[[length(rows) + 1]] <- c(rep = (tk$k - 1) * chunk + r, m = m, Rhat = sg_R(p, X),
                                                sg_stats(y, p, X, m))
  }
  out <- data.frame(n = tk$n, do.call(rbind, rows))
  write.csv(out, paste0(tk$file, ".tmp"), row.names = FALSE); file.rename(paste0(tk$file, ".tmp"), tk$file)
}
cl <- makePSOCKcluster(min(16, detectCores() - 2))
clusterExport(cl, c("tasks", "ms", "chunk", "wd"))
clusterEvalQ(cl, { setwd(wd); source("sg_core.R") })
invisible(parLapplyLB(cl, seq_len(nrow(tasks)), run_chunk))
stopCluster(cl)
suppressMessages(library(data.table))
d <- rbindlist(lapply(list.files("out/raw_wide", "\\.csv$", full.names = TRUE), fread))
s <- d[, .(R = .N, Rhat = mean(Rhat), sd_HL = sd(HL), pred_sd_HL = sqrt(mean(sdHL^2)), sd_EF = sd(EF),
           pred_sd_EF = sqrt(mean(sdEF^2)), var_ratio = var(HL) / var(EF),
           pred_ratio = mean(sdHL^2 / sdEF^2), simple_ratio = mean(1 + Rhat / (2 * (m - 1))),
           mean_HL_minus_G = mean(HL - G), pred_mu_HL = mean(mu - G), mean_EF_minus_G = mean(EF - G),
           pred_mu_EF = mean(muEF - G), size_zEF = mean(zEF > qnorm(.95)), size_zHL = mean(zHL > qnorm(.95)),
           size_zHL2 = mean(abs(zHL) > qnorm(.975)), size_HLchisq = mean(pHL_chisq < .05)), by = .(n, m)]
fwrite(s, "out/wide_summary.csv"); options(width = 250); print(s, digits = 3)
