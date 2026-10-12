## sim.R -- Monte Carlo check of Theorems 1-2 (sparse-group regime).
## Every chunk of replicates is written to its own CSV as soon as it finishes; finished chunks are
## skipped on restart (crash-safe, resumable). Seeds are per chunk, so a rerun reproduces the same data.
library(parallel)
wd <- normalizePath(".")
dir.create("out/raw", showWarnings = FALSE, recursive = TRUE)
ms <- c(2, 5, 10, 25, 50); ns <- c(1000, 4000, 16000)
reps  <- c(null = 4000, cloglog = 1000, loglog = 1000, quadpos = 1000, quadneg = 1000)
chunk <- 50
tasks <- do.call(rbind, lapply(names(reps), function(s) do.call(rbind, lapply(ns, function(n)
  data.frame(scen = s, n = n, k = seq_len(reps[[s]] / chunk))))))
tasks$file <- sprintf("out/raw/%s_%d_%03d.csv", tasks$scen, tasks$n, tasks$k)
tasks <- tasks[!file.exists(tasks$file), ]
cat("chunks to run:", nrow(tasks), "\n")

run_chunk <- function(i) {
  tk <- tasks[i, ]
  seed <- 1e6 * match(tk$scen, names(reps)) + 1e3 * match(tk$n, ns) + tk$k
  set.seed(seed)
  rows <- vector("list", chunk * length(ms) * 2)
  j <- 0
  for (r in seq_len(chunk)) {
    X  <- sg_design(tk$n)
    pt <- sg_truth(X, tk$scen)
    y  <- rbinom(tk$n, 1, pt)
    fit <- glm.fit(X, y, family = binomial())
    p  <- fit$fitted.values
    eta_true <- 0.6 * X[, 2] + 0.5 * X[, 3]
    for (m in ms) {
      for (grp in c("fitted", "oracle")) {
        ord <- if (grp == "fitted") order(p) else order(eta_true)
        j <- j + 1
        rows[[j]] <- c(rep = (tk$k - 1) * chunk + r, m = m, oracle = grp == "oracle",
                       sg_stats(y, p, X, m, ord))
      }
    }
  }
  out <- data.frame(scen = tk$scen, n = tk$n, do.call(rbind, rows[seq_len(j)]))
  write.csv(out, paste0(tk$file, ".tmp"), row.names = FALSE)
  file.rename(paste0(tk$file, ".tmp"), tk$file)
  tk$file
}

cl <- makePSOCKcluster(min(16, detectCores() - 2))
clusterExport(cl, c("tasks", "ms", "ns", "reps", "chunk", "wd"))
clusterEvalQ(cl, { setwd(wd); source("sg_core.R") })
t0 <- Sys.time()
done <- parLapplyLB(cl, seq_len(nrow(tasks)), run_chunk)
stopCluster(cl)
cat("finished", length(done), "chunks in", format(Sys.time() - t0), "\n")
