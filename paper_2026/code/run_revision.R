## run_revision.R -- Studies 4-7 of declarations/ADDENDUM2_rivals_revision.md (hash c9de6124).
##   Rscript run_revision.R 4|5|6|7
## One file per cell in results/study<k>/, one row per replicate; a finished cell is skipped on restart. Each row
## carries a fingerprint of its data (number of events and the sum of the first covariate), so a replicate that
## silently repeats another (an RNG reset inside a test) is caught by checking that fingerprints are distinct.
ROOT <- "."
suppressMessages({ library(parallel); library(data.table) })
study <- as.integer(commandArgs(TRUE)[1])
OUT <- file.path(ROOT, "results", paste0("study", study)); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
src <- readLines(file.path(ROOT, "code", "run_ef_studies.R"))
gen_src <- src[grep("^design_data <- function", src):(grep("^## one replicate", src) - 1)]

ao <- function(eta, lam, side) {
  if (side == "cloglog") { if (lam == 0) 1 - exp(-exp(eta)) else 1 - (1 + lam * exp(eta))^(-1 / lam) }
  else                   { if (lam == 0) exp(-exp(-eta))    else (1 + lam * exp(-eta))^(-1 / lam) }
}

## burn-injury covariates for Study 7 (Addendum 1 construction)
burn <- local({
  suppressMessages(library(aplore3))
  d0 <- transform(burn1000, y = as.integer(death == "Dead"))
  f0 <- glm(y ~ age + tbsa + race + inh_inj + flame, binomial, data = d0)
  X0 <- model.matrix(f0); colnames(X0) <- make.names(colnames(X0))
  eta <- as.numeric(X0 %*% coef(f0)) - coef(f0)[1]
  link <- list(logit = plogis, loglog = function(e) exp(-exp(-e)), cloglog = function(e) 1 - exp(-exp(e)))
  icpt <- sapply(link, function(F) uniroot(function(a) mean(F(a + eta)) - mean(d0$y), c(-20, 20))$root)
  list(X = X0, p = lapply(names(link), function(k) link[[k]](icpt[k] + eta)) |> setNames(names(link)))
})

one <- function(rep, cell) {
  set.seed(cell$seed + rep)
  if (cell$kind == "design") {
    d <- design_data(cell$design, cell$n); y <- rbinom(cell$n, 1, true_p(cell$truth, d))
    X <- cbind(`(Intercept)` = 1, d$X); colnames(X)[-1] <- paste0("x", seq_len(ncol(X) - 1))
  } else if (cell$kind == "map") {
    x <- rnorm(cell$n); y <- rbinom(cell$n, 1, ao(cell$mu + cell$s * x, cell$lam, cell$side))
    X <- cbind(`(Intercept)` = 1, x1 = x)
  } else {                                              # burn plasmode
    X <- burn$X; y <- rbinom(nrow(X), 1, burn$p[[cell$truth]])
  }
  if (cell$k > 0) {
    bad <- sample.int(nrow(X), cell$k); col <- if (cell$kind == "burn") "tbsa" else 2
    X[bad, col] <- cell$mult * X[bad, col]
  }
  if (min(sum(y), length(y) - sum(y)) < 3) return(NULL)
  out <- tryCatch(rival_tests(y, X, giviti = cell$giviti), error = function(e) NULL)
  if (is.null(out)) return(NULL)
  c(rep = rep, fp_events = sum(y), fp_x = round(sum(X[, 2]), 6), out)
}

cells <- switch(as.character(study),
  `4` = rbind(CJ(kind = "design", design = c("D1", "D5"), truth = "logit", n = 1000L, k = 0L, mult = 1),
              CJ(kind = "design", design = c("D1", "D5"), truth = "logit", n = 1000L, k = c(1L, 2L, 5L, 10L), mult = c(4, -4))),
  `5` = CJ(kind = "design", design = c("D1", "D2", "D3", "D4", "D5"),
           truth = c("logit", "cloglog", "loglog", "probit", "quad", "inter"), n = c(500L, 1000L, 2000L), k = 0L, mult = 1),
  `6` = rbind(data.table(kind = "map", side = c("cloglog", "cloglog", "loglog", "loglog", "cloglog", "cloglog", "loglog", "cloglog"),
                         lam = c(0, 0, 0, 0, 0.25, 0.5, 0.25, 1), mu = c(2, -2, -2, 2, 1, 0, -1, 0), s = c(2, 2, 2, 2, 1, 1, 1, 2)),
              data.table(kind = "map", side = "cloglog", lam = 1, mu = c(2, -2, 1, -1), s = c(2, 2, 1, 1)))[, `:=`(n = 1000L, k = 0L, mult = 1)],
  `7` = rbind(CJ(kind = "burn", truth = c("logit", "loglog", "cloglog"), k = c(1L, 5L, 10L), mult = 4)))
cells <- unique(cells)
for (col in c("design", "truth", "side")) if (is.null(cells[[col]])) cells[[col]] <- NA_character_
for (col in c("lam", "mu", "s", "n")) if (is.null(cells[[col]])) cells[[col]] <- NA_real_
cells[, giviti := study %in% c(4L, 5L, 7L)]
cells[, id := ifelse(kind == "map", sprintf("map_%s_l%.2f_mu%+d_s%.1f", side, lam, as.integer(mu), s),
              ifelse(kind == "burn", sprintf("burn_%s_k%02d", truth, k),
                     sprintf("%s_%s_n%05d_k%02d_m%s", design, truth, n, k, gsub("-", "neg", mult))))]
cells[, seed := 6e8 + study * 1e7 + .I * 1e4]
B <- if (study == 5) 1000L else 2000L

cl <- makePSOCKcluster(20)
clusterExport(cl, c("one", "ao", "burn", "gen_src", "ROOT"))
invisible(clusterEvalQ(cl, { eval(parse(text = gen_src)); source(file.path(ROOT, "code", "rivals.R")); NULL }))
for (j in seq_len(nrow(cells))) {
  cell <- as.list(cells[j]); f <- file.path(OUT, paste0(cell$id, ".csv.gz"))
  if (file.exists(f)) next
  t0 <- Sys.time()
  R <- rbindlist(lapply(clusterApplyLB(cl, seq_len(B), one, cell = cell),
                        function(v) if (is.null(v)) NULL else as.data.table(as.list(v))))
  dup <- sum(duplicated(R[, .(fp_events, fp_x)]))
  fwrite(R, f)
  cat(sprintf("%-38s B=%d dup=%d | EF %.3f HL %.3f PH %.3f Ts %.3f Xie %.3f Stk %.3f Cub %.3f Giv %.3f [%.0fs]\n",
              cell$id, nrow(R), dup, mean(R$EF < .05, na.rm = TRUE), mean(R$HL < .05, na.rm = TRUE),
              mean(R$PH < .05, na.rm = TRUE), mean(R$Tsiatis < .05, na.rm = TRUE), mean(R$Xie < .05, na.rm = TRUE),
              mean(R$Stukel < .05, na.rm = TRUE), mean(R$Cubic < .05, na.rm = TRUE), mean(R$GiViTI < .05, na.rm = TRUE),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
stopCluster(cl)
cat("study", study, "done\n")
