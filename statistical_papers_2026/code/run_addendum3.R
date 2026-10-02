## run_addendum3.R -- Studies 8 and 9 of declarations/ADDENDUM3_inrange_calibrated.md (hash 8594a3e6).
##   Rscript run_addendum3.R 8|9 [B]
## One file per cell in results/study<k>/, one row per replicate; a finished cell is skipped on restart. Each row
## carries fingerprints of its data (events, sum of the first covariate, the HL statistic) so that a replicate that
## silently repeats another is caught. A second argument B runs a smoke test into results/smoke<k>/.
ROOT <- "."
suppressMessages({ library(parallel); library(data.table) })
args  <- commandArgs(TRUE)
study <- as.integer(args[1])
smoke <- length(args) > 1
B     <- if (smoke) as.integer(args[2]) else 2000L
OUT   <- file.path(ROOT, "results", paste0(if (smoke) "smoke" else "study", study))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
src <- readLines(file.path(ROOT, "code", "run_ef_studies.R"))
gen_src <- src[grep("^design_data <- function", src):(grep("^## one replicate", src) - 1)]
eval(parse(text = gen_src))

## burn-injury covariates, logit truth (Addendum 1 construction)
burn <- local({
  suppressMessages(library(aplore3))
  d0 <- transform(burn1000, y = as.integer(death == "Dead"))
  f0 <- glm(y ~ age + tbsa + race + inh_inj + flame, binomial, data = d0)
  X0 <- model.matrix(f0); colnames(X0) <- make.names(colnames(X0))
  list(X = X0, p = as.numeric(fitted(f0)))
})

## Study 9: intercept shift that keeps each departure's population event rate at the logit null's
shift_file <- file.path(ROOT, "results", "study9_shifts.csv")
if (study == 9 && !file.exists(shift_file)) {
  SH <- rbindlist(lapply(c("D1", "D2", "D3", "D4", "D5"), function(design) {
    set.seed(777); d <- design_data(design, 1e6); target <- mean(plogis(d$eta))
    rbindlist(lapply(c("logit", "cloglog", "loglog", "probit", "quad", "inter"), function(truth) {
      f <- function(c0) { dd <- d; dd$eta <- d$eta + c0; mean(true_p(truth, dd)) - target }
      c0 <- if (truth == "logit") 0 else uniroot(f, c(-6, 6), tol = 1e-10)$root
      data.table(design, truth, shift = c0, event_rate = target)
    }))
  }))
  fwrite(SH, shift_file)
}

one <- function(rep, cell) {
  set.seed(cell$seed + rep)
  if (cell$kind == "burn") {
    X <- burn$X; y <- rbinom(nrow(X), 1, burn$p)
  } else {
    d <- design_data(cell$design, cell$n); d$eta <- d$eta + cell$shift
    y <- rbinom(cell$n, 1, true_p(cell$truth, d))
    X <- cbind(`(Intercept)` = 1, d$X); colnames(X)[-1] <- paste0("x", seq_len(ncol(X) - 1))
  }
  if (cell$k > 0) {
    bad <- sample.int(nrow(X), cell$k)
    if (cell$err == "mult") X[bad, 2] <- cell$mult * X[bad, 2]
    if (cell$err == "yflip") y[bad] <- 1 - y[bad]
    if (cell$err == "cap") X[bad, "tbsa"] <- pmin(4 * X[bad, "tbsa"], 100)
  }
  if (min(sum(y), length(y) - sum(y)) < 3) return(NULL)
  out <- tryCatch(rival_tests(y, X, giviti = cell$giviti), error = function(e) NULL)
  if (is.null(out)) return(NULL)
  c(rep = rep, fp_events = sum(y), fp_x = round(sum(X[, 2]), 6), out)
}

if (study == 8) {
  des <- function(err, mult, k) CJ(kind = "design", design = c("D1", "D5"), err = err, mult = mult, k = k)
  cells <- rbind(des("none", 1, 0L),
                 des("mult", c(4, -4), c(1L, 2L, 5L, 10L)),
                 des("mult", -1, c(1L, 2L, 5L, 10L, 20L)),
                 des("yflip", 1, c(5L, 10L, 20L, 50L)),
                 data.table(kind = "burn", design = NA_character_, err = "cap", mult = 4, k = c(0L, 1L, 5L, 10L)))
  cells[, `:=`(truth = "logit", n = 1000L, shift = 0, giviti = TRUE)]
  cells[, id := ifelse(kind == "burn", sprintf("burn_cap_k%02d", k),
                       sprintf("%s_%s%s_k%02d", design, err, ifelse(err == "mult", gsub("-", "neg", mult), ""), k))]
} else {
  SH <- fread(shift_file)
  cells <- CJ(kind = "design", design = c("D1", "D2", "D3", "D4", "D5"),
              truth = c("logit", "cloglog", "loglog", "probit", "quad", "inter"), n = c(500L, 1000L, 2000L))
  cells <- merge(cells, SH[, .(design, truth, shift)], by = c("design", "truth"))
  cells[, `:=`(err = "none", mult = 1, k = 0L, giviti = FALSE)]
  cells[, id := sprintf("%s_%s_n%05d", design, truth, n)]
}
cells[, seed := 7e8 + study * 1e7 + .I * 1e4]

cl <- makePSOCKcluster(20)
clusterExport(cl, c("one", "burn", "gen_src", "ROOT"))
invisible(clusterEvalQ(cl, { eval(parse(text = gen_src)); source(file.path(ROOT, "code", "rivals.R")); NULL }))
for (j in seq_len(nrow(cells))) {
  cell <- as.list(cells[j]); f <- file.path(OUT, paste0(cell$id, ".csv.gz"))
  if (file.exists(f)) next
  t0 <- Sys.time()
  R <- rbindlist(lapply(clusterApplyLB(cl, seq_len(B), one, cell = cell),
                        function(v) if (is.null(v)) NULL else as.data.table(as.list(v))))
  dup <- sum(duplicated(R[, .(fp_events, fp_x, HL_stat)]))
  fwrite(R, f)
  r <- function(t) mean(R[[t]] < .05, na.rm = TRUE)
  cat(sprintf("%-26s B=%d dup=%d | EF %.3f HL %.3f PH %.3f Ts %.3f Xie %.3f Stk %.3f StkW %.3f Cub %.3f Giv %.3f [%.0fs]\n",
              cell$id, nrow(R), dup, r("EF"), r("HL"), r("PH"), r("Tsiatis"), r("Xie"), r("Stukel"), r("StukelW"),
              r("Cubic"), r("GiViTI"), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
stopCluster(cl)
cat("study", study, "done\n")
