## sparse_threshold.R -- the usage rule "every risk group carries enough expected events": for each Study-1 cell,
## the median over replicates of min_g min(e_g, n_g - e_g) (the smallest expected count of events or non-events
## in any of the ten groups), beside EF's realised size, to find where the conservativeness starts.
ROOT <- "."
source(file.path(ROOT, "code", "ef_exact.R"))
suppressMessages(library(data.table))
src <- readLines(file.path(ROOT, "code", "run_ef_studies.R"))
eval(parse(text = src[grep("^design_data <- function", src):(grep("^## one replicate", src) - 1)]))
out <- list()
for (design in c("D1", "D2", "D3", "D4", "D5")) for (n in c(100, 200, 500, 1000, 5000)) {
  mins <- replicate(200, {
    d <- design_data(design, n); y <- rbinom(n, 1, plogis(d$eta)); X <- cbind(1, d$X)
    p <- suppressWarnings(glm.fit(X, y, family = binomial()))$fitted.values
    z <- ef_parts(y, p, X, 10); ng <- tabulate(ef_groups(p, 10), 10)
    min(pmin(z$e, ng - z$e))
  })
  S <- fread(file.path(ROOT, "results", "study1", sprintf("%s_logit_n%05d_k00_m1.csv.gz", design, n)))
  out[[length(out) + 1]] <- data.table(design, n, min_expected = median(mins),
                                       size_EF = mean(S$p_EF_chisq < .05), size_HL = mean(S$p_HL_chisq < .05))
}
O <- rbindlist(out)[order(min_expected)]
fwrite(O, file.path(ROOT, "results", "sparse_threshold.csv"))
print(O, digits = 3)
