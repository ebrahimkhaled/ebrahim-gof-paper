## smoke_revision.R -- one replicate of every kind of cell in run_revision.R, run serially, before the full studies.
ROOT <- "."
src <- readLines(file.path(ROOT, "code", "run_revision.R"))
src[grep("^study <- ", src)] <- "study <- 4L"
i <- grep("^cl <- makePSOCKcluster", src)
eval(parse(text = src[1:(i - 1)]))
eval(parse(text = gen_src)); source(file.path(ROOT, "code", "rivals.R"))
cat("study 4 cells:", nrow(cells), "\n")
show <- function(cell) { r <- one(1, cell); if (is.null(r)) cat("NULL\n") else print(round(r, 3)) }
show(as.list(cells[1])); show(as.list(cells[nrow(cells)]))
show(list(kind = "burn", truth = "loglog", k = 5L, mult = 4, giviti = TRUE, seed = 1))
show(list(kind = "map", side = "cloglog", lam = 0, mu = 2, s = 2, n = 1000L, k = 0L, mult = 1, giviti = FALSE, seed = 1))
show(list(kind = "design", design = "D2", truth = "cloglog", n = 500L, k = 0L, mult = 1, giviti = TRUE, seed = 7))
