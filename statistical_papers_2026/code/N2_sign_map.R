## N2_sign_map.R -- where does EF beat HL? The alignment A(delta) and the second-order gain (projection form)
## over a two-sided Aranda-Ordaz family of links and over the location and spread of the linear predictor.
##   true link: p = 1 - (1 + lam e^eta)^(-1/lam)            (lam = 1 logit, lam -> 0 cloglog), "cloglog side"
##              p = (1 + lam e^(-eta))^(-1/lam)             (its mirror image; lam -> 0 loglog), "loglog side"
##   eta = mu + s x, x ~ N(0, 1); the fitted model is the linear logit model in x.
## A does not depend on n; the gain is reported at n = 1000 with G = 10.
ROOT <- "."
source(file.path(ROOT, "code", "ef_exact.R"))
suppressMessages(library(data.table))
ao <- function(eta, lam, side) {
  if (side == "cloglog") { if (lam == 0) 1 - exp(-exp(eta)) else 1 - (1 + lam * exp(eta))^(-1 / lam) }
  else                   { if (lam == 0) exp(-exp(-eta))    else (1 + lam * exp(-eta))^(-1 / lam) }
}
G <- 10; N <- 2e5; n <- 1000
grid <- CJ(side = c("cloglog", "loglog"), lam = c(0, 0.25, 0.5, 0.75), mu = c(-2, -1, 0, 1, 2), s = c(0.5, 1, 2))
res <- vector("list", nrow(grid))
for (j in seq_len(nrow(grid))) {
  gj <- grid[j]
  set.seed(4242 + j)
  x <- rnorm(N); pt <- ao(gj$mu + gj$s * x, gj$lam, gj$side); y <- rbinom(N, 1, pt); X <- cbind(1, x)
  f <- suppressWarnings(glm.fit(X, y, family = binomial())); ps <- f$fitted.values
  g <- ef_groups(ps, G); ng0 <- tabulate(g, G)
  dg <- as.numeric(rowsum(pt - ps, g, reorder = TRUE)) / ng0
  pb <- as.numeric(rowsum(ps, g, reorder = TRUE)) / ng0
  A  <- sum((1 - 2 * pb) * dg / (pb * (1 - pb)))
  ## second-order gain at n, projection form with k = G - 2 and the chi-square(G - 2) 5% point
  V <- (n / G) * pb * (1 - pb); m <- (n / G) * dg / sqrt(V); mm <- sqrt(sum(m^2)); k <- G - 2
  q <- qchisq(.95, k); lam_nc <- mm^2
  fQ <- dchisq(q, k, ncp = lam_nc)
  rho <- if (mm > 1e-8) sqrt(q) / mm * besselI(mm * sqrt(q), k / 2, TRUE) / besselI(mm * sqrt(q), k / 2 - 1, TRUE) else 0
  res[[j]] <- data.table(gj, events = mean(pt), A = A, gain_n1000 = -rho * fQ * A,
                         power_HL = pchisq(q, k, ncp = lam_nc, lower.tail = FALSE))
}
M <- rbindlist(res)
fwrite(M, file.path(ROOT, "results", "N2_sign_map.csv"))
M[, A := round(A, 2)]
options(width = 200)
print(dcast(M, side + lam ~ mu + s, value.var = "A"), digits = 2)
cat("\nshare of the grid where EF gains (A < 0):",
    M[, .(gain_share = round(mean(A < 0), 2)), by = .(side)][, paste(side, gain_share)], "\n")
print(M[, .(gain_share = round(mean(A < 0), 2), mean_gain = round(mean(gain_n1000), 4)), by = .(side, lam)])
print(M[, .(gain_share = round(mean(A < 0), 2)), by = .(side, mu)])
