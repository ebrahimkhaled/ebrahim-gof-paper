## verify_N1_second_order.R -- the second-order power theorem (N1) against the simulated gains of Study 2.
## Gaussian model r ~ N(m, Sigma) from the pseudo-true fit (as in theorem2_predictions.R); Q = r'r, L = b'r.
##   exact (Gaussian):   P(Q - L > q) - P(Q > q)                  by Monte Carlo with 2e6 draws
##   second order:       -f_Q(q) E[L | Q = q]                     (kernel estimate on the same draws)
##   projection form:    -rho(q) f_Q(q) A,  rho = sqrt(q)/|m| I_{k/2}(|m| sqrt q) / I_{k/2-1}(|m| sqrt q), k = round(sum lam)
## q is the 5% point of the null law of Q. Compared with the simulated size-adjusted gain EF - HL (study2_summary.csv).
ROOT <- "."
source(file.path(ROOT, "code", "ef_exact.R"))
suppressMessages(library(data.table))
src <- readLines(file.path(ROOT, "code", "run_ef_studies.R"))
eval(parse(text = src[grep("^design_data <- function", src):(grep("^## one replicate", src) - 1)]))
SIM <- fread(file.path(ROOT, "results", "study2_summary.csv"))
G <- 10; N <- 1e6; R <- 2e6
out <- list()
for (design in c("D1", "D3", "D4")) for (truth in c("cloglog", "loglog", "probit", "quad", "inter")) {
  set.seed(777)
  d <- design_data(design, N); pt <- true_p(truth, d); X <- cbind(1, d$X); y <- rbinom(N, 1, pt)
  ps <- glm.fit(X, y, family = binomial())$fitted.values
  g  <- ef_groups(ps, G); ng0 <- tabulate(g, G)
  dg <- as.numeric(rowsum(pt - ps, g, reorder = TRUE)) / ng0
  pb <- as.numeric(rowsum(ps, g, reorder = TRUE)) / ng0
  Sig <- ef_parts(y, ps, X, G)$Sigma
  es <- eigen(Sig, symmetric = TRUE); Lh <- es$vectors %*% diag(sqrt(pmax(es$values, 0))) %*% t(es$vectors)
  set.seed(1)
  Z0 <- matrix(rnorm(R * G), R) %*% Lh                        # N(0, Sigma) draws, reused for every n
  qQ <- quantile(rowSums(Z0^2), .95)                          # 5% point of the null law of Q
  for (n in c(200, 500, 1000, 2000)) {
    dsg <- design; trt <- truth; nn <- n
    V <- (n / G) * pb * (1 - pb); m <- (n / G) * dg / sqrt(V); b <- (1 - 2 * pb) / sqrt(V)
    Rm <- sweep(Z0, 2, m, "+"); Q <- rowSums(Rm^2); L <- as.numeric(Rm %*% b)
    exact <- mean(Q - L > qQ) - mean(Q > qQ)
    h <- 0.05 * sd(Q); w <- abs(Q - qQ) < h
    fQ <- mean(w) / (2 * h); EL <- mean(L[w])
    second <- -fQ * EL
    A <- sum(b * m); mm <- sqrt(sum(m^2)); k <- round(sum(es$values))
    rho <- if (mm > 1e-8) sqrt(qQ) / mm * besselI(mm * sqrt(qQ), k / 2, TRUE) / besselI(mm * sqrt(qQ), k / 2 - 1, TRUE) else 0
    proj <- -rho * fQ * A
    s <- SIM[SIM$design == dsg & SIM$truth == trt & SIM$n == nn]
    out[[length(out) + 1]] <- data.table(design, truth, n, A = A, gauss_exact = exact, second_order = second,
                                         projection_form = proj, simulated = s$gain_adj)
  }
  cat(design, truth, "\n")
}
O <- rbindlist(out)
fwrite(O, file.path(ROOT, "results", "verify_N1_second_order.csv"))
options(width = 200); print(O, digits = 3)
cat("\nsign agreement (simulated vs second-order) where |simulated| > 0.01:",
    O[abs(simulated) > .01, sum(sign(simulated) == sign(second_order))], "of", O[abs(simulated) > .01, .N],
    "\ncorrelation simulated ~ second_order:", round(cor(O$simulated, O$second_order), 3),
    "  ~ gauss_exact:", round(cor(O$simulated, O$gauss_exact), 3), "\n")
