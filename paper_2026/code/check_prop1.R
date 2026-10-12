## Proof-check numerics for the influence bound (Proposition 3 of the TEST version): |Delta EF| at fixed coefficients
## on the scale G (1+R)^2 Vmin^(-1/2), R = max_g |r_g|. The paper's constant is 11, so the printed maximum 0.615 is
## about 6% of the stated bound. (Earlier drafts used the constant 16.)
source("./code/ef_exact.R")

stats_at <- function(y, p, G) {
  g <- ef_groups(p, G); ng <- tabulate(g, G)
  o <- as.numeric(rowsum(y, g, reorder = TRUE)); e <- as.numeric(rowsum(p, g, reorder = TRUE))
  pb <- e / ng; V <- ng * pb * (1 - pb); r <- (o - e) / sqrt(V); b <- (1 - 2 * pb) / sqrt(V)
  list(HL = sum(r^2), C = sum(b * r), EF = sum(r^2) - sum(b * r), r = r, V = V)
}

## (1) realistic: D1, n = 1000, MLE coefficients, replacement by extreme covariates
set.seed(11)
n <- 1000; G <- 10
x1 <- runif(n, -3, 3); x2 <- rbinom(n, 1, .5); y <- rbinom(n, 1, plogis(.6 * x1 + .5 * x2))
bh <- coef(glm(y ~ x1 + x2, family = binomial)); X <- cbind(1, x1, x2)
p0 <- plogis(X %*% bh)[, 1]; s0 <- stats_at(y, p0, G)
R <- max(abs(s0$r)); bound <- 16 * G * (1 + R)^2 / sqrt(min(s0$V))
worst <- 0
for (i in sample(n, 200)) for (mult in c(-1000, -50, -4, 4, 50, 1000)) for (yt in 0:1) {
  xt <- c(1, x1[i] * mult, x2[i]); pp <- p0; pp[i] <- plogis(sum(xt * bh)); yy <- y; yy[i] <- yt
  s1 <- stats_at(yy, pp, G); worst <- max(worst, abs(s1$EF - s0$EF))
}
cat(sprintf("D1 n=1000: R=%.2f Vmin=%.1f  max|dEF|=%.3f  bound=%.1f  ratio=%.4f\n", R, min(s0$V), worst, bound, worst / bound))

## (2) adversarial: small groups (V near 2), arbitrary p and y chosen to make R large, every insertion rank
set.seed(12)
best <- 0; bestHLg <- 0; viol6 <- 0
for (trial in 1:3000) {
  G <- sample(3:6, 1); ng <- sample(9:20, 1); n <- G * ng
  p <- sort(runif(n, .15, .85))
  y <- rbinom(n, 1, p)
  if (runif(1) < .5) { y[1:(n %/% 3)] <- rbinom(n %/% 3, 1, .95) }   # misfit to inflate R
  s0 <- stats_at(y, p, G)
  if (min(s0$V) < 2) next
  R <- max(abs(s0$r)); tmin <- 1 / sqrt(min(s0$V))
  for (i in sample(n, min(n, 15))) for (pt in c(1e-9, runif(3), 1 - 1e-9)) for (yt in 0:1) {
    pp <- p; pp[i] <- pt; yy <- y; yy[i] <- yt
    s1 <- stats_at(yy, pp, G)
    best <- max(best, abs(s1$EF - s0$EF) / (G * (1 + R)^2 * tmin))
    dHLg <- abs(s1$r^2 - s0$r^2) / ((1 + R)^2 / sqrt(s0$V))  # per-group, per own V
    bestHLg <- max(bestHLg, max(dHLg))
  }
}
cat(sprintf("adversarial: max |dEF| / (G (1+R)^2 Vmin^-1/2) = %.3f  (paper constant 16, corrected 11)\n", best))
cat(sprintf("adversarial: max per-group |d r_g^2| / ((1+R)^2 V_g^-1/2) = %.3f (paper 6*sqrt2=%.2f, corrected 4*sqrt2=%.2f)\n",
            bestHLg, 6 * sqrt(2), 4 * sqrt(2)))

## (3) the invalid algebraic step: d(2R+d) <= 6 sqrt2 (1+R)^2 Vmin^-1/2 with d = 2 sqrt2 (1+R) Vmin^-1/2
for (R in c(0.5, 1, 3, 10)) {
  Vm <- 2; d <- 2 * sqrt(2) * (1 + R) / sqrt(Vm)
  cat(sprintf("R=%4.1f Vmin=2: d(2R+d)=%.2f vs 6sqrt2(1+R)^2/sqrt(Vmin)=%.2f\n", R, d * (2 * R + d), 6 * sqrt(2) * (1 + R)^2 / sqrt(Vm)))
}

## (4) Stukel: movement orthogonal to beta-hat keeps S bounded
set.seed(13)
n <- 1000; x1 <- runif(n, -3, 3); x2 <- rnorm(n); y <- rbinom(n, 1, plogis(.6 * x1))
X <- cbind(1, x1, x2); bh <- coef(glm(y ~ x1 + x2, family = binomial))
st <- function(Xm, yv) { e <- (Xm %*% bh)[, 1]; stukel_joint(yv, plogis(e), Xm, e)[1] }
S0 <- st(X, y)
for (L in c(10, 1e3, 1e6)) {
  xt <- c(1, 0, 0); xt[3] <- L; v <- xt; v[2] <- -(bh[1] + bh[3] * L) / bh[2] * 0 # keep eta bounded:
  ## choose x1 so that eta stays 0: x1 = -(b0 + b2 L)/b1
  xt[2] <- -(bh[1] + bh[3] * L) / bh[2] + 1
  Xm <- X; Xm[1, ] <- xt; yy <- y; yy[1] <- 0
  cat(sprintf("||x||=%.1e eta=%.2f  dStukel=%.3f\n", sqrt(sum(xt^2)), sum(xt * bh), st(Xm, yy) - S0))
}
