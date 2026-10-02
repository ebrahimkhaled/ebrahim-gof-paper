# Study specification — EF for the Egyptian Statistical Journal

Written 2026-10-02, before any replicate of the studies below ran. Hashed with SHA-256.
Code: `code/ef_exact.R` (statistics), `code/run_ef_studies.R` (runner).

## Tests (all at ten groups, equal-size groups by rank of the fitted risk)
- HL-chi2: Hosmer–Lemeshow, chi-square(G − 2).            - HL-exact: HL, its limiting weighted chi-square.
- EF-chi2: EF = HL − C, chi-square(G − 2) (the earlier EF). - EF-exact: EF, Theorem 1 (noncentral weighted chi-square minus a constant, Imhof).
- Stukel: Stukel's joint score test, chi-square(2).

## Designs (covariates; linear predictor eta; the fitted model is always the linear logit model in the covariates)
- D1 base: x1 ~ U(−3, 3), x2 ~ Bernoulli(0.5); eta = 0.6 x1 + 0.5 x2.
- D2 rare: as D1, eta = −2.8 + 0.6 x1 + 0.5 x2 (events about 10%).
- D3 skewed: x1 ~ Exp(1) − 1, x2 ~ Bernoulli(0.5); eta = −0.3 + 0.9 x1 + 0.5 x2.
- D4 strong: as D1, eta = 1.6 x1 + 0.8 x2 (high discrimination).
- D5 six: x1..x4 ~ N(0, 1), x5, x6 ~ Bernoulli(0.3); eta = 0.4 (x1 + … + x6).

## Study 1 — null calibration
Logit truth, designs D1–D5, n ∈ {100, 200, 500, 1000, 5000}, B = 2000 a cell.
Read: realised size at 5% and 1% of each test, with Monte Carlo SE; the KS distance of the p-values.
Claim to be tested (could fail): EF-exact holds 5% within 3 MCSE (±0.015) in every cell with n ≥ 200;
EF-chi2 and HL-chi2 depart from 5% by more than 3 MCSE in at least one cell (the reason the exact law matters).

## Study 2 — power and the alignment theorem
Truths: cloglog (p = 1 − exp(−exp(eta))), loglog (p = exp(−exp(−eta))), probit (p = Phi(eta)),
omitted quadratic (logit, eta + 0.15 x1²), omitted interaction (logit, eta + 0.5 x1 x2). Designs D1, D3, D4;
n ∈ {200, 500, 1000, 2000}; B = 2000. Power raw and size-adjusted (critical values from the matched Study-1 cell; n = 2000 nulls run here with the same design).
Theorem 2 check: in each cell, the simulated mean of C against the predicted A(delta) computed from a pseudo-true fit
on 10^6 records; sign agreement required in every cell where |A| exceeds 2 MCSE of mean(C).
Expected (could fail): EF ≥ HL (size-adjusted) on cloglog and loglog, ≈ HL on probit, ≤ HL on the omitted quadratic.

## Study 3 — corrupted records
Logit truth, design D1, n = 1000; k ∈ {1, 2, 5, 10} records with x1 multiplied by 4 (exaggerated) or by −4 (sign error)
before the fit; B = 2000. Read: false-alarm rates. Claim (could fail): EF-exact stays ≤ 0.10 up to k = 10 exaggerated;
Stukel exceeds 0.10 at k = 1.
Plus a deterministic influence illustration: one record's x1 multiplied by 2…1000, change in EF, HL and Stukel,
with and without refitting, against the bound of Proposition 3.

All per-replicate p-values are written to disk per cell (results/), so every number can be recomputed.
