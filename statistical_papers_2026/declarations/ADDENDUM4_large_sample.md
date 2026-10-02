# Addendum 4 — the score-test property, the grouping rule of Paul et al., and a combined large-sample test

Written 2026-10-03, before any replicate below ran (a null-only check of the code at n = 5000 and a structural smoke
test with 10-20 replicates per cell were run first; their rates were not used). Hashed with SHA-256 and deposited
publicly (GitHub release of ebrahimkhaled/ebrahim-gof-paper, archived by Zenodo) before the runs.
Code: `code/large_sample.R` (statistics), `code/run_addendum4.R` (runner).

Statistics (`ls_stats`): groups of equal size (differing by at most one) by rank of the fitted risk, ties at random.
- HL-chi2, EF-chi2: chi-square(G - 2) references.
- HL-N, EF-N: one-sided normal references of Theorem 3 (plug-in moments at the fitted coefficients).
- L: the removed component C / sigma_L, two-sided (to first order the Osius-Rojek statistic).
- L-W: the same with the weights zeta = (1 - 2 pi) / {pi (1 - pi)} winsorized at the 2.5% and 97.5% quantiles of
  zeta over the records; two-sided.
- Combined: Fisher's combination of EF-N and L-W, -2 {log p_EF-N + log p_L-W} against chi-square(4).
Two groupings: G of Paul, Pennell and Lemeshow (2013), G = max(10, min(n1/2, (n - n1)/2, 2 + 8 (n/1000)^2)), and
groups of about 25 records (G = floor(n/25)). Also: EF and HL with G = 10 and chi-square(8); Stukel; Stukel-W.
B = 2000 per cell. Size-adjusted power against the matched logit-truth cell.

## Study 11 — the score-test property
Designs D1 and W (eta = 1.2 x1 + 0.5 x2, wide risk range); n = 4000. Alternative: records with similar true risk
share a deviation, logit pi_i = eta_i + u_b(i), with b the bins of 25 records by rank of eta and u_b iid N(0, tau^2),
tau in {0.3, 0.6}; matched null tau = 0.
Claim (p): EF-N has size-adjusted power at least that of HL-N, within 0.01, in every alternative cell, with groups of
25 records and with Paul's G.

## Study 12 — Paul's grouping rule and the normal references
Designs D1, D2, D4; n in {2000, 5000, 10000, 20000}; truths logit (null), cloglog, loglog, omitted square, each
departure with the event rate of its null (intercept shifts of Study 9).
Claims:
 (k) with Paul's G, HL-chi2 or EF-chi2 has size outside 0.05 +- 0.015 in at least one null cell with n >= 10000;
 (l) EF-N, HL-N, L-W and the combined test have size within 0.05 +- 0.015 in every null cell, both groupings;
 (m) with Paul's G, EF-N has higher size-adjusted power than HL-chi2 in every alternative cell with A < 0 (A as
     computed for Study 9 at G = 10);
 (n) with Paul's G, the combined test has higher size-adjusted power than HL-chi2 in at least 80% of the
     alternative cells.

## Study 13 — gross errors with many groups
Design D1, n = 5000, logit truth; k in {0, 5, 25} records with x1 multiplied by 4 or by -4 before the fit.
Claim (o): with Paul's G, EF-N and the combined test stay at or below 0.10 at k = 5 for both errors, and L (not
winsorized) exceeds 0.10 at k = 25 with x1 multiplied by -4.

Everything else is reported as found.
