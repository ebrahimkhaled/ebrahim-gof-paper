# Addendum 2 — the revision studies (rivals, partition family, map check, plasmode with corruption)

Written 2026-10-02, after the mock review (review/REFEREE_SP_GOF.md, PROOF_CHECK.md, NUMBERS_AUDIT.md) and before
any replicate below ran. Hashed with SHA-256. Code: `code/rivals.R` (tests), `code/run_revision.R` (runner).
Tests come from the EDGE battery library (simulations/_battery_tests.R, ports of ebrahim.gof 2.8.0 checked to 1e-8).
Grouped tests use G = 10 equal-size groups by fitted risk with random tie-breaking; EF and HL are also computed at
G = 5 and G = 20. All rates are reported raw (nominal 5%) AND size-adjusted (95% point of the matched null).

## Study 4 — how fragile are the directed record-level tests?
Logit truth; designs D1 (two covariates) and D5 (six covariates) of STUDY_SPEC_EF_ESJ.md; n = 1000;
k in {0, 1, 2, 5, 10} records with x1 multiplied by 4 or by -4 before the fit; B = 2000.
Tests: EF, HL, Pigeon-Heyse, Tsiatis, Xie (partition family); Stukel joint score, cubic calibration LR, GiViTI
(directed, record-level).
Claims written before the run (each can fail):
 (a) Stukel and the cubic calibration test exceed 0.10 at k = 1 (x4) in D1;
 (b) GiViTI exceeds 0.10 by k = 2 (x4) in D1;
 (c) EF, HL and Pigeon-Heyse stay at or below 0.10 up to k = 10 (x4) in both designs;
 (d) Tsiatis and Xie, which group in covariate space, exceed 0.10 by k = 10 (x4) in D1.

## Study 5 — power within the partition family, with the directed tests as a ceiling
Designs D1-D5; truths logit (null), cloglog, loglog, probit, omitted square, omitted interaction (as in Study 2);
n in {500, 1000, 2000}; B = 1000. Same tests as Study 4.
Claims: (e) among the risk-ordered partition tests (EF, HL, Pigeon-Heyse), EF has the highest size-adjusted power in the
cells where the alignment A < 0 and HL or Pigeon-Heyse where A > 0; (f) the sign of mean(C) equals the sign of A
(directional read-out) in at least 90% of cells with |A| > 0.2.

## Study 6 — the map's cells simulated
n = 1000, B = 2000, EF and HL only, raw and size-adjusted (matched null: the logit link with the same mu and s).
Cells (side, lambda, mu, s): (cloglog, 0, 2, 2), (cloglog, 0, -2, 2), (loglog, 0, -2, 2), (loglog, 0, 2, 2),
(cloglog, 0.25, 1, 1), (cloglog, 0.5, 0, 1), (loglog, 0.25, -1, 1), and lambda = 1 (logit) at (0, 2) as a null check.
Prediction: the sign of the simulated size-adjusted gain equals the sign of the map's prediction in every cell where
|prediction| > 0.02; the magnitude is reported, not claimed.

## Study 7 — plasmode with corrupted records under the asymmetric truths
The burn-injury plasmode of Addendum 1 (loglog-type and cloglog-type truths); k in {1, 5, 10} patients with tbsa x4;
B = 2000; tests EF, HL, Pigeon-Heyse, Stukel. Report raw and size-adjusted power (null: the logit-truth cell with the
same k). No directional claim beyond Addendum 1.
