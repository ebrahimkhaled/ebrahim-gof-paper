# Addendum 3 — in-range errors, response misclassification, a bounded-influence Stukel test, and event-rate-matched power

Written 2026-10-02, after the second mock review (review/REFEREE_SP_GOF_draft2.md, points M2, M4, M5) and before any
replicate below ran. Hashed with SHA-256 and deposited publicly (GitHub release of ebrahimkhaled/ebrahim-gof-paper,
archived by Zenodo) before the runs. Code: `code/rivals.R` (tests), `code/run_addendum3.R` (runner).
Grouped tests use G = 10 equal-size groups by fitted risk with random tie-breaking. Rates are at the nominal 5% level,
raw, and for power also size-adjusted (95% point of the matched null).

New test. "Stukel-W": Stukel's joint score test with the added directions built from the fitted linear predictor
winsorized at its 2.5% and 97.5% sample quantiles (a Mallows-type bound on the leverage of the added directions);
everything else as in Stukel's joint score test.

## Study 8 — robustness screen with in-range errors and response misclassification
Logit truth; designs D1 and D5; n = 1000; B = 2000. Errors applied before the fit to k records chosen at random:
 - covariate x1 multiplied by 4 or by -4, k in {1, 2, 5, 10} (repeats Study 4 with the new test);
 - sign error x1 -> -x1, k in {1, 2, 5, 10, 20} (stays inside the covariate range in both designs);
 - response misclassification y -> 1 - y, k in {5, 10, 20, 50};
 - and k = 0.
Plus the burn plasmode (logit truth, Addendum 1 construction) with tbsa x4 capped at 100 (an in-range value),
k in {0, 1, 5, 10}.
Tests: EF, HL, Pigeon-Heyse, Tsiatis, Xie, Stukel, Stukel-W, cubic calibration LR, GiViTI.
Claims written before the run (each can fail):
 (g) with in-range sign errors in D1, EF, HL and Pigeon-Heyse stay at or below 0.10 up to k = 10, and at least one
     of Stukel, cubic LR and GiViTI exceeds 0.10 by k = 10;
 (h) with response misclassification in D1, EF, HL and Pigeon-Heyse stay at or below 0.10 up to k = 20;
 (i) Stukel-W has a lower false-alarm rate than Stukel at k = 1 with x1 multiplied by -4 in D1;
 (j) in the plasmode with capped tbsa, EF stays at or below 0.10 up to k = 10.
Everything else is reported as found.

## Study 9 — power of the risk-grouped tests with the event rate held fixed
Designs D1-D5; truths logit (null), cloglog, loglog, probit, omitted square, omitted interaction, as in Study 5, but
with the intercept of each departure shifted so that its population event rate equals that of the logit null of the
same design (shift found on 10^6 draws of the covariates, seed 777). n in {500, 1000, 2000}; B = 2000.
Tests: EF, HL, Pigeon-Heyse; Tsiatis, Xie, Stukel, Stukel-W, cubic LR as references (GiViTI omitted for run time).
The alignment A of each calibrated departure is computed as in Study 5 (pseudo-true fit on 10^6 records, ten groups).
Claims:
 (e') among EF, HL and Pigeon-Heyse, EF has the highest size-adjusted power in every cell with A < 0 in which
      |EF - HL| > 0.01;
 (f') the sign of the size-adjusted EF - HL difference is opposite to the sign of A in at least 90% of the cells with
      |EF - HL| > 0.01.
Cells in which every risk-grouped test has size-adjusted power at or below 0.05 are reported separately.

## Study 10 — computing time (descriptive, no claim)
One data set of design D1 at n = 1000 and n = 5000: elapsed time of EF (with HL), Stukel, the le Cessie-van
Houwelingen test (ebrahim.gof, its original O(n^3) form) and BAGofT (the BAGofT package, its original implementation,
default settings; n = 1000 only), median of repeated runs where a run takes under a minute.
