# Addendum 1 to STUDY_SPEC_EF_ESJ.md (72603ab2) — the plasmode study

Written 2026-10-02 before any replicate of this study ran; hashed with SHA-256.
Code: `code/run_plasmode.R`.

## Data
The 1000 real patients of the burn-injury study (Hosmer, Lemeshow and Sturdivant 2013; `aplore3::burn1000`),
covariates kept exactly as recorded: age, tbsa (percentage of body surface burned), race, inh_inj, flame.
The working model is the textbook linear logit model death ~ age + tbsa + race + inh_inj + flame. Its fitted
linear predictor on the real outcomes, eta_i, defines the risk scores; the real outcomes are not used again.

## Truths (outcomes simulated; intercept a chosen so that the mean risk equals the observed 15%)
- logit (the null): p_i = plogis(a + eta_i)
- loglog-type:      p_i = exp(-exp(-(a + eta_i)))
- cloglog-type:     p_i = 1 - exp(-exp(a + eta_i))

## Predictions written before the run (from Theorem 3 / the sign map: risks are low, about 15%)
- loglog-type: A < 0, EF more powerful than HL (size-adjusted).
- cloglog-type: A > 0, EF less powerful than HL.
A is computed from the pseudo-true logit fit to the true risks on the same 1000 covariate rows.

## Corrupted records (logit truth)
k in {1, 2, 5, 10} patients chosen at random have tbsa multiplied by 4 before the fit (a data-entry error).
Prediction: EF and HL stay at or below 0.10; Stukel's joint score test exceeds 0.10 at k = 1.

## Settings
B = 2000 per cell; G = 10 equal-size groups by rank of the fitted risk, ties broken at random; size-adjusted
power uses the logit-truth cell as the matched null. All per-replicate results written to results/plasmode/.
