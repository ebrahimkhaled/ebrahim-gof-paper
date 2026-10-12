# Addendum 5 — a large real data set: the SUPPORT study

Written 2026-10-03, before any goodness-of-fit statistic was computed on these data by this project. Hashed with
SHA-256 and deposited publicly before the analysis. Code: `code/support_application.R`.

## Data and models (fixed in advance, taken unchanged from earlier published work of the authors)
SUPPORT (Study to Understand Prognoses and Preferences for Outcomes and Risks of Treatments), file support2.csv;
outcome in-hospital death; complete cases of age, meanbp, hrt, resp, temp, crea, sod, wblc, num.co, sex; records
permuted once with seed 20260931 (the data preparation of the DeepGOF-1 SUPPORT application).
- M0: logistic model linear in the ten covariates.
- M1: M0 with meanbp replaced by a natural cubic spline with 3 degrees of freedom (the repair reported there).

## Analysis (every number below is reported, whatever it shows)
For M0 and M1:
 - the number of groups G of Paul et al. (2013) and the records per group;
 - the design constant R_n at the fitted risks, and the variance factors 1 - 1/m (EF) and 1 - 1/m + R/(2m) (HL)
   that Theorem 3 gives for the chi-square(G - 2) reference, with the implied size of that reference at 5%;
 - p-values of HL and EF with Paul's G, each with the chi-square(G - 2) reference and with the normal reference of
   Theorem 3; HL and EF with G = 10 and chi-square(8); Stukel's test; Stukel-W;
 - a parametric bootstrap of the null distribution of HL and EF with Paul's G (B = 999 data sets simulated from the
   fitted model M0, model refitted each time), giving bootstrap p-values and the bootstrap size of the
   chi-square(G - 2) and normal references at 5%.
No claim is made in advance about whether the references agree on these data.
