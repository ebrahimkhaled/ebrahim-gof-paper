# Reproduction materials, 2026 revision (Statistical Papers)

> **A modified Hosmer–Lemeshow goodness-of-fit test for asymmetric links: second-order power and robustness**
> Ebrahim Khaled Ebrahim, Ibrahim Galal Khattab and Ahmed El-Kotory (submitted to *Statistical Papers*, 2026).

This folder reproduces every table, figure and number of the revised paper and its Online Resource 1.
Run all scripts **from this folder** (`statistical_papers_2026/`); paths are relative to it.
The test itself is in the CRAN package [`ebrahim.gof`](https://CRAN.R-project.org/package=ebrahim.gof) (`ef.gof()`).

## Requirements

R ≥ 4.4 with `data.table`, `CompQuadForm`, `parallel` (base), `aplore3` (burn and other data), `VGAM` (flour beetles),
`MASS`, `doBy`; Python 3 for `code/make_esm_tables.py`. The simulation runners use a PSOCK cluster; set the number of
workers at the top of each runner to suit your machine.

## Declarations (pre-registration)

`declarations/` holds the study specifications with the claims they test. Each was fixed by its SHA-256 hash before
the studies it covers were run:

| file | covers |
|---|---|
| `STUDY_SPEC_EF_ESJ.md` | Studies 1–3: size, power, corrupted records |
| `ADDENDUM1_plasmode.md` | the burn-injury plasmode study and its predictions |
| `ADDENDUM2_rivals_revision.md` | Studies 4–7: the robustness screen of eight tests, power of the survivors, simulated map cells, plasmode with corruption; claims (a)–(f) |
| `ADDENDUM3_inrange_calibrated.md` | Studies 8–10: in-range errors, misclassified responses and Stukel-W; power with the event rate held fixed; computing time; claims (g)–(j), (e′), (f′). Deposited publicly (release v2.1.0) before any of its replicates ran |

| `ADDENDUM4_large_sample.md` | Studies 11–13: the score-test property, the grouping rule of Paul et al. with normal references, a combined large-sample test; claims (k)–(p). Deposited publicly (release v2.3.0) before any of its replicates ran |

| `ADDENDUM5_support.md` | The SUPPORT application (declared analysis, no claim made in advance). Deposited publicly (release v2.5.0) before the analysis |

`results/claims_ledger.csv` (built by `code/analyse_claims.R`) gives every declared claim with its outcome, including
those that failed.

## Code → paper

| script | produces |
|---|---|
| `code/ef_exact.R` | the statistic: groups with random tie-breaking, `HL`, `C`, `EF`, exact and chi-squared p-values |
| `code/rivals.R` + `code/lib/` | the eight tests of the screen (validated ports of `ebrahim.gof` 2.8.0) |
| `code/run_ef_studies.R`, `analyse_study2.R` | Studies 1–3 (Table 2, the size rule) |
| `code/run_revision.R`, `analyse_revision.R`, `analyse_study6.R` | Studies 4–7 (Tables 1, 3; Figs 1, 2; claims a–f) |
| `code/run_plasmode.R` | the plasmode study (Table 4) |
| `code/run_addendum3.R` | Studies 8 and 9 (Tables 1, 3, 4a; Fig 2; Tables S9–S10) |
| `code/time_rivals.R` | Study 10, computing time (Table S8) |
| `code/support_application.R`, `support_explore.R` | the SUPPORT application (Section 6.3, Table S16); needs the public file `support2.csv` (set `SUPPORT_CSV`) |
| `theory/FITTED_GROUPING_PROOF.md` | the theory for groups formed on the fitted risk (Proposition S1 of Online Resource 1) |
| `code/large_sample.R`, `run_addendum4.R`, `check_large_sample.R` | Studies 11–13 (Table 5, Fig 5, Tables S14–S15); the score identity and its orthogonality check |
| `code/analyse_claims.R`, `summarise_study9.R` | verdicts on every declared claim (Table S13) and the Study 9 numbers quoted in the text |
| `code/verify_N1_second_order.R`, `second_order_calibration.R` | second-order theory against simulation (Fig 3) |
| `code/N2_sign_map.R` | the map of predicted gains (Fig 1) |
| `code/theorem2_predictions.R`, `check_identity.R`, `check_pigeon_heyse.R` | numerical checks of Theorem 2, Lemma 2, Proposition 1 (Table S7) |
| `code/influence_one_record.R` | the one-record influence illustration (Proposition 2) |
| `code/beetle_illustration.R`, `realdata_screen.R` | the flour-beetle analysis (Fig 5) and the seven-data-set screen (Table S6) |
| `code/explore_N3_signfree.R`, `explore_refined_null.R`, `sparse_threshold.R` | exploratory checks reported in the text |
| `theory/sparse_sim/` | many-small-groups simulations (Fig 4; Tables S1–S5); `SPARSE_GROUP_THEOREM.md` is the working note of Theorems 4–5 |
| `code/make_figures_tables_v2.R` | all figures and main tables from `results/` |
| `code/make_esm_tables.py` | Tables S1–S6 of Online Resource 1 |

`results/` holds the per-replication output of every study (one file per cell, with data fingerprints) and the
summaries the tables are built from, so the tables can be rebuilt without rerunning the simulations.

## Note on AI assistance

An AI assistant was used to write and run parts of the simulation and analysis code under the authors' direction;
the authors checked the results and are responsible for them.
