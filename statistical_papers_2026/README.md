# Reproduction materials, 2026 revision (Statistical Papers)

> **A modified Hosmer–Lemeshow goodness-of-fit test for asymmetric links: second-order power and robustness**
> Ebrahim Khaled Ebrahim and Ahmed El-Kotory (submitted to *Statistical Papers*, 2026).

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

## Code → paper

| script | produces |
|---|---|
| `code/ef_exact.R` | the statistic: groups with random tie-breaking, `HL`, `C`, `EF`, exact and chi-squared p-values |
| `code/rivals.R` + `code/lib/` | the eight tests of the screen (validated ports of `ebrahim.gof` 2.8.0) |
| `code/run_ef_studies.R`, `analyse_study2.R` | Studies 1–3 (Table 2, the size rule) |
| `code/run_revision.R`, `analyse_revision.R`, `analyse_study6.R` | Studies 4–7 (Tables 1, 3; Figs 1, 2; claims a–f) |
| `code/run_plasmode.R` | the plasmode study (Table 4) |
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
