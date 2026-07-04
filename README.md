# Reproduction materials — Ebrahim–Farrington goodness-of-fit paper

Simulation code, results, and figure scripts for:

> **Goodness-of-fit testing for logistic regression: when does a directional correction to the
> Hosmer–Lemeshow test help?**
> Ebrahim Khaled Ebrahim and Ahmed El-Kotory (submitted to *Statistics in Medicine*, 2026).

The paper studies the **Ebrahim–Farrington (EF)** test — the Hosmer–Lemeshow statistic modified by a
single directional correction term (the grouped Osius–Rojek/Farrington standardization) — and shows,
through an *alignment functional* `A(δ)`, exactly where the correction changes the test.

## The test and comparators

The EF test and the full partition family of comparators are implemented in the R package
**`ebrahim.gof`**, available on CRAN:

```r
install.packages("ebrahim.gof")
```

<https://CRAN.R-project.org/package=ebrahim.gof>

## Datasets

All real datasets are **public** and loaded from named R packages — none are redistributed here:

| dataset | R package |
|---|---|
| low birth weight (`birthwt`), Pima Indians (`Pima.tr`, `Pima.te`) | `MASS` |
| `kyphosis` | `rpart` |
| GLOW (`glow500`), ICU (`icu`) | `aplore3` |

## Contents

### `R/` — simulation scripts
| script | reproduces |
|---|---|
| `definitive_sim.R` | Type-I size study + Aranda–Ordaz link power sweep (Tables 2, 3, 5; Figs 4, 7, 8), K = 5000, 24-core parallel |
| `within_family2.R` | within-partition-family comparison |
| `tradeoff_fast.R` | EF / HL / Stukel size-adjusted trade-off (Table 7) |
| `or_tradeoff.R` | ungrouped Osius–Rojek trade-off column (via `ebrahim.gof::run.all.gof`) |
| `readout_demo.R` | directional read-out demonstration (Figure 1) |
| `real_data_analysis.R` | real-data concordance tables (Tables 6, 8) |
| `alignment_A_delta.R` | `A(δ)` for the cloglog / log-log links and the quadratic departure |
| `stukel_power.R` | Stukel-test power on the exact cloglog departure |

### `figures/` — Python plotting scripts (matplotlib) that build the figures from the CSVs

### `results/` — simulation-output CSVs (`definitive_results.csv`, `Master_full.csv`, etc.)

## How to reproduce

1. Install dependencies:
   ```r
   install.packages(c("ebrahim.gof", "MASS", "rpart", "aplore3", "ResourceSelection"))
   ```
2. Run the scripts in `R/`; each writes its results to a CSV. **Note:** some scripts set an absolute
   working path at the top — edit it to your local clone location before running.
3. Run the scripts in `figures/` (Python 3 + `numpy` + `matplotlib`) to regenerate the figures from the CSVs.

## Citation

If you use the test or this code, please cite the paper and the package:

> Ebrahim, E. K. (2026). *ebrahim.gof: Ebrahim–Farrington Goodness-of-Fit Test for Logistic Regression.*
> R package version 2.1.0. <https://doi.org/10.32614/CRAN.package.ebrahim.gof>

## License

Code is released under the MIT License (see `LICENSE`). The datasets belong to their respective
R-package maintainers.
