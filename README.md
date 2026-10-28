# Effective sample size for predictions based on the Cox survival model

Code, manuscript and supplementary material for *Effective sample size for predictions based on the Cox survival model: insights into individual uncertainty* (Hackmann, Thomassen, le Cessie, Putter, de Wreede, Steyerberg; submitted to Statistics in Medicine).

The effective sample size of a prediction is N* = S(1 − S) / Var(Ŝ): the number of patients with the same covariates in a binomial experiment that would give the same uncertainty as the prediction. For the Cox model, Var(Ŝ) is the delta-method variance from `survfit` (`coxph` with Efron ties).

## Quick start

All paths are relative to the repository root. Open the R project (or set the working directory) there, so that `here::here()` finds the root.

1. **Run the analyses once**: `source("code/run_results.R")`.
   - This runs every slow analysis and saves the output in `results/*.rds`. Expect hours (the simulation studies); the bootstrap takes minutes.
   - A step is skipped if its output file exists. Delete a file in `results/` to rerun that step.
   - If the cluster fails to start ("workers failed to connect"), set `WORKERS <- 1` at the top of `run_results.R` (sequential, slow) or use fewer workers.
2. **Make the figures**: `source("article_figures.R")` writes `images/fig1.png` to `images/fig7.png` (about a minute).
3. **Render the documents**: `quarto render manuscript.qmd` and `quarto render supplement.qmd` (Word and PDF). The supplement only reads `results/`, so it compiles in seconds to minutes.

Only want to look at the code? You do not need step 1: see [Using the code](#using-the-code).

## Requirements

- R 4.5.1 (the version used for the article) and Quarto.
- `code/compile.R` installs and loads the main packages (`survival`, `mstate`, `simsurv`, `tidyverse`, `patchwork`, ...) if they are missing.
- In addition: `future.apply`, `mice`, `Hmisc` and `quarto`. The GUSTO data in Supplement 12 are downloaded via `Hmisc::getHdata`, so that step needs internet access.

## Files

| File | Content |
|---|---|
| `manuscript.qmd` | Main article. Figures at the end, from `images/`. |
| `supplement.qmd` | Supplementary material S1–S12. Reads `results/`; no long computations. |
| `article_figures.R` | All code for the seven figures of the article. |
| `references.bib`, `american-medical-association.csl`, `reference.docx` | Bibliography, citation style and Word template. |
| `code/run_results.R` | Runs all slow analyses and saves `results/*.rds`. |
| `NOTES_supplement_restructure.md` | What changed in the supplement (old to new numbering, removed material). |

### Code (`code/`)

Estimation (sourced by `compile.R`):

| File | Content |
|---|---|
| `compile.R` | Loads packages and sources the functions below; loads the colon data (`colon_df`) and colours. |
| `survfit_n.R` | `survfit_n()`: adds the effective sample size to a `survfit` object (Kaplan-Meier or Cox). |
| `sf_to_df.R` | `sf_to_df()`: turns a `survfit_n` result into a data frame on a given time grid. |
| `calculate_ess.R`, `modified.R` | Effective sample size at time t, including the modification used when only censoring occurs. |
| `variance_cumhaz.R`, `variance_mstate.R` | Delta-method variance of the Cox survival prediction, as in Klein & Moeschberger and via `mstate`. |
| `plot_effective_n.R` | Plot helper for survival and effective sample size curves. |
| `bootstrap.R` | Earlier bootstrap code. Superseded by `bootstrap_ess.R`. |

Analyses (the output goes to `results/`):

| File | Content | Output | Used in |
|---|---|---|---|
| `ess_realworld.R` | Eight datasets of the `survival` package, backward AIC, horizons at fractions of events | `ess_realworld.rds` | Fig 4, S6 |
| `sim_guardrails.R` | Simulation 2: when is N̂* > N, and coverage of the 95% CI. Scenarios, simulation, summaries, plots | `sim_main.rds`, `sim_main_design.rds`, `sim_csens.rds`, `sim_guardrails.rds`, `sim_guardrails_shape.rds`, `sim_guardrails_cens.rds` | Figs 6 and 7, S7 |
| `sim_sensitivity_plots.R` | Tables and figures of the sensitivity analyses in Supplement 7 | – | S7 |
| `check_trim.R` | Compares Figs 6 and 7 with all 200 profiles and with the extreme profiles left out | `results/check_trim_fig*.png` | – |
| `sim_early.R` | Simulation 1: N* early in follow-up (binary and continuous predictor, no censoring) | `sim_early.rds` | Fig 5, S8 |
| `sim_convergence.R` | Convergence of N* towards the Kaplan-Meier N* | `sim_convergence.rds` | S9 |
| `bootstrap_ess.R` | N* from a non-parametric and a parametric bootstrap variance | `ess_bootstrap_colon.rds` | S11 |
| `other_results.R` | GUSTO linear and logistic regression, Cox models stratified by treatment | `other_results.rds` | S12 |
| `sim_files.R` | Building blocks for simulations (`run_one()`, `lambda()`, `extract_sf()`) | – | S9 |

## Using the code

**N\* for your own Cox model.** `survfit` returns the survival `S` and the standard error of the cumulative hazard (`std.err`, the SE of log S), so the standard error of S is `S * std.err`:

```r
library(survival)

fit <- coxph(Surv(time, status) ~ age + nodes + rx, data = colon_df, ties = "efron")
new <- colon_df[complete.cases(colon_df[, c("age", "nodes", "rx")]), ][1:5, ]   # patients to predict for
tt  <- 5 * 365.25                            # time point (days in the colon data)

sf  <- survfit(fit, newdata = new)
idx <- findInterval(tt, sf$time)             # last event time <= tt
S   <- as.numeric(sf$surv[idx, ])
seS <- S * as.numeric(sf$std.err[idx, ])     # SE of S (delta method)

Nstar <- S * (1 - S) / seS^2                 # effective sample size
data.frame(S, Nstar, N = fit$n)              # compare with the development sample size
```

`survfit_n()` in `code/survfit_n.R` does the same for one patient and also offers other variance estimators and the Kaplan-Meier version. `validate_colon_fig3()` in `ess_realworld.R` reproduces the medians of Figure 3, which is a good first check that everything works.

**The simulation (`code/sim_guardrails.R`).** `make_scen(design)` returns the scenarios and `run_all()` runs them (`workers = 1` runs without a cluster). The designs of the article are `"factorial"` (main simulation: event rate 1–99% × no / uniform censoring, Harrell's C = 0.8, N = 500 / event rate), `"csens"` (Harrell's C), `"oat"` (number of predictors and noise predictors), `"shape_only"` (Weibull shape) and `"censoring"` (random censoring). The result has one row per scenario, time point and covariate profile, with the proportion of replications with N̂* > N (`P_gt_N`) and the coverage of the 95% CI (`coverage`). To try it with a small run, without overwriting the real results:

```r
source("code/compile.R")
scen <- make_scen("factorial", e = c(0.1, 0.9))
res  <- run_all(scen, R = 50, workers = 2)    # 50 instead of 1000 replications; noisy, but shows the shape
plot_admin(res)                               # like Figure 6
saveRDS(res, "results/test_sim_main.rds")     # not a name used by run_results.R
```

The other analyses have small test options too: `run_colon_bootstrap(years = 5, B = 50)` and `run_convergence(n_sim = 10)`.

## Figures and supplements

| Article | Source |
|---|---|
| Fig 1–3 | Colon data, computed in `article_figures.R` |
| Fig 4 | `results/ess_realworld.rds` via `plot_fig4()` |
| Fig 5 | `simulate_early()` in `code/sim_early.R` |
| Fig 6, 7 | `results/sim_main.rds` via `plot_admin()` and `plot_coverage_N()` |

| Supplement | Content |
|---|---|
| S1 | 4D PICTURE collaborators |
| S2 | Variance of the Cox model prediction (delta method) |
| S3 | Efficiency of the Breslow estimator compared to Nelson-Aalen |
| S4 | Distribution of predictors in the colon data |
| S5 | Cox model has a narrower 95% CI than Kaplan-Meier |
| S6 | Effective sample size in eight open real-world datasets |
| S7 | Simulation studies, reported according to ADEMP, with sensitivity analyses |
| S8 | Heuristic exploration of N* > N early in follow-up, and of the asymmetry |
| S9 | Convergence of N* towards the Kaplan-Meier N* |
| S10 | Logistic regression has higher variance than Cox regression |
| S11 | Effective sample size from a bootstrapped variance |
| S12 | Other results (linear and logistic regression, stratified Cox models) |

## Reproducibility notes

- Seeds are fixed (1923 unless stated). The simulation studies use `future` with `future.seed = TRUE`; with a different number of workers the random streams are the same, but check the numbers after changing `WORKERS`.
- The main simulation uses 1000 replications per scenario for 14 scenarios (7 event rates from 1% to 99% × 2 censoring settings). The scenarios with an event rate of 1% have N = 50,000, so use fewer workers if R runs out of memory.
- `TRIM_PROFILES` at the top of `sim_guardrails.R` sets whether the summaries and figures of Simulation 2 leave out the extreme covariate profiles (default `FALSE`: all 200 profiles). It does not change the saved simulation results, so nothing needs to be rerun.
- `ess_realworld.R` and `sim_guardrails.R` both define `FRACS`, `PAL`, `run_all` and `theme_art`. They are sourced in a fixed order in `run_results.R` and the qmd files. Do not source both at the top of a script: whichever is sourced last wins.
- The Word output of `supplement.qmd` has a table of contents and a page break before each supplement. Add page numbers in Word after rendering.
- The figures in the submitted article were edited in Adobe Illustrator (legends only).

## Contact

Toby Hackmann, Department of Biomedical Data Sciences, Leiden University Medical Center.

## Funding

Funded by the European Union under Horizon Europe Work Programme 101057332. Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Health and Digital Executive Agency (HaDEA). The UK team are funded under the Innovate UK Horizon Europe Guarantee Programme, UKRI Reference Number: 10041120.
