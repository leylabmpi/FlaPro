# r_helper_lib_flapro

Helper functions for the differential-abundance and compositional (Nearest Balance) analyses
in this repository. This is a trimmed copy of the shared lab library
[leylabmpi/r_helper_lib](https://github.com/leylabmpi/r_helper_lib), pinned at commit
`858d5b3` and reduced to the 14 functions this project uses.

Load it by sourcing both files:

```r
source(file.path(helper_lib_dir, "nb_helpers.R"))
source(file.path(helper_lib_dir, "functions_mb.R"))
```

Function bodies are verbatim from the pinned upstream commit, so they stay diffable against it.
If you need a function that isn't here, take it from upstream and re-pin rather than editing
these files.

## Data conventions

Most functions take **long-format** feature tables with three columns: `Sample`, `feature`,
`value`. The column names for the latter two are configurable via `arg_feature_col` /
`arg_value_col`, but `Sample` is assumed.

Nearest Balance functions instead work on a **CLR matrix** (`arg_clr`): samples in rows,
features in columns, with row names as sample IDs. A "balance" object follows the layout
returned by `NearestBalance::find_nearest_balance_clr`:

```r
list(
  b1     = list(num = <character>,   # features in the numerator group
                den = <character>),  # features in the denominator group
  impact = <numeric>,                # share of variance explained
  sbp    = <data.frame>,             # one column `b1`, one row per feature,
                                     #   +1 = numerator, -1 = denominator, 0 = excluded
                                     #   row names are feature IDs
  coord  = <numeric>                 # balance coordinate
)
```

---

## `functions_mb.R` — linear models and adjustment

### `check_full(x, arg_feature_col)`

Asserts that a long table is *complete* — that it holds one row per sample-feature pair,
including explicit zeros — by checking `n_samples * n_features == nrow(x)`. Raises a
`testthat` failure otherwise. Call this before computing prevalence, since a table with
implicit (dropped) zeros silently inflates prevalence.

### `adjust_via_residuals(arg_meta_samples, arg_features_long, arg_factor_to_adj_for, arg_feature_col = "feature", arg_value_col = "value")`

Removes the effect of one or more covariates from every feature by replacing values with
residuals. Pivots the long table to wide, fits `lm(as.matrix(features) ~ <factors>)` across
all features at once, and returns the residuals back in long format with any extra columns
from the input preserved.

`arg_factor_to_adj_for` is the right-hand side of a formula as a string, e.g. `"Age + Sex"`.
`arg_meta_samples` must have a `Sample` column and cover every sample in the feature table.
Calls `check_full` first.

### `do_lm_tidy(arg_lm_in, arg_formula, arg_response_col = "value", arg_feature_col = "feature")`

Fits an independent `lm` per feature and returns one tidy table of coefficients. Nests by
feature, fits `<arg_response_col> ~ <arg_formula>`, tidies via `broom::tidy`, drops the
intercept rows, and sorts by term then p-value.

`arg_formula` is the right-hand side as a string, e.g. `"DiseaseScore + Age"`. No multiple-testing
correction is applied — adjust p-values yourself downstream.

### `do_lmer_tidy(arg_lm_in, arg_formula, arg_response_col = "value", arg_feature_col = "feature")`

As `do_lm_tidy`, but fits a mixed model with `lme4::lmer` for repeated measures, and drops
random-effect rows (`effect != "ran_pars"`) so only fixed effects are returned. Include the
random term in the formula string, e.g. `"DiseaseScore + Age + (1|Participant_ID)"`.

Requires `lmerTest` and `broom.mixed` to be attached for p-values and `tidy` support.

---

## `nb_helpers.R` — Nearest Balance

### `do_clr_lm_get_fac_coeffs(arg_lm_nb_formula, arg_clr, arg_meta_df, arg_sel_factor_coef)`

Fits a multivariate `lm` of the CLR matrix against the model formula and extracts the
coefficient vector for one factor of interest. Returns a data frame with columns `name`
(feature) and `lm_coef`.

Coefficient row names are sanitised to alphanumerics-plus-underscore, so
`arg_sel_factor_coef` must be given in that sanitised form. For a continuous predictor the
coefficient name is just the variable name; for a categorical one it is the variable
concatenated with the non-reference level (e.g. `GroupCD`).

### `make_nb_formula_to_1factor(arg_train, arg_test, arg_meta_train, arg_formula, arg_factor_coef)`

Runs one Nearest Balance fit: multivariate `lm` on the training CLR matrix, extract the
coefficient row for `arg_factor_coef`, then `find_nearest_balance_clr` on it. Returns a
balance object.

`arg_test` is accepted but unused; pass `NULL` if you have nothing for it. Mainly called
internally by `parallel_run_nb`.

### `parallel_run_nb(arg_lm_nb_formula, arg_sel_factor_coef, arg_clr, arg_meta_df, arg_splits, arg_n_sim, arg_num_rparallel_cores)`

Runs `arg_n_sim` cross-validation iterations of `make_nb_formula_to_1factor` in parallel over
a PSOCK cluster of `arg_num_rparallel_cores` workers, and returns the list of balance objects.

`arg_splits` is a list of at least `arg_n_sim` integer vectors of training-row indices into
`arg_clr` — generate it yourself (e.g. with `caret::createDataPartition`) so the resampling
scheme stays under your control.

Iterations that error are caught and **dropped**, with a message reporting how many. The
returned list can therefore be shorter than `arg_n_sim`; check `length()` before assuming a
one-to-one mapping with your splits.

The worker body is deterministic, so results do not depend on `arg_num_rparallel_cores`. The
cluster is not seeded — if you ever make the worker body stochastic, add explicit parallel RNG
handling.

### `aggregate_balance_iterations(arg_nb_res, arg_reproducibility_threshold)`

Combines balances from many iterations into a consensus. Column-binds the per-iteration `sbp`
vectors, computes for each feature how often it appeared in the numerator (`freq_num`) and the
denominator (`freq_den`), and assigns a consensus sign: `+1` if `freq_num` reaches the
threshold, `-1` if `freq_den` does, `0` otherwise.

Returns a list of two elements:

- `sbp_iters` — per-feature table with all iteration columns plus `freq_num`, `freq_den`,
  `b1_consensus`, and `reprod` (the higher of the two frequencies, i.e. how reproducible that
  feature's assignment was).
- `sbp_consensus` — a minimal `sbp`-shaped data frame of just the consensus signs, ready for
  `format_consensus_balance`.

`arg_reproducibility_threshold` is a proportion in `[0, 1]`, e.g. `0.8`.

> **Requires a variable named `n_sim` in the calling environment.** The frequency denominator
> is taken from `n_sim` rather than from a parameter or the actual number of iteration columns.
> Because failed iterations are dropped upstream, set `n_sim` to the number of iterations
> actually present (`length(arg_nb_res)`) if any failed, or frequencies will read low.

### `get_sbp_consensus_with_reprod(arg_sbp_iters)`

Reduces the `sbp_iters` table to the three columns usually needed downstream: `taxName`,
`b1` (renamed from `b1_consensus`), and `reprod`.

### `format_consensus_balance(arg_sbp_consensus)`

Converts a consensus `sbp` data frame into a balance object with the same shape as
`find_nearest_balance_clr` output, so consensus balances work with the functions below.
Populates `b1$num`, `b1$den`, and `sbp`; `impact` and `coord` are not set.

### `compute_balance(arg_abundance, arg_nb)`

Computes the balance value per sample. Returns a tibble with `Sample` and `NB_Value`.

`arg_abundance` is an abundance matrix, not CLR values — samples in rows, features in columns.
Pass zero-replaced abundances: the underlying `balance::balance.fromSBP` substitutes any
remaining zeros with the smallest non-zero value in the matrix and emits a message, which is
rarely what you want for a log-ratio.

Its column names must match the row names of `arg_nb$sbp` exactly, or `balance.fromSBP` stops
with a component-mismatch error. Subset the matrix to the balance's features first if needed:

```r
keep <- colnames(abundance) %in% c(nb$b1$num, nb$b1$den)
compute_balance(abundance[, keep, drop = FALSE], nb)
```

### `balance_size(arg_nb)`

Returns `list(n_num = <int>, n_den = <int>)`, the number of features on each side of the
balance.

### `join_coefs_and_sbp_reprod(arg_sbp_consensus_reprod, arg_coef_oa_bacs)`

Joins the consensus table (by `taxName`) to a coefficient table (by `name`, as returned by
`do_clr_lm_get_fac_coeffs`) and keeps only features actually in the balance, i.e. those with
`b1` of `1` or `-1`. Useful for relating balance membership to per-feature effect sizes.

### `inverse_clr(arg_clr)`

Inverts a CLR transformation: exponentiates and normalises each row to sum to 1. Returns
closed compositions (relative abundances), not the original counts.

---

## R environment

Developed and tested with **R 4.3.3**.

| Package | Version | Used for |
|---|---|---|
| `dplyr` | 1.1.4 | data manipulation throughout |
| `tibble` | 3.2.1 | row-name/column conversions |
| `tidyr` | 1.3.1 | `pivot_*`, `nest`, `unnest` |
| `purrr` | 1.0.4 | `map`, `list_cbind` |
| `rlang` | 1.1.6 | `as_name` for coefficient selection |
| `broom` | 1.0.8 | `tidy` for `lm` fits |
| `broom.mixed` | 0.2.9.6 | `tidy` for `lmer` fits |
| `lme4` | 1.1.37 | mixed models |
| `lmerTest` | 3.1.3 | p-values for mixed models |
| `testthat` | 3.2.3 | assertions in `check_full` |
| `foreach` | 1.5.2 | parallel iteration |
| `doParallel` | 1.0.17 | parallel backend |
| `NearestBalance` | 0.2.0 | `find_nearest_balance_clr` |
| `balance` | 0.2.5 | `balance.fromSBP` (attached via `NearestBalance`) |

`parallel` is part of base R. Both files call `require()` for their dependencies on source, so
a missing package surfaces as a warning at load time rather than an error mid-analysis.

`NearestBalance` is not on CRAN. Install it from its source repository before using
`nb_helpers.R`:

<https://bitbucket.org/knomics/nearestbalance/src/master/>

```r
remotes::install_bitbucket("knomics/nearestbalance")
```

Installing it also attaches `balance`, which provides `balance.fromSBP`. Refer to the
repository above for current installation instructions if the command changes.
