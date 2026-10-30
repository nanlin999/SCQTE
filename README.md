# SCQTE reproducibility code

This repository provides the numerical routines used for sequential collaborative
quantile treatment effect (SCQTE) estimation and a configurable example for one
simulation scenario at a time. The C++ source in `src/scqte_core.cpp` contains
the numerical optimization routines; the R files separate estimation,
propensity-score updates, scenario generation, population targets, and benchmarks.

## Requirements

- R with a working C++ compiler
- R packages `Rcpp`, `RcppArmadillo`, and `quantreg`

Install the packages once in R with
`install.packages(c("Rcpp", "RcppArmadillo", "quantreg"))`.
The example compiles the C++ routines automatically. On macOS, the loader applies
a temporary Fortran linker setting for systems with a stale `FLIBS` configuration;
the setting is restored after compilation.

## Run an example

From this repository's root directory, run `--help` to see all options. For
example, to generate 10,000 observations split among 20 equal-size sites:

```sh
Rscript examples/run_scenario.R --n 10000 --sites 20 --covariates uniform --error normal --taus 0.1,0.5,0.9 --seed 667 --replications 1 --intervals true --output results/equal_normal.csv
```

The default is one equal-site, normal-error replication at quantile level 0.1.
Use `--replications 100` to run 100 successive seeds starting at `--seed`.
For unequal sites, replace `--n` and `--sites` with `--site-sizes`; for example,
`60x100,40x100` means 100 sites of size 60 followed by 100 of size 40.
`--reverse-sites true` reverses complete site blocks after generating data.
Intervals are optional because they add computation; request them with
`--intervals true`.

Each output row identifies the covariate model, error model, site sizes, seed,
method, estimand, replication, quantile level, estimate, population target,
absolute error, convergence status, and, when requested, interval limits and
coverage. Check the `converged` column before summarizing results.

Other scenario examples:

```sh
Rscript examples/run_scenario.R --n 10000 --sites 50 --covariates uniform --error pareto --seed 100667 --taus 0.9 --intervals true --output results/equal_pareto.csv
Rscript examples/run_scenario.R --n 10000 --sites 20 --covariates shift --error normal --seed 200667 --taus 0.1 --output results/covariate_shift.csv
Rscript examples/run_scenario.R --site-sizes 70x100,30x100 --reverse-sites true --covariates uniform --error normal --seed 667 --taus 0.1 --output results/unequal_30-70.csv
```

The covariate-shift option generates site-specific beta-distributed covariates.
The outcome model is
`Y = 1 + A + X1 + X2 + (1 + X1 + 1.5 A) epsilon` throughout. For the equal-site
uniform-covariate example, `X1` and `X2` are independent uniform variables on
`[0,1]` and `[0,2]`. The Pareto error is positive by construction.

The conditional population effect is `1 + 1.5 q_epsilon(tau)`. The marginal
effect is the contrast between the two *marginal potential-outcome quantiles*;
it is computed by numerical integration over the specified covariate
distribution, not by substituting the conditional coefficient. The one-round
unconditional algorithm uses bandwidth `0.05`; the two-round and other smoothed
fits use `((2 + log(n))/n)^0.4`.

Generated CSV files under `results/` are ignored by Git.
