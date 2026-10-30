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

From this repository's root directory:

```sh
Rscript examples/run_scenario.R --design equal --error normal --sites 20 --taus 0.1,0.5,0.9 --output results/equal_normal.csv
```

The default is one equal-site, normal-error replication at quantile level 0.1.
Use `--first-replication 1 --replications 100` to run 100 replications of the
*selected* scenario. Each output row identifies the design, method, estimand,
replication, quantile level, estimate, population target, absolute error, and
whether the numerical fit converged. Equal-site designs also report interval
limits and coverage for the methods with intervals. Check the `converged` column
before summarizing results.

Other scenario examples:

```sh
Rscript examples/run_scenario.R --design equal --error pareto --sites 50 --taus 0.9 --output results/equal_pareto.csv
Rscript examples/run_scenario.R --design shift --error normal --sites 20 --taus 0.1 --output results/covariate_shift.csv
Rscript examples/run_scenario.R --design unequal --error normal --sites 200 --pattern 30-70 --taus 0.1 --output results/unequal_30-70.csv
```

The unequal-site patterns are `50-50`, `60-40`, `70-30`, `40-60`, and `30-70`.
The covariate-shift and unequal-site examples mirror their respective study
designs and do not calculate intervals. The outcome model is
`Y = 1 + A + X1 + X2 + (1 + X1 + 1.5 A) epsilon` throughout. For the equal-site
design, `X1` and `X2` are independent uniform variables on `[0,1]` and `[0,2]`;
the covariate-shift design uses site-specific beta distributions. The Pareto
error is positive by construction.

The conditional population effect is `1 + 1.5 q_epsilon(tau)`. The marginal
effect is the contrast between the two *marginal potential-outcome quantiles*;
it is computed by numerical integration over the specified covariate
distribution, not by substituting the conditional coefficient. The one-round
unconditional algorithm uses bandwidth `0.05`; the two-round and other smoothed
fits use `((2 + log(n))/n)^0.4`.

## Verify

```sh
Rscript tests/test_public_code.R
```

The tests check estimates, interval limits, population targets, and convergence
reporting for normal, Pareto, covariate-shift, and unequal-size settings.
Generated CSV files under `results/` are ignored by Git.
