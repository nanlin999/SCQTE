#!/usr/bin/env Rscript

# Run one user-specified simulation setup and write method-level results.
arguments <- commandArgs(trailingOnly = TRUE)
usage <- paste(
  "Usage: Rscript examples/run_scenario.R [options]",
  "  --n N                    Total sample size (default: 10000)",
  "  --sites M                Number of equal-size sites (default: 20)",
  "  --site-sizes SPEC        Explicit sizes, e.g. 60x100,40x100",
  "  --covariates MODEL       uniform or shift (default: uniform)",
  "  --error MODEL            normal or pareto (default: normal)",
  "  --taus LIST              Comma-separated quantile levels (default: 0.1)",
  "  --seed INTEGER           Seed for the first replication (default: 667)",
  "  --replications N         Number of replications (default: 1)",
  "  --reverse-sites BOOL     Reverse complete site blocks (default: false)",
  "  --intervals BOOL         Compute confidence intervals (default: false)",
  "  --output FILE            Output CSV (default: results/example.csv)",
  "  --help                   Show these options",
  sep = "\n"
)
if (identical(arguments, "--help")) {
  cat(usage, "\n")
  quit(save = "no", status = 0L)
}
keys <- arguments[seq_along(arguments) %% 2L == 1L]
values <- arguments[seq_along(arguments) %% 2L == 0L]
if (length(keys) != length(values)) stop("Options need values; see --help.")
allowed <- paste0("--", c("n", "sites", "site-sizes", "covariates", "error",
                          "taus", "seed", "replications", "reverse-sites",
                          "intervals", "output"))
if (any(!keys %in% allowed)) stop("Unknown option; see --help.")
options <- setNames(as.list(values), sub("^--", "", keys))
option <- function(name, default) {
  if (is.null(options[[name]])) default else options[[name]]
}
parse_sizes <- function(specification) {
  unlist(lapply(strsplit(specification, ",", fixed = TRUE)[[1L]], function(term) {
    parts <- as.integer(strsplit(term, "x", fixed = TRUE)[[1L]])
    rep(parts[1L], if (length(parts) == 1L) 1L else parts[2L])
  }), use.names = FALSE)
}
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.")
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."))
for (file in c("load_core.R", "benchmarks.R", "propensity.R",
               "conditional_qte.R", "unconditional_qte.R",
               "simulation_scenarios.R", "population_targets.R")) {
  source(file.path(root, "R", file))
}

error <- match.arg(option("error", "normal"), c("normal", "pareto"))
covariates <- match.arg(option("covariates", "uniform"), c("uniform", "shift"))
size_spec <- option("site-sizes", NULL)
if (is.null(size_spec)) {
  n <- as.integer(option("n", "10000"))
  sites <- as.integer(option("sites", "20"))
  if (n %% sites != 0L) stop("--n must be divisible by --sites.")
  sizes <- rep(n %/% sites, sites)
  size_spec <- paste0(n %/% sites, "x", sites)
} else {
  sizes <- parse_sizes(size_spec)
  n <- sum(sizes)
  sites <- length(sizes)
}
taus <- as.numeric(strsplit(option("taus", "0.1"), ",", fixed = TRUE)[[1L]])
seed <- as.integer(option("seed", "667"))
replications <- as.integer(option("replications", "1"))
reverse_sites <- match.arg(option("reverse-sites", "false"),
                           c("true", "false")) == "true"
compute_intervals <- match.arg(option("intervals", "false"),
                               c("true", "false")) == "true"
output <- option("output", file.path(root, "results", "example.csv"))

scqte_load_core(root)
truth_unconditional <- scqte_population_qte(covariates, error, sizes, taus)
truth_conditional <- 1 + 1.5 * if (error == "normal") qnorm(taus) else
  (1 - taus)^(-1 / 5)

rows <- list()
for (replication in seq_len(replications)) {
  replication_seed <- as.integer(as.double(seed) + replication - 1)
  scenario <- scqte_scenario(sizes, error, covariates,
                             replication_seed, reverse_sites)
  for (j in seq_along(taus)) {
    tau <- taus[j]
    for (estimand in c("conditional", "unconditional")) {
      fit <- if (estimand == "conditional")
        scqte_conditional(scenario$data, scenario$sizes, tau,
                          compute_intervals = compute_intervals) else
        scqte_unconditional(scenario$data, scenario$sizes, tau,
                            compute_intervals = compute_intervals)
      truth <- if (estimand == "conditional") truth_conditional[j] else
        truth_unconditional[j]
      fit$covariates <- covariates
      fit$error <- error
      fit$n <- n
      fit$sites <- sites
      fit$site_sizes <- size_spec
      fit$reverse_sites <- reverse_sites
      fit$seed <- replication_seed
      fit$replication <- replication
      fit$tau <- tau
      fit$estimand <- estimand
      fit$population_qte <- truth
      fit$absolute_error <- abs(fit$estimate - truth)
      fit$covered <- ifelse(is.finite(fit$lower) & is.finite(fit$upper),
                            fit$lower <= truth & truth <= fit$upper, NA)
      rows[[length(rows) + 1L]] <- fit
    }
  }
}
results <- do.call(rbind, rows)
results <- results[c("covariates", "error", "n", "sites", "site_sizes",
                     "reverse_sites", "seed", "replication",
                     "tau", "estimand", "method", "estimate", "population_qte",
                     "absolute_error", "lower", "upper", "covered",
                     "final_iterations", "converged")]
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write.csv(results, output, row.names = FALSE)
cat("Wrote ", nrow(results), " method-level results to ",
    normalizePath(output), "\n", sep = "")
