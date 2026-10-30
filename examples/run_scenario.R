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
allowed <- c("--n", "--sites", "--site-sizes", "--covariates", "--error",
             "--taus", "--seed", "--replications", "--reverse-sites",
             "--intervals", "--output")
names <- arguments[seq_along(arguments) %% 2L == 1L]
values <- arguments[seq_along(arguments) %% 2L == 0L]
if (length(arguments) %% 2L || any(!names %in% allowed) ||
    anyDuplicated(names)) {
  stop("Invalid options. Run with --help for usage.")
}
option <- function(name, default = NULL) {
  hit <- match(name, names)
  if (is.na(hit)) default else values[hit]
}
integer_option <- function(name, default, minimum = 1L) {
  value <- option(name, as.character(default))
  number <- suppressWarnings(as.integer(value))
  if (length(number) != 1L || is.na(number) ||
      number < minimum || as.character(number) != value) {
    stop(name, " must be an integer of at least ", minimum, ".")
  }
  number
}
boolean_option <- function(name, default = "false") {
  value <- option(name, default)
  if (!value %in% c("true", "false")) stop(name, " must be true or false.")
  value == "true"
}
parse_sizes <- function(specification) {
  terms <- strsplit(specification, ",", fixed = TRUE)[[1L]]
  if (!length(terms) || any(!grepl("^[1-9][0-9]*(x[1-9][0-9]*)?$", terms))) {
    stop("--site-sizes must look like 60x100,40x100.")
  }
  sizes <- unlist(lapply(terms, function(term) {
    parts <- strsplit(term, "x", fixed = TRUE)[[1L]]
    count <- if (length(parts) == 1L) 1L else suppressWarnings(as.integer(parts[2L]))
    size <- suppressWarnings(as.integer(parts[1L]))
    if (is.na(size) || is.na(count)) stop("Site sizes and counts must fit in R integers.")
    rep(size, count)
  }), use.names = FALSE)
  if (length(sizes) < 2L || anyNA(sizes) || !is.finite(sum(sizes))) {
    stop("--site-sizes must specify at least two valid sites.")
  }
  as.integer(sizes)
}
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.")
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."))
for (file in c("load_core.R", "benchmarks.R", "propensity.R",
               "conditional_qte.R", "unconditional_qte.R",
               "simulation_scenarios.R", "population_targets.R")) {
  source(file.path(root, "R", file))
}

error <- option("--error", "normal")
covariates <- option("--covariates", "uniform")
if (!error %in% c("normal", "pareto") ||
    !covariates %in% c("uniform", "shift")) {
  stop("Choose --error normal|pareto and --covariates uniform|shift.")
}
size_spec <- option("--site-sizes")
if (is.null(size_spec)) {
  n <- integer_option("--n", 10000L)
  sites <- integer_option("--sites", 20L, 2L)
  if (n %% sites != 0L) stop("--n must be divisible by --sites.")
  sizes <- rep(n %/% sites, sites)
  size_spec <- paste0(n %/% sites, "x", sites)
} else {
  sizes <- parse_sizes(size_spec)
  n <- sum(sizes)
  sites <- length(sizes)
  if (!is.null(option("--n")) && integer_option("--n", n) != n) {
    stop("--n must equal the sum of --site-sizes.")
  }
  if (!is.null(option("--sites")) && integer_option("--sites", sites, 2L) != sites) {
    stop("--sites must equal the number of --site-sizes.")
  }
}
taus <- as.numeric(strsplit(option("--taus", "0.1"), ",", fixed = TRUE)[[1L]])
seed <- integer_option("--seed", 667L, 0L)
replications <- integer_option("--replications", 1L)
reverse_sites <- boolean_option("--reverse-sites")
compute_intervals <- boolean_option("--intervals")
output <- option("--output", file.path(root, "results", "example.csv"))
if (!length(taus) || anyNA(taus) || any(taus <= 0 | taus >= 1) ||
    anyDuplicated(taus) || as.double(seed) + replications - 1 > .Machine$integer.max) {
  stop("Invalid quantile levels or seed range.")
}

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
