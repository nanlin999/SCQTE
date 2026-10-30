#!/usr/bin/env Rscript

# Run a selected simulation design and write method-level results.
arguments <- commandArgs(trailingOnly = TRUE)
option <- function(name, default) {
  hit <- which(arguments == name)
  if (!length(hit)) return(default)
  if (length(hit) != 1L || hit == length(arguments)) stop("Expected one value after ", name)
  arguments[hit + 1L]
}
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.")
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."))
for (file in c("load_core.R", "benchmarks.R", "propensity.R",
               "conditional_qte.R", "unconditional_qte.R",
               "simulation_scenarios.R", "population_targets.R")) {
  source(file.path(root, "R", file))
}

design <- option("--design", "equal")
error <- option("--error", "normal")
sites <- as.integer(option("--sites", if (design == "unequal") "200" else "20"))
pattern <- option("--pattern", if (design == "unequal") "60-40" else "equal")
taus <- as.numeric(strsplit(option("--taus", "0.1"), ",", fixed = TRUE)[[1L]])
first_replication <- as.integer(option("--first-replication", "1"))
replications <- as.integer(option("--replications", "1"))
output <- option("--output", file.path(root, "results", "example.csv"))
if (anyNA(c(sites, taus, first_replication, replications)) ||
    replications < 1L || first_replication < 1L ||
    any(taus <= 0 | taus >= 1) || anyDuplicated(taus)) {
  stop("Invalid sites, tau levels, or replication range.")
}

scqte_load_core(root)
truth_unconditional <- scqte_population_qte(design, error, sites, taus)
truth_conditional <- 1 + 1.5 * if (error == "normal") qnorm(taus) else
  (1 - taus)^(-1 / 5)

rows <- list()
for (replication in seq.int(first_replication,
                            length.out = replications)) {
  scenario <- scqte_scenario(design, error, sites, pattern, replication)
  for (j in seq_along(taus)) {
    tau <- taus[j]
    for (estimand in c("conditional", "unconditional")) {
      fit <- if (estimand == "conditional")
        scqte_conditional(scenario$data, scenario$sizes, tau,
                          compute_intervals = design == "equal") else
        scqte_unconditional(scenario$data, scenario$sizes, tau,
                            compute_intervals = design == "equal")
      truth <- if (estimand == "conditional") truth_conditional[j] else
        truth_unconditional[j]
      fit$design <- design
      fit$error <- error
      fit$sites <- sites
      fit$pattern <- pattern
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
results <- results[c("design", "error", "sites", "pattern", "replication",
                     "tau", "estimand", "method", "estimate", "population_qte",
                     "absolute_error", "lower", "upper", "covered",
                     "final_iterations", "converged")]
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write.csv(results, output, row.names = FALSE)
cat("Wrote ", nrow(results), " method-level results to ",
    normalizePath(output), "\n", sep = "")
