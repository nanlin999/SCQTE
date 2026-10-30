#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.")
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."))
for (file in c("load_core.R", "benchmarks.R", "propensity.R",
               "conditional_qte.R", "unconditional_qte.R",
               "simulation_scenarios.R", "population_targets.R")) {
  source(file.path(root, "R", file))
}
scqte_load_core(root)

assert_close <- function(actual, expected, tolerance = 1e-8) {
  if (length(actual) != length(expected) || anyNA(actual) ||
      max(abs(actual - expected)) > tolerance) {
    stop("Reference mismatch: expected ", paste(expected, collapse = ", "),
         "; got ", paste(actual, collapse = ", "))
  }
}

check_fit <- function(fit, expected) {
  actual <- setNames(fit$estimate, fit$method)
  assert_close(actual[names(expected)], expected)
}

# Reference values for the supported example designs.
equal_normal <- scqte_scenario(rep(500L, 20L), "normal", "uniform", 667L)
stopifnot(sum(equal_normal$sizes) == 10000L)
conditional <- scqte_conditional(equal_normal$data, equal_normal$sizes, 0.1)
unconditional <- scqte_unconditional(equal_normal$data, equal_normal$sizes, 0.1)
check_fit(conditional, c(IP = -0.940271962058658,
                         conquer = -0.929633311274172,
                         FedAvg = -0.931560314860604,
                         `SCQTE-1R` = -0.932382273439628))
check_fit(unconditional, c(IP = -0.832956959153672,
                           conquer = -0.827323489153382,
                           FedAvg = -0.856375474942910,
                           `SCQTE-1R` = -0.838612173776605,
                           `SCQTE-2R` = -0.839772933276069))
assert_close(conditional$lower[conditional$method == "SCQTE-1R"],
             -1.0822884623499)
assert_close(conditional$upper[conditional$method == "SCQTE-1R"],
             -0.782476084529359)
assert_close(unconditional$lower[unconditional$method == "SCQTE-1R"],
             -1.01683649084197)
assert_close(unconditional$upper[unconditional$method == "SCQTE-1R"],
             -0.660387856711242)
assert_close(scqte_population_qte("uniform", "normal", rep(500L, 20L), 0.1),
             -0.851324555137758, 2e-6)

equal_pareto <- scqte_scenario(rep(200L, 50L), "pareto", "uniform", 100667L)
conditional <- scqte_conditional(equal_pareto$data, equal_pareto$sizes, 0.9)
unconditional <- scqte_unconditional(equal_pareto$data, equal_pareto$sizes, 0.9)
check_fit(conditional, c(IP = 3.32885589145778,
                         conquer = 3.32835754172834,
                         FedAvg = 3.31642875363058,
                         `SCQTE-1R` = 3.33013913492285))
check_fit(unconditional, c(IP = 3.03948429240956,
                           conquer = 3.04885653784949,
                           FedAvg = 3.09636389578865,
                           `SCQTE-1R` = 3.06611220735245,
                           `SCQTE-2R` = 3.04905614068988))
assert_close(scqte_population_qte("uniform", "pareto", rep(200L, 50L), 0.9),
             3.0857796051125, 2e-6)

shift <- scqte_scenario(rep(500L, 20L), "normal", "shift", 200715L)
conditional <- scqte_conditional(shift$data, shift$sizes, 0.6,
                                  compute_intervals = FALSE)
unconditional <- scqte_unconditional(shift$data, shift$sizes, 0.6,
                                      compute_intervals = FALSE)
check_fit(conditional, c(IP = 1.44811508976792,
                         conquer = 1.45199326007178,
                         FedAvg = 1.44068877738743,
                         `SCQTE-1R` = 1.45213884740012))
check_fit(unconditional, c(IP = 1.41189493895856,
                           conquer = 1.40903682227767,
                           FedAvg = 1.42148033051935,
                           `SCQTE-1R` = 28.7992059200958,
                           `SCQTE-2R` = 1.46980436030453))
unstable <- unconditional[unconditional$method == "SCQTE-1R", ]
stopifnot(identical(unstable$converged, FALSE),
          unstable$final_iterations == 5001L,
          all(is.na(unconditional$lower)))

unequal <- scqte_scenario(c(rep(70L, 100L), rep(30L, 100L)),
                          "normal", "uniform", 667L, reverse_sites = TRUE)
stopifnot(length(unequal$sizes) == 200L,
          identical(unequal$sizes[1L], 30L),
          identical(unequal$sizes[200L], 70L))
conditional <- scqte_conditional(unequal$data, unequal$sizes, 0.1,
                                  compute_intervals = FALSE)
unconditional <- scqte_unconditional(unequal$data, unequal$sizes, 0.1,
                                      compute_intervals = FALSE)
check_fit(conditional, c(FedAvg = -0.760241055155428,
                         `SCQTE-1R` = -0.896035841712526))
check_fit(unconditional, c(FedAvg = -0.934643691895402,
                           `SCQTE-1R` = -0.767495055492192,
                           `SCQTE-2R` = -0.781591367583534))
cat("Public-code regression tests passed for equal-normal, equal-Pareto, ",
    "covariate-shift, and unequal-site designs.\n", sep = "")
