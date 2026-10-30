# Compile and expose the C++ optimization routines to the R example.
scqte_load_core <- function(root = ".",
                            macos_linker_fix = Sys.info()[["sysname"]] == "Darwin") {
  cpp <- file.path(root, "src", "scqte_core.cpp")
  if (macos_linker_fix) {
    # A temporary Makevars avoids stale macOS Fortran linker paths while
    # preserving the user's compiler settings outside this compilation.
    previous <- Sys.getenv("R_MAKEVARS_USER", unset = NA_character_)
    existing <- if (is.na(previous)) path.expand("~/.R/Makevars") else previous
    temporary <- tempfile("scqte-Makevars-")
    writeLines(c(if (file.exists(existing)) readLines(existing, warn = FALSE),
                 "FLIBS ="), temporary)
    Sys.setenv(R_MAKEVARS_USER = temporary)
    on.exit({
      if (is.na(previous)) Sys.unsetenv("R_MAKEVARS_USER") else
        Sys.setenv(R_MAKEVARS_USER = previous)
      unlink(temporary)
    }, add = TRUE)
  }
  Rcpp::sourceCpp(cpp, env = .GlobalEnv)
  invisible(cpp)
}

# Smoothing bandwidth used for the pooled and two-round fits.
scqte_bandwidth <- function(n, p = 2L) ((p + log(n)) / n)^0.4

scqte_interval <- function(X, Y, coefficient, tau, h, index) {
  # Return only the coefficient interval for the requested estimand.
  limits <- asymptoticCI(X, drop(Y - X %*% coefficient), coefficient,
                         tau, nrow(X), h, qnorm(0.975))
  unname(limits[index, ])
}

scqte_result <- function(method, estimate, limits = c(NA_real_, NA_real_),
                         final_iterations = NA_integer_) {
  # The C++ solvers report 5001 iterations when the 5000-step limit is hit.
  data.frame(method = method, estimate = unname(estimate),
             lower = limits[1L], upper = limits[2L],
             final_iterations = final_iterations,
             converged = if (is.na(final_iterations)) NA else
               final_iterations <= 5000L)
}
