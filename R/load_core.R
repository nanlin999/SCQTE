scqte_load_core <- function(root = ".",
                            macos_linker_fix = Sys.info()[["sysname"]] == "Darwin") {
  for (package in c("Rcpp", "RcppArmadillo", "quantreg")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Install the R package '", package, "' before running this example.")
    }
  }
  cpp <- file.path(root, "src", "scqte_core.cpp")
  if (!file.exists(cpp)) stop("Cannot find ", cpp)
  if (macos_linker_fix) {
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

scqte_validate <- function(data, sizes, tau) {
  if (!is.data.frame(data) || !all(c("x1", "x2", "A", "Y") %in% names(data))) {
    stop("data must contain x1, x2, A, and Y columns.")
  }
  if (anyNA(data[c("x1", "x2", "A", "Y")]) ||
      any(!is.finite(as.matrix(data[c("x1", "x2", "A", "Y")])))) {
    stop("Data must be finite and complete.")
  }
  if (!all(data$A %in% c(0, 1))) stop("A must be binary.")
  if (!is.numeric(sizes) || length(sizes) < 2L ||
      any(!is.finite(sizes)) || any(sizes < 1) ||
      any(sizes != as.integer(sizes)) || sum(sizes) != nrow(data)) {
    stop("sizes must be positive integer block sizes summing to nrow(data).")
  }
  if (length(tau) != 1L || !is.finite(tau) || tau <= 0 || tau >= 1) {
    stop("tau must be a single number strictly between zero and one.")
  }
  invisible(TRUE)
}

scqte_bandwidth <- function(n, p = 2L) ((p + log(n)) / n)^0.4

scqte_interval <- function(X, Y, coefficient, tau, h, index) {
  limits <- asymptoticCI(X, drop(Y - X %*% coefficient), coefficient,
                         tau, nrow(X), h, qnorm(0.975))
  unname(limits[index, ])
}

scqte_result <- function(method, estimate, limits = c(NA_real_, NA_real_),
                         final_iterations = NA_integer_) {
  data.frame(method = method, estimate = unname(estimate),
             lower = limits[1L], upper = limits[2L],
             final_iterations = final_iterations,
             converged = if (is.na(final_iterations)) NA else
               final_iterations <= 5000L)
}
