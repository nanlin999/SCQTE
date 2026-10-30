# Reorder complete site blocks when comparing different arrival orders.
scqte_reverse_sites <- function(data, sizes) {
  ends <- cumsum(sizes)
  starts <- c(1L, head(ends, -1L) + 1L)
  rows <- unlist(rev(Map(seq.int, starts, ends)), use.names = FALSE)
  list(data = data[rows, , drop = FALSE], sizes = rev(sizes))
}

# Generate one outcome-model dataset with contiguous, user-specified sites.
scqte_scenario <- function(sizes, error = "normal", covariates = "uniform",
                           seed = 667L, reverse_sites = FALSE) {
  if (!is.numeric(sizes) || length(sizes) < 2L ||
      any(!is.finite(sizes)) || any(sizes < 1L) ||
      any(sizes != as.integer(sizes)) ||
      !error %in% c("normal", "pareto") ||
      !covariates %in% c("uniform", "shift") ||
      length(seed) != 1L || !is.finite(seed) ||
      seed < 0L || seed != as.integer(seed) ||
      length(reverse_sites) != 1L || is.na(reverse_sites) ||
      !is.logical(reverse_sites)) {
    stop("Invalid simulation parameters.")
  }
  sizes <- as.integer(sizes)
  n <- sum(sizes)
  set.seed(as.integer(seed))

  if (covariates == "uniform") {
    x1 <- runif(n)
    x2 <- runif(n, 0, 2)
  } else {
    # Site-specific beta distributions produce covariate shift across sites.
    position <- seq(0, 1, length.out = length(sizes))
    x1 <- unlist(Map(function(m, t) rbeta(m, 1 + 4 * t, 5 - 4 * t),
                     sizes, position), use.names = FALSE)
    x2 <- 2 * unlist(Map(function(m, t) rbeta(m, 5 - 4 * t, 1 + 4 * t),
                         sizes, position), use.names = FALSE)
  }
  A <- rbinom(n, 1, plogis(-0.5 + 0.5 * x1 + 0.5 * x2))
  epsilon <- if (error == "normal") rnorm(n) else runif(n)^(-1 / 5)
  data <- data.frame(x1 = x1, x2 = x2, A = A,
                     Y = 1 + A + x1 + x2 + (1 + x1 + 1.5 * A) * epsilon)

  if (reverse_sites) {
    # Reverse complete site blocks, leaving observations within each block intact.
    ordered <- scqte_reverse_sites(data, sizes)
    data <- ordered$data
    sizes <- ordered$sizes
  }
  list(data = data, sizes = sizes)
}
