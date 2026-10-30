# Gauss-Legendre nodes and weights on [0, 1] for covariate integration.
scqte_quadrature <- function(n) {
  k <- seq_len(n - 1L)
  off_diagonal <- k / sqrt(4 * k^2 - 1)
  matrix <- matrix(0, n, n)
  matrix[cbind(k, k + 1L)] <- off_diagonal
  matrix[cbind(k + 1L, k)] <- off_diagonal
  eig <- eigen(matrix, symmetric = TRUE)
  indices <- order(eig$values)
  list(nodes = (eig$values[indices] + 1) / 2,
       weights = eig$vectors[1L, indices]^2)
}

# Evaluate each marginal potential-outcome CDF on a fixed quadrature grid.
scqte_truth_at_nodes <- function(covariates, error, sizes, taus, n_nodes) {
  quadrature <- scqte_quadrature(n_nodes)
  u <- quadrature$nodes
  if (covariates == "shift") {
    # Average the site-specific joint densities using their population shares.
    t <- seq(0, 1, length.out = length(sizes))
    density1 <- vapply(t, function(v) dbeta(u, 1 + 4 * v, 5 - 4 * v), numeric(n_nodes))
    density2 <- vapply(t, function(v) dbeta(u, 5 - 4 * v, 1 + 4 * v), numeric(n_nodes))
    joint_density <- tcrossprod(sweep(density1, 2L, sizes / sum(sizes), `*`),
                               density2)
  } else {
    joint_density <- matrix(1, n_nodes, n_nodes)
  }
  mass <- joint_density * outer(quadrature$weights, quadrature$weights)
  mass <- mass / sum(mass)
  x1 <- matrix(rep(u, n_nodes), nrow = n_nodes)
  x2 <- matrix(rep(2 * u, each = n_nodes), nrow = n_nodes)
  quantiles <- matrix(NA_real_, length(taus), 2L)
  for (a in 0:1) {
    location <- 1 + a + x1 + x2
    scale <- 1 + x1 + 1.5 * a
    # Integrate the conditional outcome CDF to get each marginal CDF.
    cdf <- function(q) {
      standardized <- (q - location) / scale
      if (error == "normal") {
        probability <- pnorm(standardized)
      } else {
        probability <- matrix(0, n_nodes, n_nodes)
        above_support <- standardized >= 1
        probability[above_support] <- 1 - standardized[above_support]^(-5)
      }
      sum(mass * probability)
    }
    for (j in seq_along(taus)) {
      # Invert each marginal CDF separately, then take their difference.
      quantiles[j, a + 1L] <- uniroot(function(q) cdf(q) - taus[j],
                                      lower = -100, upper = 1000,
                                      tol = 1e-10)$root
    }
  }
  quantiles[, 2L] - quantiles[, 1L]
}

# Return the difference between the two marginal potential-outcome quantiles.
scqte_population_qte <- function(covariates, error, sizes, taus,
                                  tolerance = 2e-6) {
  if (!covariates %in% c("uniform", "shift") ||
      !error %in% c("normal", "pareto") ||
      !is.numeric(sizes) || length(sizes) < 2L ||
      any(!is.finite(sizes)) || any(sizes < 1L) ||
      any(!is.finite(taus)) || any(taus <= 0 | taus >= 1)) {
    stop("Invalid population-target arguments.")
  }
  previous <- NULL
  # Refine the quadrature until the treatment contrast stabilizes.
  for (n_nodes in c(64L, 128L, 256L, 512L)) {
    current <- scqte_truth_at_nodes(covariates, error, sizes, taus, n_nodes)
    if (!is.null(previous) && max(abs(current - previous)) < tolerance) {
      return(current)
    }
    previous <- current
  }
  stop("Population quantile calculation did not converge at 512 quadrature nodes.")
}
