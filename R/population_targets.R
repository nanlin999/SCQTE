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

scqte_truth_at_nodes <- function(design, error, sites, taus, n_nodes) {
  quadrature <- scqte_quadrature(n_nodes)
  u <- quadrature$nodes
  if (design == "shift") {
    t <- seq(0, 1, length.out = sites)
    density1 <- vapply(t, function(v) dbeta(u, 1 + 4 * v, 5 - 4 * v), numeric(n_nodes))
    density2 <- vapply(t, function(v) dbeta(u, 5 - 4 * v, 1 + 4 * v), numeric(n_nodes))
    joint_density <- tcrossprod(density1, density2) / sites
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
      quantiles[j, a + 1L] <- uniroot(function(q) cdf(q) - taus[j],
                                      lower = -100, upper = 1000,
                                      tol = 1e-10)$root
    }
  }
  quantiles[, 2L] - quantiles[, 1L]
}

scqte_population_qte <- function(design, error, sites, taus,
                                  tolerance = 2e-6) {
  if (!design %in% c("equal", "shift", "unequal") ||
      !error %in% c("normal", "pareto") ||
      any(!is.finite(taus)) || any(taus <= 0 | taus >= 1)) {
    stop("Invalid population-target arguments.")
  }
  previous <- NULL
  for (n_nodes in c(64L, 128L, 256L, 512L)) {
    current <- scqte_truth_at_nodes(design, error, sites, taus, n_nodes)
    if (!is.null(previous) && max(abs(current - previous)) < tolerance) {
      return(current)
    }
    previous <- current
  }
  stop("Population quantile calculation did not converge at 512 quadrature nodes.")
}
