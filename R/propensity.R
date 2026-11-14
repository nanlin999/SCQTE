# Construct the inverse-propensity-weighted quantile-regression data.
scqte_weighted_data <- function(X, A, Y, propensity_coefficient) {
  # The weighted two-column design implements marginal QTE as a coefficient
  # contrast: column 1 is the intercept, column 2 is the treatment indicator.
  probability <- plogis(drop(X %*% propensity_coefficient))
  weights <- A / probability + (1 - A) / (1 - probability)
  # Weight both the design and response for the C++ weighted-QR routines.
  list(X = cbind(weights, weights * A), Y = weights * Y,
       weights = weights)
}

scqte_central_propensity <- function(X, A) {
  # Logistic score solved once using pooled observations.
  as.numeric(Newtoncpp(X, A, rep(0, ncol(X))))
}

scqte_renewable_propensity <- function(X, A, sizes) {
  # Sequential Newton update reuses accumulated site information.
  as.numeric(seqNewtoncpp(X, A, rep(0, ncol(X)), sizes))
}

scqte_fedavg_propensity <- function(X, A, sizes) {
  hessian <- matrix(0, ncol(X), ncol(X))
  rhs <- numeric(ncol(X))
  ends <- cumsum(sizes)
  starts <- c(1L, head(ends, -1L) + 1L)
  for (site in seq_along(sizes)) {
    rows <- seq.int(starts[site], ends[site])
    Xi <- X[rows, , drop = FALSE]
    beta <- Newtoncpp(Xi, A[rows], rep(0, ncol(X)))
    Hi <- Loghess(Xi, beta)
    hessian <- hessian + Hi
    rhs <- rhs + drop(Hi %*% beta)
  }
  drop(solve(hessian, rhs))
}
