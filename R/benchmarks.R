# Hessian-weighted aggregation of independent local quantile-regression fits.
scqte_fedavg <- function(X, Y, sizes, tau) {
  hessian <- matrix(0, ncol(X), ncol(X))
  rhs <- numeric(ncol(X))
  ends <- cumsum(sizes)
  starts <- c(1L, head(ends, -1L) + 1L)
  for (site in seq_along(sizes)) {
    rows <- seq.int(starts[site], ends[site])
    Xi <- X[rows, , drop = FALSE]
    beta <- quantreg::rq.fit.fnb(Xi, Y[rows], tau = tau,
                                  eps = 1e-6)$coefficients
    Hi <- crossprod(Xi)
    hessian <- hessian + Hi
    rhs <- rhs + drop(Hi %*% beta)
  }
  drop(solve(hessian, rhs))
}
