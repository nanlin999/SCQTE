scqte_weighted_data <- function(X, A, Y, propensity_coefficient) {
  probability <- plogis(drop(X %*% propensity_coefficient))
  if (any(!is.finite(probability)) || any(probability <= 0 | probability >= 1)) {
    stop("Propensity estimates must be finite and strictly between zero and one.")
  }
  weights <- A / probability + (1 - A) / (1 - probability)
  list(X = cbind(weights, weights * A), Y = weights * Y,
       weights = weights)
}

scqte_central_propensity <- function(X, A) {
  as.numeric(Newtoncpp(X, A, rep(0, ncol(X))))
}

scqte_renewable_propensity <- function(X, A, sizes) {
  as.numeric(seqNewtoncpp(X, A, rep(0, ncol(X)), sizes))
}
