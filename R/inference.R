# Joint propensity and smoothed-quantile score summaries for one site.
scqte_joint_site_summary <- function(X, A, Y, gamma, beta, tau, h) {
  propensity <- plogis(drop(X %*% gamma))
  weight <- A / propensity + (1 - A) / (1 - propensity)
  z <- cbind(1, A)
  t <- weight * drop(z %*% beta - Y) / h
  density <- dnorm(t)
  quantile_score <- pnorm(t) - tau
  score_gamma <- X * (A - propensity)
  score_beta <- z * (weight * quantile_score)
  list(
    G = crossprod(X, X * (propensity * (1 - propensity))),
    C = crossprod(z, X * (weight * (propensity - A) *
                          (quantile_score + t * density))),
    D = crossprod(z, z * (weight^2 * density / h)),
    Bgg = crossprod(score_gamma),
    Bqg = crossprod(score_beta, score_gamma),
    Bqq = crossprod(score_beta)
  )
}

# Two-round inference uses the final renewable propensity fit and the
# quantile estimate available when each site was visited in round two.
scqte_2r_joint_interval <- function(X, A, Y, sizes, gamma, fit, tau, h) {
  H <- matrix(0, 5L, 5L)
  B <- matrix(0, 5L, 5L)
  ends <- cumsum(sizes)
  starts <- c(1L, head(ends, -1L) + 1L)
  for (site in seq_along(sizes)) {
    rows <- seq.int(starts[site], ends[site])
    part <- scqte_joint_site_summary(X[rows, , drop = FALSE], A[rows], Y[rows],
                                  gamma, fit$record[, site], tau, h)
    H[1:3, 1:3] <- H[1:3, 1:3] + part$G
    H[4:5, 1:3] <- H[4:5, 1:3] - part$C
    H[4:5, 4:5] <- H[4:5, 4:5] - part$D
    B[1:3, 1:3] <- B[1:3, 1:3] + part$Bgg
    B[4:5, 1:3] <- B[4:5, 1:3] + part$Bqg
    B[1:3, 4:5] <- B[1:3, 4:5] + t(part$Bqg)
    B[4:5, 4:5] <- B[4:5, 4:5] + part$Bqq
  }
  inverse <- solve(H)
  se <- sqrt((inverse %*% B %*% t(inverse))[5L, 5L])
  half_width <- qnorm(.975) * se
  c(fit$coeff[2L] - half_width, fit$coeff[2L] + half_width)
}
