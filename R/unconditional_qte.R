scqte_unconditional <- function(data, sizes, tau,
                                h_two_round = scqte_bandwidth(nrow(data)),
                                h_one_round = 0.05,
                                compute_intervals = TRUE) {
  scqte_validate(data, sizes, tau)
  if (!is.finite(h_two_round) || h_two_round <= 0 ||
      !is.finite(h_one_round) || h_one_round <= 0) stop("Bandwidths must be positive.")
  X <- cbind(1, data$x1, data$x2)
  A <- data$A
  Y <- data$Y
  central_data <- scqte_weighted_data(X, A, Y,
                                      scqte_central_propensity(X, A))
  sequential_data <- scqte_weighted_data(X, A, Y,
                                         scqte_renewable_propensity(X, A, sizes))

  ip <- quantreg::rq.fit.fnb(central_data$X, central_data$Y,
                             tau = tau, eps = 1e-6)$coefficients
  central <- SQRcppW(central_data$X, central_data$Y, tau = tau,
                     h = h_two_round)
  fedavg <- scqte_fedavg(central_data$X, central_data$Y, sizes, tau)
  two_round <- seqSQRcppW(sequential_data$X, sequential_data$Y, sizes,
                           tau = tau, h = h_two_round)
  one_round <- seqSQRcppW1R(X, A, Y, sizes, rep(0, ncol(X)),
                             tau = tau, h = h_one_round)

  interval <- function(coefficient) if (compute_intervals)
    scqte_interval(central_data$X, central_data$Y, coefficient,
                   tau, h_two_round, 2L) else c(NA_real_, NA_real_)
  rbind(
    scqte_result("IP", ip[2L]),
    scqte_result("conquer", central$coeff[2L],
                 interval(central$coeff), central$ite),
    scqte_result("FedAvg", fedavg[2L]),
    scqte_result("SCQTE-1R", one_round$coeff[2L],
                 interval(one_round$coeff), tail(one_round$ite, 1L)),
    scqte_result("SCQTE-2R", two_round$coeff[2L],
                 interval(two_round$coeff), tail(two_round$ite, 1L))
  )
}
