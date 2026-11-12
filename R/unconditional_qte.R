# Estimate the marginal QTE through inverse-propensity-weighted quantile fits.
scqte_unconditional <- function(data, sizes, tau,
                                h_two_round = 2 * scqte_bandwidth(nrow(data)),
                                h_one_round = 2 * scqte_bandwidth(nrow(data)),
                                compute_intervals = TRUE) {
  X <- cbind(1, data$x1, data$x2)
  A <- data$A
  Y <- data$Y
  # Central and renewable propensity fits define the corresponding IPW
  # pseudo-data; the one-round routine estimates propensity internally.
  central_data <- scqte_weighted_data(X, A, Y,
                                      scqte_central_propensity(X, A))
  sequential_data <- scqte_weighted_data(X, A, Y,
                                         scqte_renewable_propensity(X, A, sizes))
  fedavg_data <- scqte_weighted_data(X, A, Y,
                                    scqte_fedavg_propensity(X, A, sizes))

  ip <- quantreg::rq.fit.fnb(central_data$X, central_data$Y,
                             tau = tau, eps = 1e-6)$coefficients
  central <- SQRcppW(central_data$X, central_data$Y, tau = tau,
                     h = h_two_round)
  fedavg <- scqte_fedavg(fedavg_data$X, fedavg_data$Y, sizes, tau)
  two_round <- seqSQRcppW(sequential_data$X, sequential_data$Y, sizes,
                           tau = tau, h = h_two_round)
  one_round <- seqSQRcppW1R(X, A, Y, sizes, rep(0, ncol(X)),
                             tau = tau, h = h_one_round)

  # Column 2 of the weighted design is the marginal treatment contrast.
  # Intervals use the pooled weighted data and the two-round bandwidth.
  interval <- function(coefficient, bandwidth) if (compute_intervals)
    scqte_interval(central_data$X, central_data$Y, coefficient,
                   tau, bandwidth, 2L) else c(NA_real_, NA_real_)
  rbind(
    scqte_result("IP", ip[2L]),
    scqte_result("conquer", central$coeff[2L],
                 interval(central$coeff, h_two_round), central$ite),
    scqte_result("FedAvg", fedavg[2L]),
    scqte_result("SCQTE-1R", one_round$coeff[2L],
                 interval(one_round$coeff, h_one_round), tail(one_round$ite, 1L)),
    scqte_result("SCQTE-2R", two_round$coeff[2L],
                 interval(two_round$coeff, h_two_round), tail(two_round$ite, 1L))
  )
}
