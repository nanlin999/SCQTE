scqte_conditional <- function(data, sizes, tau,
                              h = scqte_bandwidth(nrow(data)),
                              compute_intervals = TRUE) {
  scqte_validate(data, sizes, tau)
  if (!is.finite(h) || h <= 0) stop("h must be positive.")
  X <- cbind(1, data$x1, data$x2, data$A)
  Y <- data$Y
  ip <- quantreg::rq.fit.fnb(X, Y, tau = tau, eps = 1e-6)$coefficients
  central <- SQRcpp(X[, -1L, drop = FALSE], Y, tau = tau, h = h)
  fedavg <- scqte_fedavg(X, Y, sizes, tau)
  renewable <- seqSQRcpp(X[, -1L, drop = FALSE], Y, sizes,
                         tau = tau, h = h)
  interval <- function(coefficient) if (compute_intervals)
    scqte_interval(X, Y, coefficient, tau, h, 4L) else c(NA_real_, NA_real_)
  rbind(
    scqte_result("IP", ip[4L]),
    scqte_result("conquer", central$coeff[4L],
                 interval(central$coeff), central$ite),
    scqte_result("FedAvg", fedavg[4L]),
    scqte_result("SCQTE-1R", renewable$coeff[4L],
                 interval(renewable$coeff),
                 tail(renewable$ite, 1L))
  )
}
