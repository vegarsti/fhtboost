FHT_minus_loglikelihood_with_all_parameters <- function(optimization_parameters, X, Z, times, delta) {
  p <- dim(X)[2]
  d <- dim(Z)[2]
  parameter_list <- parameter_vector_to_list(optimization_parameters, p, d)
  beta_ <- parameter_list$beta
  gamma_ <- parameter_list$gamma
  y0 <- exp(X %*% beta_)
  mu <- Z %*% gamma_
  log_f <- log(y0) - 0.5*log(2*pi*times^3) - (y0 + mu*times)^2/(2*times)
  log_S <- log(pnorm((mu*times + y0)/sqrt(times)) - exp(-2*y0*mu)*pnorm((mu*times - y0))/sqrt(times))
  loglikelihood_value <- sum(delta*log_f) + sum((1-delta)*log_S)
  return(-loglikelihood_value)
}
