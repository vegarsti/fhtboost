FHT_componentwise_minus_loglikelihood_with_parameters <- function(beta_, gamma_, X, Z, times, delta) {
  p <- dim(X)[2]
  d <- dim(Z)[2]
  y0 <- exp(X %*% beta_)
  mu <- Z %*% gamma_
  loglikelihood_vector <- FHT_loglikelihood_with_y0_mu(y0, mu, times, delta)
  return(-loglikelihood_vector)
}
