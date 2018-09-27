FHT_componentwise_loss_function_derivative_mu <- function(beta_, gamma_, X, Z, times, delta) {
  p <- dim(X)[2]
  d <- dim(Z)[2]
  y0 <- exp(X %*% beta_)
  mu <- Z %*% gamma_
  sigma2 <- 1
  negative_gradient <- loss_function_derivative_mu(y0, mu, sigma2, times, delta)
  return(negative_gradient)
}
