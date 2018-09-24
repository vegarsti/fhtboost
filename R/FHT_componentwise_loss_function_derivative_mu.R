FHT_componentwise_loss_function_derivative_mu <- function(parameter_list, X, Z, times, delta) {
  p <- dim(X)[2]
  d <- dim(Z)[2]
  y0 <- exp(X %*% parameter_list$beta)
  mu <- Z %*% parameter_list$gamma
  sigma2 <- 1
  gradient <- - loss_function_derivative_mu(y0, mu, sigma2, times, delta)
  return(gradient)
}
