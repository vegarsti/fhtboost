FHT_componentwise_loss_function_derivative_y0 <- function(beta_, gamma_, X, Z, times, delta) {
  p <- dim(X)[2]
  d <- dim(Z)[2]
  y0 <- exp(X %*% beta_)
  #y0 <- exp(X %*% beta_) / (exp(X %*% beta_) + 1) # <-- logistic
  mu <- Z %*% gamma_
  sigma2 <- 1
  negative_gradient <- loss_function_derivative_y0(y0, mu, sigma2, times, delta)
  return(negative_gradient)
}
