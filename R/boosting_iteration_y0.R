boosting_iteration_y0 <- function(nu, X, Z, u, beta_hat_m1, gamma_hat_m1, d, ds, times, delta) {
  beta_hat_addition <- nu*best_least_squares_update(X, u, d, ds)
  beta_hat_m <- beta_hat_m1 + beta_hat_addition
  u_m <- FHT_componentwise_loss_function_derivative_y0(
    beta_hat_m, gamma_hat_m1, X, Z, times, delta
  )
  return(list(u_m=u_m, beta_hat_m=beta_hat_m, beta_hat_addition=beta_hat_addition))
}
