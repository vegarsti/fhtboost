boosting_iteration_mu <- function(nu, X, Z, u, beta_, p, ps, gamma_hat_m1, times, delta) {
  gamma_hat_addition <- nu*best_least_squares_update(Z, u, p, ps)
  gamma_hat_m <- gamma_hat_m1 + gamma_hat_addition
  u_m <- FHT_componentwise_loss_function_derivative_mu(
    beta_, gamma_hat_m, X, Z, times, delta
  )
  loss_m <- FHT_minus_loglikelihood_with_all_parameters(
    beta_, gamma_hat_m, X, Z, times, delta
  )
  return(list(u_m=u_m, loss_m=loss_m, gamma_hat_m=gamma_hat_m, gamma_hat_addition=gamma_hat_addition))
}
