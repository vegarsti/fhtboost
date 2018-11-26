boosting_iteration_both <- function(nu, X, Z, u_y0, u_mu, beta_hat_m1, gamma_hat_m1, d, ds, p, ps, times, delta) {
  # d corresponds to X, p to Z
  y0 <- exp(X %*% beta_hat_m1)
  result_y0 <- best_least_squares_update(X, u_y0, d, ds)
  result_mu <- best_least_squares_update(Z, u_mu, p, ps)
  rsses <- c(result_y0$rss, result_mu$rss)
  print(rsses)
  best_rss <- which.min(rsses)
  boosted_mu <- best_rss == 2
  if (boosted_mu) {
    # mu; gamma
    gamma_hat_addition <- nu*result_mu$parameter_updates
    gamma_hat_m <- gamma_hat_m1 + gamma_hat_addition
    beta_hat_addition <- 0
    beta_hat_m <- beta_hat_m1
  } else {
    # y0; beta
    beta_hat_addition <- nu*result_y0$parameter_updates
    beta_hat_m <- beta_hat_m1 + beta_hat_addition
    gamma_hat_addition <- 0
    gamma_hat_m <- gamma_hat_m1
  }
  return(list(
    beta_hat_m=beta_hat_m, beta_hat_addition=beta_hat_addition,
    gamma_hat_m=gamma_hat_m, gamma_hat_addition=gamma_hat_addition,
    boosted_mu=boosted_mu, rsses=rsses
  ))
}
