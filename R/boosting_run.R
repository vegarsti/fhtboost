boosting_run <- function(times, delta, X, Z, m_stop, nlm_parameter_list, boost_mu=TRUE) {
  N <- dim(X)[1]
  N <- dim(Z)[1]
  nu <- 0.1
  d <- dim(X)[2]
  p <- dim(Z)[2]

  beta_from_nlm <- nlm_parameter_list$beta
  gamma_from_nlm <- nlm_parameter_list$gamma

  # Scaling X and Z -- should rewrite! or hide in function
  X_scale_factors <- as.numeric(apply(X, 2, sd))
  Z_scale_factors <- as.numeric(apply(Z, 2, sd))
  X_scale_factors[1] <- 1
  Z_scale_factors[1] <- 1
  for (j in 2:d) {
    X[, j] <- X[, j]/X_scale_factors[j]
  }
  for (j in 2:p) {
    Z[, j] <- Z[, j]/Z_scale_factors[j]
  }

  # "make" loss function after scaling Z
  minus_FHT_loglikelihood <- data_to_FHT_loss_function(X, Z, times, delta)
  #minus_FHT_loglikelihood(list(beta=beta_true, gamma=gamma_true))

  gamma_hat <- matrix(NA, nrow=m_stop, ncol=p)
  gamma_hat[1, ] <- rep(0, p)
  gamma_hat_cumsum <- matrix(NA, nrow=m_stop, ncol=p)

  beta_hat <- matrix(NA, nrow=m_stop, ncol=d)
  beta_hat[1, ] <- rep(0, d)
  beta_hat_cumsum <- matrix(NA, nrow=m_stop, ncol=d)

  # INITIALIZE
  if (mu_has_intercept) {
    gamma_hat[1, 1] <- gamma_from_nlm[1] # true value -- must estimate better (unseen)
  }
  if (y0_has_intercept) {
    beta_hat[1, 1] <- beta_from_nlm[1] # true value -- must estimate better (unseen)
  }
  gamma_hat_cumsum[1, ] <- gamma_hat[1, ]
  beta_hat_cumsum[1, ] <- beta_hat[1, ]


  negative_gradient <- matrix(NA, nrow=m_stop, ncol=N)
  negative_gradient[1, ] <- FHT_componentwise_loss_function_derivative_mu(beta_from_nlm, gamma_hat[1, ], X, Z, times, delta)
  #negative_gradient[1, ] <- FHT_componentwise_loss_function_derivative_mu(beta_hat[1, ], gamma_from_nlm, X, Z, times, delta)
  loss <- rep(NA, m_stop)
  if (boost_mu) { parameter_list <- list(beta=beta_from_nlm, gamma_hat[1, ]) }
  else { parameter_list <- list(beta=beta_hat[1, ], gamma_from_nlm) }
  #loss[1] <- minus_FHT_loglikelihood(parameter_list)

  gamma_hats <- rep(NA, p)
  u_hats_rss <- rep(NA, p)
  u_hats <- matrix(NA, nrow=p, ncol=N)
  beta_hats <- rep(NA, d)


  if (mu_has_intercept) {
    ps <- 2:p
  } else {
    ps <- 1:p
  }
  if (y0_has_intercept) {
    ds <- 2:d
  } else {
    ds <- 1:p
  }

  for (m in 2:m_stop) {
    u <- negative_gradient[(m-1), ]
    ## BOOST MU
    if (boost_mu) {
      result <- boosting_iteration_mu(
        nu, X, Z, u, beta_from_nlm, p, ps, gamma_hat_cumsum[m-1, ], times, delta
      )
      gamma_hat[m, ] <- result$gamma_hat_addition
      gamma_hat_cumsum[m, ] <- result$gamma_hat_m
      negative_gradient[m, ] <- result$u_m
      loss[m] <- result$loss_m
    } else {
      result <- boosting_iteration_y0(
        nu, X, Z, u, gamma_from_nlm, d, ds, beta_hat_cumsum[m-1, ], times, delta
      )
      beta_hat[m, ] <- result$beta_hat_addition
      beta_hat_cumsum[m, ] <- result$beta_hat_m
      negative_gradient[m, ] <- result$u_m
      loss[m] <- result$loss_m
    }
  }
  if (boost_mu) {
    parameters <- gamma_hat_cumsum
    final_parameters <- gamma_hat_final <- gamma_hat_cumsum[m, ] * Z_scale_factors
  }
  else {
    final_parameters <- beta_hat_final <- beta_hat_cumsum[m, ] * X_scale_factors
    parameters <- beta_hat_cumsum
  }
  gradient_sums <- rowSums(abs(negative_gradient))
  return(list(
    final_parameters=final_parameters,
    gradient_sums=gradient_sums,
    parameters=parameters
  ))
}
