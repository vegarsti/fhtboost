#' @export

boosting_run <- function(times, delta, X, Z, m_stop, boost_intercepts_continually, should_print=FALSE) {
  N_X <- dim(X)[1]
  N_Z <- dim(Z)[1]
  N_t <- length(times)
  N_delta <- length(delta)
  dimensions <- c(N_X, N_Z, N_t, N_delta)
  all_dimensions_are_equal <- var(dimensions) == 0
  if (! all_dimensions_are_equal) {
    cat('Dimensions:\n')
    cat('X: ')
    cat(N_X)
    cat('Z: ')
    cat(N_Z)
    cat('times: ')
    cat(N_t)
    cat('delta: ')
    cat(N_delta)
    stop("Not all dimensions are equal!")
  } else { # all dimensions are equal
    N <- N_X
  }
  d <- dim(X)[2]
  p <- dim(Z)[2]

  ps <- 2:p
  ds <- 2:d

  # Nuisance parameter: shrinkage factor
  nu <- 0.1

  # "make" loss function after scaling Z
  minus_FHT_loglikelihood <- data_to_FHT_loss_function(X, Z, times, delta)

  # initialize with nlm
  optimize_beta_0 <- function(beta0, beta_, gamma_, X, Z, times, delta) {
    beta_[1] <- beta0
    return(FHT_minus_loglikelihood_with_all_parameters(beta_, gamma_, X, Z, times, delta))
  }

  optimize_gamma_0 <- function(gamma0, beta_, gamma_, X, Z, times, delta) {
    gamma_[1] <- gamma0
    return(FHT_minus_loglikelihood_with_all_parameters(beta_, gamma_, X, Z, times, delta))
  }

  nlm_result_estimate <- maximum_likelihood_intercepts(times, delta)
  # result <- c()
  # for (i in 1:100) {
  #
  #   result <- cbind(nlm_result_estimate, result)
  # }

  #nlm_result_estimate <- c(2, -1)

  beta_0 <- nlm_result_estimate[1]
  gamma_0 <- nlm_result_estimate[2]

  gamma_hat <- matrix(NA, nrow=m_stop, ncol=p)
  gamma_hat[1, ] <- rep(0, p)
  gamma_hat_cumsum <- matrix(NA, nrow=m_stop, ncol=p)

  beta_hat <- matrix(NA, nrow=m_stop, ncol=d)
  beta_hat[1, ] <- rep(0, d)
  beta_hat_cumsum <- matrix(NA, nrow=m_stop, ncol=d)

  # INITIALIZE
  gamma_hat[1, 1] <- gamma_0
  beta_hat[1, 1] <- beta_0
  gamma_hat_cumsum[1, ] <- gamma_hat[1, ]
  beta_hat_cumsum[1, ] <- beta_hat[1, ]

  negative_gradient_mu <- matrix(NA, nrow=m_stop, ncol=N)
  negative_gradient_y0 <- matrix(NA, nrow=m_stop, ncol=N)
  negative_gradient_mu[1, ] <- FHT_componentwise_loss_function_derivative_mu(beta_hat[1, ], gamma_hat[1, ], X, Z, times, delta)
  negative_gradient_y0[1, ] <- FHT_componentwise_loss_function_derivative_y0(beta_hat[1, ], gamma_hat[1, ], X, Z, times, delta)
  loss <- rep(NA, m_stop)
  parameter_list <- list(beta=beta_hat[1, ], gamma_hat[1, ])

  gamma_hats <- rep(NA, p)
  beta_hats <- rep(NA, d)
  u_hats <- matrix(NA, nrow=p, ncol=N)
  u_hats_rss <- rep(NA, p)
  loss <- rep(0, m_stop)

  for (m in 2:m_stop) {
    # gamma/mu
    if (should_print) { cat('iteration: ', m, '\n') }
    u_y0 <- negative_gradient_y0[m-1, ]
    u_mu <- negative_gradient_mu[m-1, ]
    result <- boosting_iteration_both(
      nu, X, Z, u_y0, u_mu, beta_hat_cumsum[m-1, ], gamma_hat_cumsum[m-1, ], d, ds, p, ps, times, delta, should_print=should_print,
      iteration_number=m
    )
    gamma_hat[m, ] <- result$gamma_hat_addition
    gamma_hat_cumsum[m, ] <- result$gamma_hat_m
    beta_hat[m, ] <- result$beta_hat_addition
    beta_hat_cumsum[m, ] <- result$beta_hat_m
    if (boost_intercepts_continually) {
      gamma0 <- nlm(optimize_gamma_0, gamma_hat_cumsum[m, 1], beta_hat_cumsum[m, ], gamma_hat_cumsum[m, ], X, Z, times, delta)$estimate
      gamma_hat_cumsum[m, 1] <- gamma0
      beta0 <- nlm(optimize_beta_0, beta_hat_cumsum[m, 1], beta_hat_cumsum[m, ], gamma_hat_cumsum[m, ], X, Z, times, delta)$estimate
    } else {
      gamma0 <- gamma_hat_cumsum[m-1, 1]
      beta0 <- beta_hat_cumsum[m-1, 1]
    }
    gamma_hat_cumsum[m, 1] <- gamma0
    beta_hat_cumsum[m, 1] <- beta0
    negative_gradient_mu[m, ] <- FHT_componentwise_loss_function_derivative_mu(
      beta_hat_cumsum[m, ], gamma_hat_cumsum[m, ], X, Z, times, delta
    )
    negative_gradient_y0[m, ] <- FHT_componentwise_loss_function_derivative_y0(
      beta_hat_cumsum[m, ], gamma_hat_cumsum[m, ], X, Z, times, delta
    )

    if (should_print) {
      cat('gamma0: ', gamma_hat_cumsum[m, 1], '\n')
      cat('beta0: ', beta_hat_cumsum[m, 1], '\n')
    }

  }
  # Scale back
  gamma_hat_cumsum_scaled_back <- gamma_hat_cumsum
  beta_hat_cumsum_scaled_back <- beta_hat_cumsum
  # gamma_hat_cumsum_scaled_back <- destandardize(gamma_hat_cumsum, Z_scale_factors, Z_means)
  # beta_hat_cumsum_scaled_back <- destandardize(beta_hat_cumsum, X_scale_factors, X_means)
  parameters <- list(
    gamma_hats=gamma_hat_cumsum_scaled_back, beta_hats=beta_hat_cumsum_scaled_back
  )
  gamma_hat_final <- gamma_hat_cumsum_scaled_back[m, ]
  beta_hat_final <- beta_hat_cumsum_scaled_back[m, ]
  final_parameters <- list(gamma_hat_final=gamma_hat_final, beta_hat_final=beta_hat_final)
  for (m in 1:m_stop) {
    loss[m] <- FHT_minus_loglikelihood_with_all_parameters(
      beta=beta_hat_cumsum_scaled_back[m, ],
      gamma=gamma_hat_cumsum_scaled_back[m, ],
      X_original, Z_original, times, delta
    )
  }
  return(list(final_parameters=final_parameters, parameters=parameters, loss=loss))
}
