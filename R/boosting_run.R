boosting_run <- function(times, delta, X, Z, m_stop_mu, m_stop_y0, beta_0_from_nlm, gamma_0_from_nlm, give_intercepts=TRUE) {
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

  m_stop <- max(m_stop_mu, m_stop_y0)

  # Nuisance parameter: shrinkage factor
  nu <- 0.1

  # Scaling X and Z -- should rewrite! or hide in function
  X_original <- X
  Z_original <- Z
  X_scale_factors <- as.numeric(apply(X, 2, sd))
  Z_scale_factors <- as.numeric(apply(Z, 2, sd))
  X_means <- as.numeric(apply, X, 2, mean)
  Z_means <- as.numeric(apply, Z, 2, mean)
  X_scale_factors[1] <- 1
  Z_scale_factors[1] <- 1
  X_means[1] <- 0
  Z_means[1] <- 0
  for (j in 2:d) {
    X[, j] <- X[, j]/X_scale_factors[j] - X_means[j]
  }
  for (j in 2:p) {
    Z[, j] <- Z[, j]/Z_scale_factors[j] - Z_means[j]
  }

  # "make" loss function after scaling Z
  minus_FHT_loglikelihood <- data_to_FHT_loss_function(X, Z, times, delta)
  #minus_FHT_loglikelihood(list(beta=beta_true, gamma=gamma_true))

  # initialize with nlm
  function_to_optimize <- function(parameters) {
    y0_baseline <- exp(parameters[1])
    mu_baseline <- parameters[2]
    return(-sum(FHT_loglikelihood_with_y0_mu(y0_baseline, mu_baseline, times, delta)))
  }

  optimize_beta_0 <- function(beta0, beta_, gamma_, X, Z, times, delta) {
    beta_[1] <- beta0
    return(FHT_minus_loglikelihood_with_all_parameters(beta_, gamma_, X, Z, times, delta))
  }

  optimize_gamma_0 <- function(gamma0, beta_, gamma_, X, Z, times, delta) {
    gamma_[1] <- gamma0
    return(FHT_minus_loglikelihood_with_all_parameters(beta_, gamma_, X, Z, times, delta))
  }

  nlm_result <- nlm(function_to_optimize, c(0.1, 0.1))

  gamma_hat <- matrix(NA, nrow=m_stop, ncol=p)
  gamma_hat[1, ] <- rep(0, p)
  gamma_hat_cumsum <- matrix(NA, nrow=m_stop, ncol=p)

  beta_hat <- matrix(NA, nrow=m_stop, ncol=d)
  beta_hat[1, ] <- rep(0, d)
  beta_hat_cumsum <- matrix(NA, nrow=m_stop, ncol=d)

  if (give_intercepts) {
    gamma_0 <- gamma_0_from_nlm
    beta_0 <- beta_0_from_nlm
  } else {
    beta_0 <- nlm_result$estimate[1]
    gamma_0 <- nlm_result$estimate[2]
  }

  # INITIALIZE
  gamma_hat[1, 1] <- gamma_0
  beta_hat[1, 1] <- beta_0
  gamma_hat_cumsum[1, ] <- gamma_hat[1, ]
  beta_hat_cumsum[1, ] <- beta_hat[1, ]

  negative_gradient_mu <- matrix(NA, nrow=m_stop, ncol=N)
  negative_gradient_y0 <- matrix(NA, nrow=m_stop, ncol=N)
  negative_gradient_mu[1, ] <- FHT_componentwise_loss_function_derivative_mu(beta_hat[1, ], gamma_hat[1, ], X, Z, times, delta)
  negative_gradient_y0[1, ] <- FHT_componentwise_loss_function_derivative_mu(beta_hat[1, ], gamma_hat[1, ], X, Z, times, delta)
  loss <- rep(NA, m_stop)
  parameter_list <- list(beta=beta_hat[1, ], gamma_hat[1, ])
  #loss[1] <- minus_FHT_loglikelihood(parameter_list)

  gamma_hats <- rep(NA, p)
  beta_hats <- rep(NA, d)
  u_hats <- matrix(NA, nrow=p, ncol=N)
  u_hats_rss <- rep(NA, p)
  loss <- rep(0, m_stop)

  for (m in 2:m_stop) {
    #print(m)
    # gamma/mu
    if (m <= m_stop_mu) {
      #print("mu")
      u <- negative_gradient_mu[(m-1), ]
      result <- boosting_iteration_mu(
        nu, X, Z, u, beta_hat_cumsum[m-1, ], gamma_hat_cumsum[m-1, ], p, ps, times, delta
      )
      gamma_hat[m, ] <- result$gamma_hat_addition
      gamma_hat_cumsum[m, ] <- result$gamma_hat_m
      negative_gradient_mu[m, ] <- result$u_m
    } else {
      gamma_hat[m, ] <- 0
      gamma_hat_cumsum[m, ] <- gamma_hat_cumsum[m-1, ]
    }
    gamma0 <- nlm(optimize_gamma_0, gamma_hat_cumsum[m-1, 1], beta_hat_cumsum[m-1, ], gamma_hat_cumsum[m, ], X, Z, times, delta)$estimate
    gamma_hat_cumsum[m, 1] <- gamma0
    # beta/y0
    if (m <= m_stop_y0) {
      u <- negative_gradient_y0[(m-1), ]
      result <- boosting_iteration_y0(
        nu, X, Z, u, beta_hat_cumsum[m-1, ], gamma_hat_cumsum[m, ], d, ds, times, delta
      )
      beta_hat[m, ] <- result$beta_hat_addition
      beta_hat_cumsum[m, ] <- result$beta_hat_m
      negative_gradient_y0[m, ] <- result$u_m
    } else {
      beta_hat[m, ] <- 0
      beta_hat_cumsum[m, ] <- beta_hat_cumsum[m-1, ]
    }
    beta0 <- nlm(optimize_beta_0, beta_hat_cumsum[m-1, 1], beta_hat_cumsum[m, ], gamma_hat_cumsum[m, ], X, Z, times, delta)$estimate
    beta_hat_cumsum[m, 1] <- beta0
  }
  # Scale back
  gamma_hat_cumsum_scaled_back <- destandardize(gamma_hat_cumsum, Z_scale_factors, Z_means)
  beta_hat_cumsum_scaled_back <- destandardize(beta_hat_cumsum, X_scale_factors, X_means)

  gamma_hat_final <- gamma_hat_cumsum[m, ] * Z_scale_factors
  beta_hat_final <- beta_hat_cumsum[m, ] * X_scale_factors
  parameters <- list(
    gamma_hats=gamma_hat_cumsum_scaled_back, beta_hats=beta_hat_cumsum_scaled_back
  )

  for (m in 1:m_stop) {
    loss[m] <- FHT_minus_loglikelihood_with_all_parameters(
      beta=beta_hat_cumsum_scaled_back[m, ],
      gamma=gamma_hat_cumsum_scaled_back[m, ],
      X_original, Z_original, times, delta
    )
  }

  final_parameters <- list(gamma_hat_final=gamma_hat_final, beta_hat_final=beta_hat_final)
  return(list(final_parameters=final_parameters, parameters=parameters, loss=loss))
}
