run_boosting <- function() {
  boost_mu <- FALSE # if FALSE, boost y0 instead
  simulated_data <- simulate_FHT_data()
  times <- simulated_data$observations$survival_times
  delta <- simulated_data$observations$delta
  X <- simulated_data$design_matrices$X
  Z <- simulated_data$design_matrices$Z

  mu_has_intercept <- TRUE
  y0_has_intercept <- TRUE

  # center Z for correct boosting
  #Z <- scale(Z, center=T, scale=F)

  # Our loss function
  # usage: minus_FHT_loglikelihood(list(beta=beta, gamma=gamma))

  beta_true <- simulated_data$true_parameters$beta
  gamma_true <- simulated_data$true_parameters$gamma

  ### SANITY CHECK -> does nlm recover the parameters? ###
  minus_FHT_loglikelihood_nlm <- data_to_optimizable_function(X, Z, times, delta)

  d <- dim(X)[2]
  p <- dim(Z)[2]

  # Run optimization
  initial_parameters <- runif(p+d, min=0.1, max=0.5)
  nlm_result <- nlm(minus_FHT_loglikelihood_nlm, initial_parameters)
  #sum(abs(nlm_result$estimate - c(beta_true, gamma_true)))
  ### END SANITY CHECK ###


  m_stop <- 100

  N <- dim(X)[1]

  nu <- 0.1
  d <- dim(X)[2]
  p <- dim(Z)[2]
  N <- dim(Z)[1]
  # check dim(Z)[1] == dim(X)[1]

  nlm_parameter_list <- parameter_vector_to_list(nlm_result$estimate, d, p)
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
  #loss[1] <- minus_FHT_loglikelihood(parameter_list)

  gamma_hats <- rep(NA, p)
  u_hats_rss <- rep(NA, p)
  u_hats <- matrix(NA, nrow=p, ncol=N)
  beta_hats <- rep(NA, d)


  # DIVIDE INTO K FOLDS
  K <- 5
  folds <- create_folds(N, K)

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
  CV_means <- rep(0, m_stop)
  for (m in 2:m_stop) {
    u <- negative_gradient[(m-1), ]
    ## BOOST MU
    if (boost_mu) {
      gamma_hat[m, ] <- nu*best_least_squares_update(Z, u, p, ps)
      gamma_hat_cumsum[m, ] <- gamma_hat_cumsum[m-1, ] + gamma_hat[m, ]
      negative_gradient[m, ] <- FHT_componentwise_loss_function_derivative_mu(
        beta_from_nlm,
        gamma_hat_cumsum[m, ],
        X,
        Z,
        times,
        delta
      )
      parameter_list <- list(beta=beta_from_nlm, gamma=gamma_hat_cumsum[m, ])

      # # CV
      # for (k in 1:K) {
      #   X_k <- X[-folds[k], ]
      #   Z_k <- Z[-folds[k, ]]
      #   times_k <- times(-folds[k])
      #   delta_k <- delta(-folds[k])
      #   FHT_componentwise_minus_loglikelihood_with_parameters(
      #     beta_,
      #     gamma_,
      #     X,
      #     Z,
      #     times,
      #     delta
      #   )
      # }

    } else {
      beta_hat[m, ] <- nu*best_least_squares_update(X, u, d, ds)
      beta_hat_cumsum[m, ] <- beta_hat_cumsum[m-1, ] + beta_hat[m, ]
      negative_gradient[m, ] <- FHT_componentwise_loss_function_derivative_y0(
        beta_hat_cumsum[m, ],
        gamma_from_nlm,
        X,
        Z,
        times,
        delta
      )
      parameter_list <- list(beta=beta_hat_cumsum[m, ], gamma=gamma_from_nlm)
    }
    loss[m] <- minus_FHT_loglikelihood(parameter_list)
  }
  #gamma_hat_final <- gamma_hat_cumsum[m, ] * Z_scale_factors
  beta_hat_final <- beta_hat_cumsum[m, ] * X_scale_factors
  gradient_sums <- rowSums(abs(negative_gradient))
  plot(gradient_sums, typ='l')
}
