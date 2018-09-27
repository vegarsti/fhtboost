run_boosting <- function() {
  boost_mu <- FALSE # boost y0 instead
  simulated_data <- simulate_FHT_data()
  times <- simulated_data$observations$survival_times
  delta <- simulated_data$observations$delta
  X <- simulated_data$design_matrices$X
  Z <- simulated_data$design_matrices$Z

  # center Z for correct boosting
  #Z <- scale(Z, center=T, scale=F)

  # Our loss function
  # usage: minus_FHT_loglikelihood(list(beta=beta, gamma=gamma))

  beta_true <- c(0.5, 1)       # the true beta!
  gamma_true <- c(-0.2, -0.2)

  ### SANITY CHECK -> does nlm recover the parameters? ###
  minus_FHT_loglikelihood_nlm <- data_to_optimizable_function(X, Z, times, delta)

  p <- dim(X)[2]
  d <- dim(Z)[2]

  # Run optimization
  initial_parameters <- runif(p+d, min=0.1, max=0.5)
  nlm_result <- nlm(minus_FHT_loglikelihood_nlm, initial_parameters)
  sum(abs(nlm_result$estimate - c(beta_true, gamma_true)))
  ### END SANITY CHECK ###

  beta_from_nlm <- nlm_result$estimate[1:2]
  gamma_from_nlm <- nlm_result$estimate[3:4]

  m_stop <- 100

  N <- dim(X)[1]

  nu <- 0.1

  N <- dim(Z)[1]
  p <- dim(Z)[2]

  # Standardize X & Z -- use scale instead!

  # "make" loss function after scaling Z
  minus_FHT_loglikelihood <- data_to_FHT_loss_function(X, Z, times, delta)
  #minus_FHT_loglikelihood(list(beta=beta_true, gamma=gamma_true))

  gamma_hat <- matrix(NA, nrow=m_stop, ncol=p)
  gamma_hat[1, ] <- c(0, 0)
  gamma_hat_cumsum <- matrix(NA, nrow=m_stop, ncol=p)
  gamma_hat_cumsum[1, ] <- c(0, 0)

  beta_hat <- matrix(NA, nrow=m_stop, ncol=p)
  beta_hat[1, ] <- c(0, 0)
  beta_hat_cumsum <- matrix(NA, nrow=m_stop, ncol=p)
  beta_hat_cumsum[1, ] <- c(0, 0)


  negative_gradient <- matrix(NA, nrow=m_stop, ncol=N)
  #parameter_list <- list(beta=beta_from_nlm, gamma=gamma_hat[1, ])
  parameter_list <- list(beta=beta_hat[1, ], gamma=gamma_from_nlm)
  negative_gradient[1, ] <- FHT_componentwise_loss_function_derivative_mu(parameter_list, X, Z, times, delta)
  loss <- rep(NA, m_stop)
  #loss[1] <- minus_FHT_loglikelihood(parameter_list)

  u_hats_rss <- rep(NA, p)
  u_hats <- matrix(NA, nrow=p, ncol=N)
  gamma_hats <- rep(NA, p)
  beta_hats <- rep(NA, p)
  for (m in 2:m_stop) {
    u <- negative_gradient[(m-1), ]
    for (j in 1:p) {

      ## BOOST MU
      if (boost_mu == TRUE) {
        Z_j <- Z[,j]
        gamma_hats[j] <- solve(t(Z_j) %*% Z_j) %*% t(Z_j) %*% u
        u_hats[j, ] <- Z_j * gamma_hats[j]
        u_hats_rss[j] <- sum((u_hats[j, ] - u)^2) / N
      }

      ## BOOST y0
      else {
        X_j <- X[,j]
        beta_hats[j] <- solve(t(X_j) %*% X_j) %*% t(X_j) %*% u
        u_hats[j, ] <- X_j * beta_hats[j]
        u_hats_rss[j] <- sum((u_hats[j, ] - u)^2) / N
      }
    }
    best_p <- which.min(u_hats_rss)

    if (boost_mu == TRUE) {
      # BOOST MU
      gamma_hat[m, best_p] <- nu*gamma_hats[best_p]
      gamma_hat[m, -best_p] <- 0
      gamma_hat_cumsum[m, ] <- gamma_hat_cumsum[m-1, ] + gamma_hat[m, ]
      parameter_list <- list(beta=beta_from_nlm, gamma=gamma_hat_cumsum[m, ])
      negative_gradient[m, ] <- FHT_componentwise_loss_function_derivative_mu(
        beta_from_nlm,
        gamma_hat_cumsum[m, ],
        Z,
        times,
        delta
      )
    }
    # BOOST y0
    else {
      beta_hat[m, best_p] <- nu*beta_hats[best_p]
      beta_hat[m, -best_p] <- 0
      beta_hat_cumsum[m, ] <- beta_hat_cumsum[m-1, ] + beta_hat[m, ]
      negative_gradient[m, ] <- FHT_componentwise_loss_function_derivative_y0(
        beta_hat_cumsum[m, ],
        gamma_from_nlm,
        Z,
        times,
        delta
      )
    }
    #loss[m] <- minus_FHT_loglikelihood(parameter_list)
  }
  gradient_sums <- rowSums(abs(negative_gradient))
  plot(gradient_sums, typ='l')
}
