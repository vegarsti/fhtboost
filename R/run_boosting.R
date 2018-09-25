run_boosting <- function() {
  simulated_data <- simulate_FHT_data()
  times <- simulated_data$observations$survival_times
  delta <- simulated_data$observations$delta
  X <- simulated_data$design_matrices$X
  Z <- simulated_data$design_matrices$Z

  # center Z for correct boosting
  Z <- scale(Z, center=T, scale=F)

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
  #sum(abs(nlm_result$estimate - c(beta_true, gamma_true)))
  ### END SANITY CHECK ###

  beta_from_nlm <- nlm_result$estimate[1:2]
  gamma_from_nlm <- nlm_result$estimate[3:4]

  m_stop <- 10

  N <- dim(X)[1]

  # index_set <- 1:N
  # #index_set <- index_set[-c(146, 298)]
  #
  # X <- X[index_set, ]
  # Z <- Z[index_set, ]
  # times <- times[index_set]
  # delta <- delta[index_set]

  nu <- 0.1

  N <- dim(Z)[1]
  p <- dim(Z)[2]

  # Standardize X -- use scale instead!! (?)

  # X_original <- X
  # X_colmeans <- colMeans(X)
  # X_sds <- apply(X, 2, sd)
  # for (j in 1:p) {
  #   X[, j] <- (X[, j] - X_colmeans[j])/X_sds[j]
  # }

  # "make" loss function after scaling Z
  minus_FHT_loglikelihood <- data_to_FHT_loss_function(X, Z, times, delta)
  #minus_FHT_loglikelihood(list(beta=beta_true, gamma=gamma_true))

  gamma_hat <- matrix(NA, nrow=m_stop, ncol=p)
  gamma_hat[1, ] <- c(0, 0)
  negative_gradient <- matrix(NA, nrow=m_stop, ncol=N)
  parameter_list <- list(beta=beta_from_nlm, gamma=gamma_hat[1, ])
  #losses[1, ] <- FHT_componentwise_minus_loglikelihood_with_parameters(parameter_list, X, Z, times, delta)
  negative_gradient[1, ] <- - FHT_componentwise_loss_function_derivative_mu(parameter_list, X, Z, times, delta)
  loss <- rep(NA, m_stop)
  loss[1] <- minus_FHT_loglikelihood(parameter_list)

  u_hats_rss <- rep(NA, p)
  u_hats <- matrix(NA, nrow=p, ncol=N)
  gamma_hats <- rep(NA, p)
  for (m in 2:m_stop) {
    u <- negative_gradient[(m-1), ]
    for (j in 1:p) {
      Z_j <- Z[,j]
      gamma_hats[j] <- solve(t(Z_j) %*% Z_j) %*% t(Z_j) %*% u
      u_hats[j, ] <- Z_j * gamma_hats[j]
      u_hats_rss[j] <- sum((u_hats[j, ] - u)^2) / N
    }
    best_p <- which.min(u_hats_rss)
    gamma_hat[m, best_p] <- nu*gamma_hats[best_p]
    gamma_hat[m, -best_p] <- 0
    parameter_list <- list(beta=beta_from_nlm, gamma=colSums(gamma_hat[1:m, ]))
    negative_gradient[m, ] <- - FHT_componentwise_loss_function_derivative_mu(parameter_list, X, Z, times, delta)
    loss[m] <- minus_FHT_loglikelihood(parameter_list)
    m <- m + 1
  }

}
