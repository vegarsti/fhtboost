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
  N <- dim(X)[1]

  # Run optimization
  initial_parameters <- runif(p+d, min=0.1, max=0.5)
  nlm_result <- nlm(minus_FHT_loglikelihood_nlm, initial_parameters)
  print(nlm_result$estimate)
  print(c(beta_true, gamma_true))
  print("Difference between true and nlm:")
  print(sum(abs(nlm_result$estimate - c(beta_true, gamma_true))))
  ### END SANITY CHECK ###


  m_stop <- 100

  nlm_parameter_list <- parameter_vector_to_list(nlm_result$estimate, d, p)
  beta_from_nlm <- nlm_parameter_list$beta
  gamma_from_nlm <- nlm_parameter_list$gamma

  # DIVIDE INTO K FOLDS
  K <- 5
  folds <- create_folds(N, K)
  M <- 100
  CV_error_matrix <- matrix(NA, nrow=m_stop, ncol=K)
  for (k in 1:K) {
    subset_without_k <- folds[-k, ]
    subset_k <- folds[k, ]
    X_without_k <- X[subset_without_k, ]
    Z_without_k <- Z[subset_without_k, ]
    X_k <- X[subset_k, ]
    Z_k <- Z[subset_k, ]
    times_without_k <- times[subset_without_k]
    delta_without_k <- delta[subset_without_k]
    times_k <- times[subset_k]
    delta_k <- times[subset_k]
    result <- boosting_run(
      times_without_k, delta_without_k, X_without_k, Z_without_k, M, nlm_parameter_list, boost_mu
    )
    if (boost_mu) {
      beta_ <- beta_from_nlm
      gamma_ <- result$parameters
      for (m in 1:m_stop) {
        gamma_m <- gamma_[m, ]
        CV_error_matrix[m, k] <- FHT_minus_loglikelihood_with_all_parameters(
          beta_, gamma_m, X_k, Z_k, times_k, delta_k
        )
      }
    } else {
      beta_ <- result$parameters
      gamma_ <- gamma_from_nlm
      for (m in 1:m_stop) {
        beta_m <- beta_[m, ]
        CV_error_matrix[m, k] <- FHT_minus_loglikelihood_with_all_parameters(
          beta_m, gamma_, X_k, Z_k, times_k, delta_k
        )
      }
    }
  }
  CV_errors <- rowSums(CV_error_matrix) / N
  plot(CV_errors, typ='l')

  lines(CV_error_matrix[, 1] / length(folds[1, ]), col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
  lines(CV_error_matrix[, 2] / length(folds[2, ]), col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
  lines(CV_error_matrix[, 3] / length(folds[3, ]), col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
  lines(CV_error_matrix[, 4] / length(folds[4, ]), col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))
  lines(CV_error_matrix[, 5] / length(folds[5, ]), col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))

  m_stop <- which.min(CV_errors)
  result <- boosting_run(times, delta, X, Z, m_stop, nlm_parameter_list, boost_mu)
}
