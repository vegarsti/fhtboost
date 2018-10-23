run_boosting <- function() {
  # Extract simulated data
  simulated_data <- simulate_FHT_data(dense=FALSE)
  times <- simulated_data$observations$survival_times
  delta <- simulated_data$observations$delta
  non_para <- non_parametric_estimates(times, delta, continuous = TRUE)
  plot(non_para$times_sequence, non_para$kaplan_meiers, typ='s')
  X <- simulated_data$design_matrices$X
  Z <- simulated_data$design_matrices$Z

  # dimensions
  d <- dim(X)[2]
  p <- dim(Z)[2]
  N <- dim(X)[1]

  beta_true <- simulated_data$true_parameters$beta
  gamma_true <- simulated_data$true_parameters$gamma

  ### SANITY CHECK -> does nlm recover the parameters? ###
  minus_FHT_loglikelihood_nlm <- data_to_optimizable_function(X, Z, times, delta)

  # Run optimization
  initial_parameters <- runif(p+d, min=0.1, max=0.5)
  nlm_result <- nlm(minus_FHT_loglikelihood_nlm, initial_parameters)
  print(nlm_result$estimate)

  # y0 <- exp(nlm_result$estimate[1])
  # mu <- nlm_result$estimate[d+1]
  # parametric_times <- seq(0.8, max(times), by=0.01)
  # parametric_S <- FHT_parametric_survival(parametric_times, mu, y0)
  # plot(non_para$times_sequence, non_para$kaplan_meiers, typ='s', ylim=c(0, 1))
  # lines(parametric_times, parametric_S, col='red')
  #
  # parametric_A <- FHT_parametric_cumulative_hazard(parametric_times, mu, y0)
  # plot(non_para$times_sequence, non_para$nelson_aalens, typ='s')
  # lines(parametric_times, parametric_A, col='red')
  # abline(-0.7, 0.5)

  print(c(beta_true, gamma_true))
  print("Difference between true and nlm:")
  print(sum(abs(nlm_result$estimate - c(beta_true, gamma_true))))

  nlm_parameter_list <- parameter_vector_to_list(nlm_result$estimate, d, p)
  beta_from_nlm <- nlm_parameter_list$beta
  beta_0_from_nlm <- beta_from_nlm[1]
  gamma_from_nlm <- nlm_parameter_list$gamma
  gamma_0_from_nlm <- gamma_from_nlm[1]

  # DIVIDE INTO K FOLDS
  K <- 5
  K_fold_repetitions <- 10
  M <- m_stop <- 50 ### M STOP
  CV_errors_mu_K <- matrix(NA, nrow=m_stop, ncol=K_fold_repetitions)
  CV_errors_y0_K <- matrix(NA, nrow=m_stop, ncol=K_fold_repetitions)
  loss_K <- matrix(NA, nrow=m_stop, ncol=K_fold_repetitions)
  for (repeated_K_fold_iteration in 1:K_fold_repetitions) {
    folds <- create_folds(N, K)
    CV_error_matrix_mu <- matrix(NA, nrow=m_stop, ncol=K)
    CV_error_matrix_y0 <- matrix(NA, nrow=m_stop, ncol=K)
    loss <- matrix(NA, nrow=M, ncol=K)
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
        times_without_k, delta_without_k, X_without_k, Z_without_k, M, M, beta_0_from_nlm, gamma_0_from_nlm, use_nlm=TRUE
      )
      beta_hats <- result$parameters$beta_hats
      gamma_hats <- result$parameters$gamma_hats
      for (m in 2:M) {
        gamma_m <- gamma_hats[m, ]
        gamma_m1 <- gamma_hats[m-1, ]
        beta_m <- beta_hats[m, ]
        beta_m1 <- beta_hats[m-1, ]
        CV_error_matrix_mu[m, k] <- FHT_minus_loglikelihood_with_all_parameters(
          beta_m1, gamma_m, X_k, Z_k, times_k, delta_k
        )
        CV_error_matrix_y0[m, k] <- FHT_minus_loglikelihood_with_all_parameters(
          beta_m, gamma_m1, X_k, Z_k, times_k, delta_k
        )
      }
      # ???
      CV_error_matrix_mu[1, k] <- CV_error_matrix_mu[2, k]
      CV_error_matrix_y0[1, k] <- CV_error_matrix_y0[2, k]
      loss[, k] <- result$loss / length(times_k)
    }
    CV_errors_mu_K[, repeated_K_fold_iteration] <- rowSums(CV_error_matrix_mu) / N
    CV_errors_y0_K[, repeated_K_fold_iteration] <- rowSums(CV_error_matrix_y0) / N
    loss_K[, repeated_K_fold_iteration] <- rowMeans(loss)
  }
  CV_errors_mu <- rowSums(CV_errors_mu_K)
  CV_errors_y0 <- rowSums(CV_errors_y0_K)
  loss <- rowMeans(loss_K)
  plot(CV_errors_mu, typ='l')
  lines(CV_errors_y0, typ='l', col='red')

  # lines(CV_error_matrix[, 1] / length(folds[1, ]),
  #   col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))

  m_stop_mu <- which.min(CV_errors_mu)
  m_stop_y0 <- which.min(CV_errors_y0)
  result_wo_nlm <- boosting_run(times, delta, X, Z, m_stop_mu, m_stop_y0, beta_0_from_nlm, gamma_0_from_nlm, use_nlm=FALSE)
  result_w_nlm <- boosting_run(times, delta, X, Z, m_stop_mu, m_stop_y0, beta_0_from_nlm, gamma_0_from_nlm, use_nlm=TRUE)
}
