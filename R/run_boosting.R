run_boosting <- function() {
  # Extract simulated data
  dense <- TRUE
  simulated_data <- simulate_FHT_data(dense=dense)
  times <- simulated_data$observations$survival_times
  delta <- simulated_data$observations$delta
  X <- simulated_data$design_matrices$X
  Z <- simulated_data$design_matrices$Z
  beta_true <- simulated_data$true_parameters$beta
  gamma_true <- simulated_data$true_parameters$gamma

  # Lymph high-dim data
  # delta <- as.numeric(read.csv('../lymphstatus.txt', header=FALSE)[, 1])
  # times <- as.numeric(read.csv('../lymphtim.txt', header=FALSE)[, 1])
  # X <- as.numeric(read.csv('../lymphx.txt', sep=' ', header=FALSE))[1:10, ]
  # Z <- X

  # non_para <- non_parametric_estimates(times, delta, continuous = TRUE)
  # plot(non_para$times_sequence, non_para$kaplan_meiers, typ='s')

  # dimensions
  d <- dim(X)[2]
  p <- dim(Z)[2]
  N <- dim(X)[1]

  ### SANITY CHECK -> does nlm recover the parameters? ###
  minus_FHT_loglikelihood_nlm <- data_to_optimizable_function(X, Z, times, delta)
  #
  # # Run optimization
  initial_parameters <- runif(p+d, min=0.1, max=0.5)
  nlm_result <- nlm(minus_FHT_loglikelihood_nlm, initial_parameters)
  # beta_nlm <- nlm_result$estimate[1:d]
  # gamma_nlm <- nlm_result$estimate[(d+1):(d+p)]
  # nlm_loss <- FHT_minus_loglikelihood_with_all_parameters(beta_nlm, gamma_nlm, X, Z, times, delta)
  # y0 <- exp(nlm_result$estimate[1])
  # mu <- nlm_result$estimate[d+1]
  # parametric_times <- seq(0.8, max(times), by=0.01)
  # parametric_S <- FHT_parametric_survival(parametric_times, mu, y0)
  # plot(non_para$times_sequence, non_para$kaplan_meiers, typ='s', ylim=c(0, 1))
  # lines(parametric_times, parametric_S, col='red')
  # plot(parametric_times, parametric_S, col='red')
  #
  # parametric_A <- FHT_parametric_cumulative_hazard(parametric_times, mu, y0)
  # plot(non_para$times_sequence, non_para$nelson_aalens, typ='s')
  # lines(parametric_times, parametric_A, col='red')
  # abline(-0.7, 0.5)

  nlm_parameter_list <- parameter_vector_to_list(nlm_result$estimate, d, p)
  beta_from_nlm <- nlm_parameter_list$beta
  beta_0_from_nlm <- beta_from_nlm[1]
  gamma_from_nlm <- nlm_parameter_list$gamma
  gamma_0_from_nlm <- gamma_from_nlm[1]

  result <- boosting_run(
    times, delta, X, Z, m_stop_mu=100, m_stop_y0=100, beta_0_from_nlm, gamma_0_from_nlm,
    give_intercepts=FALSE, optimize_intercepts=TRUE
  )

  # DIVIDE INTO K FOLDS
  K <- 10
  K_fold_repetitions <- 10
  if (dense) {
    M <- m_stop <- 30 ### M STOP
  } else {
    M <- m_stop <- 10 ### M STOP
  }
  CV_errors_K <- matrix(NA, nrow=m_stop, ncol=K_fold_repetitions)
  CV_errors_K <- matrix(NA, nrow=m_stop, ncol=K_fold_repetitions)
  loss_K <- matrix(NA, nrow=m_stop, ncol=K_fold_repetitions)
  for (repeated_K_fold_iteration in 1:K_fold_repetitions) {
    folds <- create_folds_stratified(delta, K)
    CV_error_matrix <- matrix(NA, nrow=m_stop, ncol=K)
    loss <- matrix(NA, nrow=M, ncol=K)
    for (k in 1:K) {
      subset_without_k <- get_all_but_kth_fold(folds, k, K)
      subset_k <- get_kth_fold(folds, k)
      X_without_k <- X[subset_without_k, ]
      Z_without_k <- Z[subset_without_k, ]
      X_k <- X[subset_k, ]
      Z_k <- Z[subset_k, ]
      times_without_k <- times[subset_without_k]
      delta_without_k <- delta[subset_without_k]
      times_k <- times[subset_k]
      delta_k <- times[subset_k]
      result <- boosting_run(
        times_without_k, delta_without_k, X_without_k, Z_without_k, M, M, beta_0_from_nlm, gamma_0_from_nlm,
        give_intercepts=FALSE, optimize_intercepts=TRUE
      )
      beta_hats <- result$parameters$beta_hats
      gamma_hats <- result$parameters$gamma_hats
      for (m in 2:M) {
        gamma_m <- gamma_hats[m, ]
        gamma_m1 <- gamma_hats[m-1, ]
        beta_m <- beta_hats[m, ]
        beta_m1 <- beta_hats[m-1, ]
        CV_error_matrix[m, k] <- FHT_minus_loglikelihood_with_all_parameters(
          beta_m1, gamma_m, X_k, Z_k, times_k, delta_k
        )
      }
      # ???
      CV_error_matrix[1, k] <- CV_error_matrix[2, k]
      loss[, k] <- result$loss
    }
    CV_errors_K[, repeated_K_fold_iteration] <- rowSums(CV_error_matrix)
    loss_K[, repeated_K_fold_iteration] <- rowSums(loss)
  }
  CV_errors <- rowSums(CV_errors_K)
  loss <- rowSums(loss_K)
  plot(CV_errors, typ='l')
  # lines(CV_error_matrix[, 1] / length(folds[1, ]),
  #   col=rgb(red = 0, green = 0, blue = 0, alpha = 0.5))

  m_stop <- which.min(CV_errors)
  m_stop_mu <- m_stop
  m_stop_y0 <- m_stop

  result_w_nlm <- boosting_run(times, delta, X, Z, m_stop_mu+20, m_stop_y0+20, beta_0_from_nlm, gamma_0_from_nlm, give_intercepts=TRUE, optimize_intercepts=FALSE)
  result_w_nlm_w_boost <- boosting_run(times, delta, X, Z, m_stop_mu+20, m_stop_y0+20, beta_0_from_nlm, gamma_0_from_nlm, give_intercepts=TRUE, optimize_intercepts=TRUE)
  result_wo_nlm <- boosting_run(times, delta, X, Z, m_stop+20, m_stop+20, beta_0_from_nlm, gamma_0_from_nlm, give_intercepts=FALSE, optimize_intercepts=TRUE)
  result_wo_nlm_2 <- boosting_run(times, delta, X, Z, m_stop+20, m_stop+20, beta_0_from_nlm, gamma_0_from_nlm, give_intercepts=FALSE, optimize_intercepts=TRUE)
  result2 <- boosting_run(times, delta, X, Z, m_stop_mu+20, m_stop_y0+20, beta_0_from_nlm, gamma_0_from_nlm, give_intercepts=FALSE, optimize_intercepts=FALSE)

  # Plot loss functions
  plot(result_wo_nlm$loss, typ='l', lty=2, main='Loss with/without changing intercept while boosting', ylim=c(1200, 1500))
  lines(result_wo_nlm_2$loss, col='blue')
  lines(result2$loss, col='red')
  lines(result_w_nlm$loss, col='black', lty=3)
  abline(h=nlm_loss, col='blue')
  abline(v=m_stop_mu, col='blue', lty=3)
  legend("topright", legend=c("With", "Without", "ML intercept, no changing", "ML (nlm)", "m_stop"), col=c("black", "red", "black", "blue", "blue"), lty=c(2, 1, 3, 1, 3), cex=1.5)

  # Plot beta0
  joint_optimized <- nlm_result$estimate[1]
  beta0s <- result_wo_nlm$parameters$beta_hats[, 1]
  plot(beta0s, typ='l', ylim=c(0.5, 2), lty=2, main='Beta intercept (in y0)')
  abline(h=joint_optimized, col='red')
  abline(h=beta0s[1], col='blue')
  abline(v=m_stop_mu, col='blue', lty=3)
  legend("topleft", legend=c("Path", "Initial", "ML", "m_stop"), col=c("black", "red", "blue", "blue"), lty=c(2, 1, 1, 3), cex=1.5)

  # Plot gamma0
  joint_optimized <- nlm_result$estimate[d+1]
  gamma0s <- result_wo_nlm$parameters$gamma_hats[, 1]
  plot(gamma0s, typ='l', ylim=c(-1, 0), lty=2, main='Gamma intercept (in mu)')
  abline(h=joint_optimized, col='red')
  abline(h=gamma0s[1], col='blue')
  abline(v=m_stop_mu, col='blue', lty=3)
  legend("topleft", legend=c("Path", "Initial", "ML", "m_stop"), col=c("black", "red", "blue", "blue"), lty=c(2, 1, 1, 3), cex=1.5)
}
