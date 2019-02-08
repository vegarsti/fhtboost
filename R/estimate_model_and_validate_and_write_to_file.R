estimate_model_and_validate_and_write_to_file <- function(N, setup_type, add_noise, seed, directory) {
  criterion <- 'deviance'
  full_filename <- make_filename(directory, "cv", criterion, seed)
  CV_errors_K <- readr::read_csv(full_filename)
  CV_errors <- rowMeans(CV_errors_K)
  if (criterion == 'deviance') {
    m_stop <- which.max(CV_errors)
  } else if (criterion == 'loglik') {
    m_stop <- which.min(CV_errors)
  }
  simulated_data <- simulate_FHT_data(N=N, setup_type=setup_type, add_noise=add_noise, seed=seed)
  times <- simulated_data$observations$survival_times
  delta <- simulated_data$observations$delta
  X <- simulated_data$design_matrices$X
  Z <- simulated_data$design_matrices$Z
  best_intercepts <- maximum_likelihood_intercepts(times, delta)
  y0 <- best_intercepts[1]
  mu <- best_intercepts[2]
  null_y0 <- y0
  null_mu <- mu
  null_model_loglikelihood <- - sum(FHT_loglikelihood_with_y0_mu(y0, mu, times, delta))
  result <- boosting_run(times, delta, X, Z, m_stop, boost_intercepts_continually=TRUE, should_print=FALSE, run_in_parallel=FALSE)
  simulated_data_test <- simulate_FHT_data(N=N_test, setup_type=setup_type, add_noise=add_noise, seed=TEST_SEED)
  times_test <- simulated_data_test$observations$survival_times
  delta_test <- simulated_data_test$observations$delta
  X_test <- simulated_data_test$design_matrices$X
  Z_test <- simulated_data_test$design_matrices$Z
  estimated_y0s_test <- exp(X_test %*% result$final_parameters$beta_hat_final)
  estimated_mus_test <- Z_test %*% result$final_parameters$gamma_hat_final
  best_intercepts_test <- maximum_likelihood_intercepts(times_test, delta_test)
  null_y0_test <- best_intercepts[1]
  null_mu_test <- best_intercepts[2]

  # write parameters
  final_beta <- result$final_parameters$beta_hat_final
  full_filename <- make_filename(directory, "beta", seed)
  write.csv(final_beta, full_filename)
  final_gamma <- result$final_parameters$gamma_hat_final
  full_filename <- make_filename(directory, "gamma", seed)
  write.csv(final_gamma, full_filename)

  # write y_0 and mu
  y0 <- exp(X %*% final_beta)
  full_filename <- make_filename(directory, "y0", seed)
  write.csv(y0, full_filename)
  mu <- Z %*% final_gamma
  full_filename <- make_filename(directory, "mu", seed)
  write.csv(mu, full_filename)

  # calculations
  beta_null <- c(null_y0, rep(0, (dim(X_test)[2]-1)))
  gamma_null <- c(null_mu, rep(0, (dim(Z_test)[2]-1)))
  null_model_loglikelihood <- - FHT_minus_loglikelihood_with_all_parameters(beta_null, gamma_null, X_test, Z_test, times_test, delta_test)
  full_model_loglikelihood <- - FHT_minus_loglikelihood_with_all_parameters(result$final_parameters$beta_hat_final, result$final_parameters$gamma_hat_final, X_test, Z_test, times_test, delta_test)
  deviance <- 2 * (full_model_loglikelihood - null_model_loglikelihood)
  # summary (deviance and m_stop)
  full_filename <- make_filename(directory, "summary", seed)
  write.csv(data.frame('deviance'=deviance, 'm_stop'=m_stop), file=full_filename, row.names=FALSE)
  full_filename <- make_filename(directory, "m", seed)
  write.csv(as.data.frame(m_stop), file=full_filename, row.names=FALSE)

  # Brier R2
  #brier_r2_result <- brier_r2(times_test, delta_test, estimated_y0s_test, estimated_mus_test, null_y0_test, null_mu_test, number_of_time_points=100)
  #full_filename <- make_filename(directory, "brier", seed)
  #write.csv(brier_r2_result, file=full_filename, row.names=FALSE)
}
