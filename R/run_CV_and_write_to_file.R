run_CV_and_write_to_file <- function(N, setup_type, add_noise, seed, K, K_fold_repetitions, directory, M, criterion) {
  # simulate data (faster than reading...)
  simulated_data <- simulate_FHT_data(N=N, setup_type=setup_type, add_noise=add_noise, seed=seed)
  times <- simulated_data$observations$survival_times
  delta <- simulated_data$observations$delta
  X <- simulated_data$design_matrices$X
  Z <- simulated_data$design_matrices$Z
  CV_result <- run_CV(M, K_fold_repetitions, K, X, Z, times, delta)
  CV_errors_K_deviance <- CV_result$CV_errors_K_deviance
  CV_errors_K_loglik <- CV_result$CV_errors_K_loglik

  # Write to file
  descriptor <- paste("cv", 'loglik', sep="_")
  full_filename <- make_filename(directory, descriptor, seed)
  write.csv(CV_errors_K_loglik, file=full_filename, row.names=FALSE)
  # Write to file
  descriptor <- paste("cv", 'deviance', sep="_")
  full_filename <- make_filename(directory, descriptor, seed)
  write.csv(CV_errors_K_deviance, file=full_filename, row.names=FALSE)
}
