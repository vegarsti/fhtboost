rm(list=ls())

library(devtools)
library(profvis) # pause
library(foreach)
library(readr)
library(doParallel)
load_all()

# Configuration
do_CV <- FALSE
train_models <- TRUE
run_in_parallel <- TRUE
number_of_data_sets <- 100

#seeds <- 1:number_of_data_sets
#seeds <- c(1)

no_cores <- detectCores() - 1
registerDoParallel(cores=no_cores)
cl <- makeCluster(no_cores)
directory <- "../dataset/"
N <- 500
N_test <- 1000
setup_type <- 'huge_clinical'
add_noise <- FALSE
K <- 5
K_fold_repetitions <- 5
M <- 80
TEST_SEED <- 9000
criterion <- 'deviance'

make_filename <- function(directory, descriptor, seed) {
  seed_string <- formatC(seed, width=4, flag="0")
  filename <- paste(seed_string, descriptor, sep="_")
  full_filename <- paste(directory, filename,  '.csv', sep='')
  return(full_filename)
}

run_CV_and_write_to_file <- function(N, setup_type, add_noise, seed, K, K_fold_repetitions, directory, M, criterion) {
  # simulate data (faster than reading...)
  simulated_data <- simulate_FHT_data(N=N, setup_type=setup_type, add_noise=add_noise, seed=seed)
  times <- simulated_data$observations$survival_times
  delta <- simulated_data$observations$delta
  X <- simulated_data$design_matrices$X
  Z <- simulated_data$design_matrices$Z
  CV_result <- run_CV(M, K_fold_repetitions, K, X, Z, times, delta, criterion=criterion)
  CV_errors_K <- CV_result$CV_errors_K

  if (should_plot) {
    plot(rowMeans(CV_errors_K), typ='l', lty=1)
    Ks <- dim(CV_errors_K)[2]
    for (k in 1:Ks) {
      lines(CV_errors_K[, k], lty=3)
    }
  }
  # Write to file
  descriptor <- paste("cv", criterion, sep="_")
  full_filename <- make_filename(directory, descriptor, seed)
  write.csv(CV_errors_K, file=full_filename, row.names=FALSE)
}

estimate_model_and_validate_and_write_to_file <- function(N, setup_type, add_noise, seed, K, K_fold_repetitions, directory, M) {
  full_filename <- make_filename(directory, "cv", seed)
  CV_errors_K <- readr::read_csv(full_filename)
  CV_errors <- rowMeans(CV_errors_K)
  m_stop <- which.min(CV_errors)
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

  # write non-null parameters
  final_beta <- result$final_parameters$beta_hat_final
  full_filename <- make_filename(directory, "beta", seed)
  write.csv(final_beta, full_filename)
  final_gamma <- result$final_parameters$gamma_hat_final
  full_filename <- make_filename(directory, "gamma", seed)
  write.csv(final_gamma, full_filename)

  # calculations
  beta_null <- c(null_y0, rep(0, (dim(X_test)[2]-1)))
  gamma_null <- c(null_mu, rep(0, (dim(Z_test)[2]-1)))
  null_model_loglikelihood <- - FHT_minus_loglikelihood_with_all_parameters(beta_null, gamma_null, X_test, Z_test, times_test, delta_test)
  full_model_loglikelihood <- - FHT_minus_loglikelihood_with_all_parameters(result$final_parameters$beta_hat_final, result$final_parameters$gamma_hat_final, X_test, Z_test, times_test, delta_test)
  deviance <- 2 * (full_model_loglikelihood - null_model_loglikelihood)
  # deviance
  full_filename <- make_filename(directory, "deviance", seed)
  write.csv(as.data.frame(deviance), file=full_filename, row.names=FALSE)

  # Brier R2
  brier_r2_result <- brier_r2(times_test, delta_test, estimated_y0s_test, estimated_mus_test, null_y0_test, null_mu_test, number_of_time_points=100)
  full_filename <- make_filename(directory, "brier", seed)
  write.csv(brier_r2_result, file=full_filename, row.names=FALSE)
}

if (do_CV) {
  foreach(seed=seeds) %dopar% {
    run_CV_and_write_to_file(N, setup_type, add_noise, seed, K, K_fold_repetitions, directory, M)
  }
  stopCluster(cl)
}

if (train_models) {
  foreach(seed=seeds) %dopar% {
    estimate_model_and_validate_and_write_to_file(N, setup_type, add_noise, seed, K, K_fold_repetitions, directory, M)
  }
  stopCluster(cl)
}
