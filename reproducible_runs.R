# Install package
install_github("vegarsti/fhtboost")

# Parallel libraries -- can be removed if not running in parallel, see below
library(foreach)
library(doParallel)

# Load data
oberthur_filename <- 'preproc_Oberthur_data.Rdata'
load(oberthur_filename)
has_age_observations <- which(!is.na(clinicalData[, 4]))

X <- as.matrix(scale(molecularData[has_age_observations, ]))
Z <- as.matrix(scale(clinicalData[has_age_observations, c(3, 4)]))
times <- clinicalData$time[has_age_observations]
delta <- clinicalData$status[has_age_observations]

# Remove the loaded data
rm(molecularData, clinicalData)


# Options for boosting runs
K_fold_repetitions <- 10
K <- 5
boost_intercepts_continually <- FALSE
boosting_type <- "both" # means boost both parameters, i.e. the intercept and the drift

# Set up parallel settings
no_cores <- detectCores() - 1
registerDoParallel(cores=no_cores)
cl <- makeCluster(no_cores)

# Run 100 seeds in parallel
seeds <- 1:100
foreach(seed=seeds) %dopar% {
  set.seed(seed)
  seed_string <- formatC(seed, width=3, flag="0")

  # Divide into test and train. test approx 1/3
  num_folds <- 3
  folds <- create_folds_stratified(delta, num_folds)
  test_indices <- sort(folds[[1]])
  train_indices <- sort(c(folds[[2]], folds[[3]]))

  ## TRAIN
  ones_train <- rep(1, length(train_indices))
  times_train <- times[train_indices]
  delta_train <- delta[train_indices]
  X_train_rest <- X[train_indices, ]
  X_train <- as.matrix(cbind(ones_train, X_train_rest))
  Z_train_rest <- Z[train_indices, ]
  Z_train <- as.matrix(cbind(ones_train, Z_train_rest))

  ## TEST, and sort
  times_test <- times[test_indices]
  delta_test <- delta[test_indices]
  X_test_rest <- X[test_indices, ]
  Z_test_rest <- Z[test_indices, ]
  order_times <- order(times_test)
  times_test <- sort(times_test)
  delta_test <- delta_test[order_times]
  X_test_rest <- X_test_rest[order_times, ]
  Z_test_rest <- Z_test_rest[order_times, ]

  ones_test <- rep(1, length(test_indices))
  X_test <- as.matrix(cbind(ones_test, X_test_rest))
  Z_test <- as.matrix(cbind(ones_test, Z_test_rest))

  M <- 100

  # Run cross-validation, function from fhtboost package
  CV_result <- run_CV(
    M, K_fold_repetitions, K, X_train, Z_train, times_train, delta_train,
    boost_intercepts_continually=boost_intercepts_continually
  )

  # Write results to file
  directory <- "./" # or some other directory
  full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "loglik.csv")
  write.csv(CV_result$CV_errors_K_loglik, file=full_filename, row.names=FALSE)
  full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "deviance.csv")
  write.csv(CV_result$CV_errors_K_deviance, file=full_filename, row.names=FALSE)
  logliks <- CV_result$CV_errors_K_loglik
  m_stop_from_CV <- which.min(rowMeans(logliks))


  # Run the resulting boosting model on the full training data
  result <- boosting_run(
    times=times_train,
    delta=delta_train,
    X=X_train,
    Z=Z_train,
    m_stop=m_stop_from_CV,
    boost_intercepts_continually=boost_intercepts_continually,
    should_print=FALSE
  )

  # Get the resulting parameter vectors
  beta_hat <- result$final_parameters$beta_hat_final
  gamma_hat <- result$final_parameters$gamma_hat_final
  y0_hat <- exp(X_train %*% beta_hat)
  mu_hat <- Z_train %*% gamma_hat

  # Write resulting data to file
  betas <- data.frame(cbind(non_null_parameters(beta_hat) - 1, beta_hat[non_null_parameters(beta_hat)]))
  names(betas) <- c("j", "beta_j")
  full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "beta.csv")
  write.csv(betas, file=full_filename, row.names=FALSE)
  gammas <- data.frame(cbind((1:length(gamma_hat)) - 1, gamma_hat))
  names(gammas) <- c("j", "gamma_j")
  full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "gamma.csv")
  write.csv(gammas, file=full_filename, row.names=FALSE)


  # Run on test set
  beta_hat_null <- rep(0, dim(X_test)[2])
  beta_hat_null[1] <- beta_hat[1]
  gamma_hat_null <- rep(0, dim(Z_test)[2])
  gamma_hat_null[1] <- gamma_hat[1]
  test_null_loglikelihood <- FHT_minus_loglikelihood_with_all_parameters(
    beta_hat_null, gamma_hat_null, X_test, Z_test, times_test, delta_test
  )
  test_loglikelihood <- FHT_minus_loglikelihood_with_all_parameters(
    beta_hat, gamma_hat, X_test, Z_test, times_test, delta_test
  )
  test_difference_of_deviance <- 2*(test_loglikelihood - test_null_loglikelihood)

  # Calculate Brier score
  y0_hat <- as.numeric(exp(X_test %*% beta_hat))
  mu_hat <- as.numeric(Z_test %*% gamma_hat)
  y0_null <- rep(exp(beta_hat_null)[1], length(y0_hat))
  mu_null <- rep(gamma_hat_null[1], length(mu_hat))

  # Estimated probabilities for the resulting model on the test set
  estimated_probabilities <- sapply(times_test, function(current_time) {
    FHT_parametric_survival(current_time, mu_hat, y0_hat)
  })
  brier_score_df <- brier_score_with_censoring_on_times_with_probabilities(
    times=times_test, delta=delta_test,
    estimated_probabilities_matrix=estimated_probabilities
  )

  # Estimated probabilities of the null model on the test set
  estimated_probabilities_null <- sapply(times_test, function(current_time) {
    FHT_parametric_survival(current_time, mu_null, y0_null)
  })
  brier_null <- brier_score_with_censoring_on_times_with_probabilities(
    times=times_test, delta=delta_test,
    estimated_probabilities_matrix=estimated_probabilities_null
  )
  brier_df_both <- data.frame(
    times=brier_score_df$times,
    brier_scores_null=brier_null$brier_scores,
    brier_scores_model=brier_score_df$brier_scores
  )
  full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "brier_data.csv")
  write.csv(brier_df_both, file=full_filename, row.names=FALSE)
}
