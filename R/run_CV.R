run_CV <- function(M, K_fold_repetitions, K, X, Z, times, delta) {
  CV_errors_K <- matrix(NA, nrow=M, ncol=K_fold_repetitions)
  loss_K <- matrix(NA, nrow=M, ncol=K_fold_repetitions)
  for (repeated_K_fold_iteration in 1:K_fold_repetitions) {
    folds <- create_folds_stratified(delta, K)
    CV_error_matrix <- matrix(NA, nrow=M, ncol=K)
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
      delta_k <- delta[subset_k]
      result <- boosting_run(times_without_k, delta_without_k, X_without_k, Z_without_k, M, boost_intercepts_continually=TRUE)
      beta_hats <- result$parameters$beta_hats
      gamma_hats <- result$parameters$gamma_hats
      for (m in 1:M) {
        gamma_m <- gamma_hats[m, ]
        beta_m <- beta_hats[m, ]
        CV_error_matrix[m, k] <- FHT_minus_loglikelihood_with_all_parameters(
          beta_m, gamma_m, X_k, Z_k, times_k, delta_k
        )
      }
    }
    CV_errors_K[, repeated_K_fold_iteration] <- rowSums(CV_error_matrix)
  }
  CV_errors <- rowMeans(CV_errors_K)
  return(list(
    CV_errors=CV_errors,
    CV_errors_K=CV_errors_K
  ))
}
