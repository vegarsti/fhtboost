run_CV <- function(M, K_fold_repetitions, K, X, Z, times, delta, criterion="LOGLIK") {
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
      Ms <- 1:M
      if (criterion == 'deviance') {
        # find null model
        best_intercepts <- maximum_likelihood_intercepts(times_without_k, delta_without_k)
        null_y0 <- best_intercepts[1]
        null_mu <- best_intercepts[2]
        beta_null <- c(null_y0, rep(0, (dim(X_k)[2]-1)))
        gamma_null <- c(null_mu, rep(0, (dim(Z_k)[2]-1)))
        null_model_loglikelihood <- - FHT_minus_loglikelihood_with_all_parameters(beta_null, gamma_null, X_k, Z_k, times_k, delta_k)

        deviances <- sapply(Ms, function(m) {
          beta_m <- result$parameters$beta_hats[m, ]
          gamma_m <- result$parameters$gamma_hat[m, ]
          estimated_y0s_test <- exp(X_k %*% beta_m)
          estimated_mus_test <- Z_k %*% gamma_m
          full_model_loglikelihood <- - FHT_minus_loglikelihood_with_all_parameters(beta_m, gamma_m, X_k, Z_k, times_k, delta_k)
          deviance <- 2 * (full_model_loglikelihood - null_model_loglikelihood)
          return(deviance)
        })
        CV_error_matrix[, k] <- deviances
      } else { #if (criterion == 'LOGLIK') {
        for (m in Ms) {
          gamma_m <- gamma_hats[m, ]
          beta_m <- beta_hats[m, ]
          CV_error_matrix[m, k] <- FHT_minus_loglikelihood_with_all_parameters(
            beta_m, gamma_m, X_k, Z_k, times_k, delta_k
          )
        }
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
