#' @export

run_CV <- function(M, K_fold_repetitions, K, X, Z, times, delta, boost_intercepts_continually=boost_intercepts_continually) {
  CV_errors_K_deviance <- matrix(NA, nrow=M, ncol=K_fold_repetitions)
  CV_errors_K_loglik <- matrix(NA, nrow=M, ncol=K_fold_repetitions)
  all_CV_errors <- list()
  for (repeated_K_fold_iteration in 1:K_fold_repetitions) {
    folds <- create_folds_stratified(delta, K)
    CV_error_matrix_loglik <- matrix(NA, nrow=M, ncol=K)
    CV_error_matrix_deviance <- matrix(NA, nrow=M, ncol=K)
    all_CV_errors[[repeated_K_fold_iteration]] <- list()
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
      result <- boosting_run(times_without_k, delta_without_k, X_without_k, Z_without_k, M, boost_intercepts_continually=boost_intercepts_continually)
      beta_hats <- result$parameters$beta_hats
      gamma_hats <- result$parameters$gamma_hats
      Ms <- 1:M
      # find null model
      best_intercepts <- maximum_likelihood_intercepts(times_without_k, delta_without_k)
      null_y0 <- best_intercepts[1]
      null_mu <- best_intercepts[2]
      beta_null <- c(null_y0, rep(0, (dim(X_k)[2]-1)))
      gamma_null <- c(null_mu, rep(0, (dim(Z_k)[2]-1)))
      null_model_loglikelihood <- - FHT_minus_loglikelihood_with_all_parameters(beta_null, gamma_null, X_k, Z_k, times_k, delta_k)

      # deviance
      deviance <- sapply(Ms, function(m) {
        full_model_loglikelihood <- - FHT_minus_loglikelihood_with_all_parameters(beta_hats[m, ], gamma_hats[m, ], X_k, Z_k, times_k, delta_k)
        deviance <- 2 * (full_model_loglikelihood - null_model_loglikelihood)
        return(deviance)
      })
      CV_error_matrix_deviance[, k] <- deviance

      # loglik
      loglik <- sapply(Ms, function(m) {
        return(FHT_minus_loglikelihood_with_all_parameters(
          beta_hats[m, ], gamma_hats[m, ], X_k, Z_k, times_k, delta_k
        ))
      })
      CV_error_matrix_loglik[, k] <- loglik
      all_CV_errors[[repeated_K_fold_iteration]][[k]] <- list(deviance=deviance, loglik=loglik)
    }
    CV_errors_K_deviance[, repeated_K_fold_iteration] <- rowSums(CV_error_matrix_deviance)
    CV_errors_K_loglik[, repeated_K_fold_iteration] <- rowSums(CV_error_matrix_loglik)
  }
  return(list(
    CV_errors_K_loglik=CV_errors_K_loglik,
    CV_errors_K_deviance=CV_errors_K_deviance,
    all_CV_errors=all_CV_errors
  ))
}
