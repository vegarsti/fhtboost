boosting_iteration_both <- function(
  nu, X, Z, u_y0, u_mu, beta_hat_m1, gamma_hat_m1, d, ds, p, ps, times, delta,
  X_scale_factors, Z_scale_factors, X_means, Z_means, should_print=FALSE, should_destandardize=TRUE, run_in_parallel=FALSE,iteration_number=1000
) {
  # d corresponds to X, p to Z
  y0 <- exp(X %*% beta_hat_m1)
  result_y0 <- best_least_squares_update(X, u_y0, d, ds, run_in_parallel=run_in_parallel)
  result_mu <- best_least_squares_update(Z, u_mu, p, ps, run_in_parallel=run_in_parallel)

  method <- 'inner' # or 'outer'
  method <- 'outer'

  ### INNER
  if (method == 'inner') {
    rsses <- c(result_y0$rss, result_mu$rss)
    best_rss <- which.min(rsses)
    boosted_mu <- best_rss == 2
  } else if (method == 'outer') { ### OUTER
    # original matrices
    if (should_destandardize) {
      X_original <- destandardize(X, X_scale_factors, X_means)
      Z_original <- destandardize(Z, Z_scale_factors, Z_means)
    } else {
      X_original <- X
      Z_original <- Z
    }

    gamma_hat_addition <- nu*result_mu$parameter_updates
    gamma_hat_m <- gamma_hat_m1 + gamma_hat_addition
    beta_hat_addition <- 0
    beta_hat_m <- beta_hat_m1

    # scale back
    if (should_destandardize) {
      gamma_hat_m_scaled_back <- destandardize(gamma_hat_m, Z_scale_factors, Z_means)
      beta_hat_m_scaled_back <- destandardize(beta_hat_m, X_scale_factors, X_means)
    } else {
      gamma_hat_m_scaled_back <- gamma_hat_m
      beta_hat_m_scaled_back <- beta_hat_m
    }


    # evaluate loss function
    gamma_loss <- FHT_minus_loglikelihood_with_all_parameters(
      beta=beta_hat_m_scaled_back,
      gamma=gamma_hat_m_scaled_back,
      X=X_original,
      Z=Z_original,
      times=times,
      delta=delta
    )

    beta_hat_addition <- nu*result_y0$parameter_updates
    beta_hat_m <- beta_hat_m1 + beta_hat_addition
    gamma_hat_addition <- 0
    gamma_hat_m <- gamma_hat_m1

    # scale back
    if (should_destandardize) {
      gamma_hat_m_scaled_back <- destandardize(gamma_hat_m, Z_scale_factors, Z_means)
      beta_hat_m_scaled_back <- destandardize(beta_hat_m, X_scale_factors, X_means)
    } else {
      gamma_hat_m_scaled_back <- gamma_hat_m
      beta_hat_m_scaled_back <- beta_hat_m
    }

    # evaluate loss function
    beta_loss <- FHT_minus_loglikelihood_with_all_parameters(
      beta=beta_hat_m_scaled_back,
      gamma=gamma_hat_m_scaled_back,
      X=X_original,
      Z=Z_original,
      times=times,
      delta=delta
    )

    # choose best wrt loss ("outer")
    rsses <- c(beta_loss, gamma_loss)
    best_loss <- which.min(rsses)
    boosted_mu <- best_loss == 2
  } else {
    stop("Invalid method for choosing learner!")
  }

  if (should_print) {
    cat("boosted mu", boosted_mu, "\n")
    cat('rss: beta, gamma: ', rsses, '\n')
    cat("beta_hat_m", sum(abs(beta_hat_m)), '\n')
    cat("gamma_hat_m", sum(abs(gamma_hat_m)), '\n')
  }

  ### OUTER

  if (is.null(boosted_mu) || length(boosted_mu) == 0 || iteration_number == 49) {
    # print() should print some diagnostic message
    cat("beta ", beta_loss, "\n")
    cat("gamma ", gamma_loss, "\n")
    print(iteration_number)
    stop('boosted mu is NULL. something went wrong')
  }

  if (boosted_mu) {
    # mu; gamma
    gamma_hat_addition <- nu*result_mu$parameter_updates
    gamma_hat_m <- gamma_hat_m1 + gamma_hat_addition
    beta_hat_addition <- 0
    beta_hat_m <- beta_hat_m1
  } else {
    # y0; beta
    beta_hat_addition <- nu*result_y0$parameter_updates
    beta_hat_m <- beta_hat_m1 + beta_hat_addition
    gamma_hat_addition <- 0
    gamma_hat_m <- gamma_hat_m1
  }

  if (should_print) {
    cat('gamma index: ', which(gamma_hat_addition != 0), '\n')
    cat('gamma hat addition: ', gamma_hat_addition[gamma_hat_addition != 0], '\n')
    cat('beta index: ', which(beta_hat_addition != 0), '\n')
    cat('beta hat addition: ', beta_hat_addition[beta_hat_addition != 0], '\n')
  }

  return(list(
    beta_hat_m=beta_hat_m, beta_hat_addition=beta_hat_addition,
    gamma_hat_m=gamma_hat_m, gamma_hat_addition=gamma_hat_addition,
    boosted_mu=boosted_mu, rsses=rsses
  ))
}
