#' @export

cyclic_boosting_run <- function(times, delta, X, Z, m_stop_y0, m_stop_mu, boost_intercepts_continually, should_print=FALSE) {
  N_X <- dim(X)[1]
  N_Z <- dim(Z)[1]
  N_t <- length(times)
  N_delta <- length(delta)
  dimensions <- c(N_X, N_Z, N_t, N_delta)
  all_dimensions_are_equal <- var(dimensions) == 0
  if (! all_dimensions_are_equal) {
    cat('Dimensions:\n')
    cat('X: ')
    cat(N_X)
    cat('Z: ')
    cat(N_Z)
    cat('times: ')
    cat(N_t)
    cat('delta: ')
    cat(N_delta)
    stop("Not all dimensions are equal!")
  } else { # all dimensions are equal
    N <- N_X
  }
  d <- dim(X)[2]
  p <- dim(Z)[2]

  ps <- 2:p
  ds <- 2:d

  # Nuisance parameter: shrinkage factor
  nu <- 0.1

  max_m_stop <- max(m_stop_y0, m_stop_mu)

  nlm_result_estimate <- maximum_likelihood_intercepts(times, delta)

  beta_0 <- nlm_result_estimate[1]
  gamma_0 <- nlm_result_estimate[2]

  gamma_hat <- matrix(NA, nrow=max_m_stop, ncol=p)
  gamma_hat[1, ] <- rep(0, p)
  gamma_hat_cumsum <- matrix(NA, nrow=max_m_stop, ncol=p)

  beta_hat <- matrix(NA, nrow=max_m_stop, ncol=d)
  beta_hat[1, ] <- rep(0, d)
  beta_hat_cumsum <- matrix(NA, nrow=max_m_stop, ncol=d)

  # INITIALIZE
  gamma_hat[1, 1] <- gamma_0
  beta_hat[1, 1] <- beta_0
  gamma_hat_cumsum[1, ] <- gamma_hat[1, ]
  beta_hat_cumsum[1, ] <- beta_hat[1, ]

  negative_gradient_mu <- matrix(NA, nrow=max_m_stop, ncol=N)
  negative_gradient_y0 <- matrix(NA, nrow=max_m_stop, ncol=N)
  negative_gradient_mu[1, ] <- FHT_componentwise_loss_function_derivative_mu(beta_hat[1, ], gamma_hat[1, ], X, Z, times, delta)
  negative_gradient_y0[1, ] <- FHT_componentwise_loss_function_derivative_y0(beta_hat[1, ], gamma_hat[1, ], X, Z, times, delta)
  loss <- rep(NA, max_m_stop)
  parameter_list <- list(beta=beta_hat[1, ], gamma_hat[1, ])

  gamma_hats <- rep(NA, p)
  beta_hats <- rep(NA, d)
  u_hats <- matrix(NA, nrow=p, ncol=N)
  u_hats_rss <- rep(NA, p)
  loss <- rep(0, max_m_stop)

  for (m in 2:max_m_stop) {
    # gamma/mu
    if (should_print) { cat('iteration: ', m, '\n') }

    if (m <= m_stop_y0) {
      u_y0 <- negative_gradient_y0[m-1, ]
      # BOOST y0
      result_y0 <- boosting_iteration_y0(
        nu=nu,
        X=X,
        Z=Z,
        u=u_y0,
        beta_hat_m1=beta_hat_cumsum[m-1, ],
        gamma_hat_m1=gamma_hat_cumsum[m-1, ],
        d=d,
        ds=ds,
        times=times,
        delta=delta
      )
      beta_hat[m, ] <- result_y0$beta_hat_addition
      beta_hat_cumsum[m, ] <- result_y0$beta_hat_m
    } else {
      beta_hat[m, ] <- beta_hat[m-1, ]
      beta_hat_cumsum[m, ] <- beta_hat_cumsum[m-1, ]
    }
    if (m <= m_stop_mu) {
      u_mu <- FHT_componentwise_loss_function_derivative_mu(beta_hat[m, ], gamma_hat[m-1, ], X, Z, times, delta)
      # BOOST MU
      result_mu <- boosting_iteration_mu(
        nu=nu,
        X=X,
        Z=Z,
        u=u_mu,
        beta_hat_m1=beta_hat_cumsum[m-1, ],
        gamma_hat_m1=gamma_hat_cumsum[m-1, ],
        p=p,
        ps=ps,
        times=times,
        delta=delta
      )
      gamma_hat[m, ] <- result_mu$gamma_hat_addition
      gamma_hat_cumsum[m, ] <- result_mu$gamma_hat_m
    } else {
      gamma_hat[m, ] <- gamma_hat[m-1, ]
      gamma_hat_cumsum[m, ] <- gamma_hat_cumsum[m-1, ]
    }

    negative_gradient_mu[m, ] <- FHT_componentwise_loss_function_derivative_mu(
      beta_hat_cumsum[m, ], gamma_hat_cumsum[m, ], X, Z, times, delta
    )
    negative_gradient_y0[m, ] <- FHT_componentwise_loss_function_derivative_mu(
      beta_hat_cumsum[m, ], gamma_hat_cumsum[m, ], X, Z, times, delta
    )

    if (should_print) {
      cat('gamma0: ', gamma_hat_cumsum[m, 1], '\n')
      cat('beta0: ', beta_hat_cumsum[m, 1], '\n')
    }

  }

  parameters <- list(
    gamma_hats=gamma_hat_cumsum, beta_hats=beta_hat_cumsum
  )
  gamma_hat_final <- gamma_hat_cumsum[m, ]
  beta_hat_final <- beta_hat_cumsum[m, ]
  final_parameters <- list(gamma_hat_final=gamma_hat_final, beta_hat_final=beta_hat_final)
  for (m in 1:max_m_stop) {
    loss[m] <- FHT_minus_loglikelihood_with_all_parameters(
      beta=beta_hat_cumsum[m, ],
      gamma=gamma_hat_cumsum[m, ],
      X, Z, times, delta
    )
  }
  return(list(final_parameters=final_parameters, parameters=parameters, loss=loss))
}
