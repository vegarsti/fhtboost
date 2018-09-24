run_boosting <- function() {
  simulated_data <- simulate_FHT_data()
  times <- simulated_data$observations$survival_times
  delta <- simulated_data$observations$delta
  X <- simulated_data$design_matrices$X
  Z <- simulated_data$design_matrices$Z

  # Our loss function
  minus_FHT_loglikelihood <- data_to_FHT_loss_function(X, Z, times, delta)
  # usage: minus_FHT_loglikelihood(list(beta=beta, gamma=gamma))

  m_stop <- 100

  N <- dim(X)[1]

  # index_set <- 1:N
  # #index_set <- index_set[-c(146, 298)]
  #
  # X <- X[index_set, ]
  # Z <- Z[index_set, ]
  # times <- times[index_set]
  # delta <- delta[index_set]

  nu <- 0.1

  N <- dim(X)[1]
  p <- dim(X)[2]

  # Standardize X -- use scale instead!! (?)

  # X_original <- X
  # X_colmeans <- colMeans(X)
  # X_sds <- apply(X, 2, sd)
  # for (j in 1:p) {
  #   X[, j] <- (X[, j] - X_colmeans[j])/X_sds[j]
  # }

  beta_true <- c(0.5, 1)       # the true beta!
  gamma_true <- c(-0.2, -0.2)

  gamma_hat <- matrix(NA, nrow=m_stop, ncol=p)
  gamma_hat[1, ] <- c(0, 0)
  gradient <- matrix(NA, nrow=m_stop, ncol=N)
  parameter_list <- list(beta=beta_true, gamma=gamma_hat[1, ])
  #losses[1, ] <- FHT_componentwise_minus_loglikelihood_with_parameters(parameter_list, X, Z, times, delta)
  gradient[1, ] <- - FHT_componentwise_loss_function_derivative_mu(parameter_list, X, Z, times, delta)
  loss <- rep(NA, m_stop)
  loss[1] <- minus_FHT_loglikelihood(parameter_list)

  u_hats_rss <- rep(NA, p)
  u_hats <- matrix(NA, nrow=p, ncol=N)
  gamma_hats <- rep(NA, p)
  for (m in 2:m_stop) {
    u <- gradient[(m-1), ]
    for (j in 1:p) {
      X_j <- X[,j]
      S <- X_j %*% solve(t(X_j) %*% X_j) %*% t(X_j)
      gamma_hats[j] <- solve(t(X_j) %*% X_j) %*% t(X_j) %*% u
      u_hats[j, ] <- X_j * gamma_hats[j]
      #u_hats[j, ] <- S %*% u
      u_hats_rss[j] <- sum((u_hats[j, ] - u)^2) / N
    }
    best_p <- which.min(u_hats_rss)
    gamma_hat[m, best_p] <- nu*gamma_hats[best_p]
    gamma_hat[m, -best_p] <- 0
    parameter_list <- list(beta=beta_true, gamma=colSums(gamma_hat[1:m, ]))
    gradient[m, ] <- - FHT_componentwise_loss_function_derivative_mu(parameter_list, X, Z, times, delta)
    loss[m] <- minus_FHT_loglikelihood(parameter_list)
  }



  link_function <- FHT_link_function(X, Z)
  FHT_parameters <- link_function(beta_hat, gamma_hat)
  y0s <- FHT_parameters$y0
  y0 <- 1
  mus <- FHT_parameters$mu
  sigma2 <- 1
  p <- 2
  N <- 1000

  u <- - sapply(mus, function(mu) { loss_function_derivative_mu(y0, mu, sigma2, times, delta) })

  u_hats <- rep(NA, p)
  for (column in 1:p) {
    X_subset <- X[,column]
    u_hat <- X_subset %*% solve(t(X_subset) %*% X_subset) %*% t(X_subset) %*% u
    rss <- sum((u_hat - u)^2) / N
    u_hats[column] <- rss
  }
  best_p <- which.min(u_hats)

  M <- 1000
  return(loss)
}
