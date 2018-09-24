run_boosting <- function() {
  simulated_data <- simulate_FHT_data()
  times <- simulated_data$observations$survival_times
  delta <- simulated_data$observations$delta
  X <- simulated_data$design_matrices$X
  Z <- simulated_data$design_matrices$Z

  # Our loss function
  minus_FHT_loglikelihood <- data_to_FHT_loss_function(X, Z, times, delta)
  # usage: minus_FHT_loglikelihood(list(beta=beta, gamma=gamma))

  m_stop <- 1000

  nu <- 0.1

  N <- dim(X)[1]
  p <- dim(X)[2]

  beta_true <- c(0.5, 1) # the true beta!

  gamma_hat <- matrix(NA, nrow=m_stop, ncol=p)
  y_hat <- matrix(NA, nrow=m_stop, ncol=N)
  y_hat[1, ] <- rep(mean(y), N)
  gamma_hat[1, ] <- c(0, 0)
  losses <- matrix(NA, nrow=m_stop, ncol=N)
  losses[1, ] <- minus_FHT_loglikelihood(list(beta=beta_true, gamma=gamma_hat[1, ]))
  sum_loss <- rep(NA, m_stop)
  sum_loss[1] <- sum(losses[1,]^2)

  # Standardize X -- use scale instead!!
  X_original <- X
  X_colmeans <- colMeans(X)
  X_sds <- apply(X, 2, sd)
  for (j in 1:p) {
    X[, j] <- (X[, j] - X_colmeans[j])/X_sds[j]
  }

  u_hats_rss <- rep(NA, p)
  u_hats <- matrix(NA, nrow=p, ncol=N)
  gamma_hats <- rep(NA, p)
  for (m in 2:m_stop) {
    u <- residuals[(m-1), ]
    for (j in 1:p) {
      X_j <- X[,j]
      S <- X_j %*% solve(t(X_j) %*% X_j) %*% t(X_j)
      beta_hats[j] <- solve(t(X_j) %*% X_j) %*% t(X_j) %*% u
      u_hats[j, ] <- X_j * beta_hats[j]
      #u_hats[j, ] <- S %*% u
      u_hats_rss[j] <- sum((u_hats[j, ] - u)^2) / N
    }
    best_p <- which.min(u_hats_rss)
    beta_hat[m, best_p] <- nu*beta_hats[best_p]
    beta_hat[m, -best_p] <- 0
    y_hat[m, ] <- y_hat[(m-1), ] + X %*% beta_hat[m, ]
    residuals[m, ] <- y - y_hat[m, ]
    residual_sums[m] <- sum(residuals[m,]^2)
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
