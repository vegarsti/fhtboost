run_normal_boosting <- function() {
  simulated_data <- simulate_normal_data()
  y <- simulated_data$y
  X <- simulated_data$X

  m_stop <- 100

  nu <- 0.1

  N <- dim(X)[1]
  p <- dim(X)[2]

  beta_hat <- matrix(NA, nrow=m_stop, ncol=p)
  y_hat <- matrix(NA, nrow=m_stop, ncol=N)
  y_hat[1, ] <- rep(mean(y), N)
  beta_hat[1, ] <- c(0, 0)
  residuals <- matrix(NA, nrow=m_stop, ncol=N)
  residuals[1, ] <- y - y_hat[1, ]
  residual_sums <- rep(NA, m_stop)
  residual_sums[1] <- sum(residuals[1,]^2)

  # Standardize X -- use scale instead!!
  X_original <- X
  X_colmeans <- colMeans(X)
  X_sds <- apply(X, 2, sd)
  for (j in 1:p) {
    X[, j] <- (X[, j] - X_colmeans[j])/X_sds[j]
  }

  u_hats_rss <- rep(NA, p)
  u_hats <- matrix(NA, nrow=p, ncol=N)
  beta_hats <- rep(NA, p)
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
  plot(residual_sums)
  beta_hat_final <- apply(beta_hat, 2, sum)
  #cat("Estimated beta hat: \n")
  beta_hat_final
  #cat("True beta: \n")
  #simulated_data$beta
}
