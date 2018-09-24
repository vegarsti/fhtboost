simulate_normal_data <- function() {
  set.seed(2)
  N <- 100
  beta_true <- c(0.5, 1)
  p <- length(beta_true)
  X1 <- cbind(c(rep(2, 50), rep(1, 50)))
  X2 <- cbind(c(rep(0, 50), rep(2, 50)))
  X <- cbind(X1, X2) + rnorm(n=p*N)
  y_true <- X %*% beta_true
  y_observations <- rnorm(length(y_true), mean=y_true, sd=1)

  return(list(
    y=y_observations,
    beta=beta_true,
    X=X
  ))
}
