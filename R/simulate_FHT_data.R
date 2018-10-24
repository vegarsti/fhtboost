#' Simulate FHT data
#'
#' Simulate FHT data
#'
#' @return \code{observations} A list of \code{survival_times} and corresponding \code{delta}
#' @return \code{true_parameters} A list of \code{beta} and \code{gamma}
#' @return \code{design_matrices} A list of \code{X} and \code{Z}
#'
#' @keywords keywords
#'
#' @export
#'
#' @examples
#' result <- generate_test_data()
#' observations <- result$observations
#' times <- observations$survival_times
#' delta <- observations$delta
#' # make

simulate_FHT_data <- function(dense=TRUE) {
  library(statmod)
  set.seed(2)
  N <- 1000

  if (dense) {
    # y0, beta, X
    beta_ <- c(1.5, 0.1, 0.2)
    d <- length(beta_)
    X0 <- rep(1, N)
    X1 <- cbind(c(rep(1, 300), rep(2, 300), rep(-0.5, 400)))
    X1 <- scale(X1)
    X2 <- rnorm(N)
    X_design_matrix <- cbind(X0, X1, X2)
    X_design_matrix[, 2:d] <- X_design_matrix[, 2:d] + rnorm(prod(dim(X_design_matrix[, 2:d])))

    # mu, gamma, Z
    # with intercept and normalization
    gamma_ <- c(-1.0, -0.2, 0.1)
    p <- length(gamma_)
    Z0 <- rep(1, N)
    Z1 <- cbind(c(rep(1, 500), rep(-1, 500)))
    Z2 <- rnorm(N)
    Z_design_matrix <- cbind(Z0, Z1, Z2)
    # add noise
    Z_design_matrix[, 2:p] <- Z_design_matrix[, 2:p] + rnorm(prod(dim(Z_design_matrix[, 2:p])))
  } else {
    beta_ <- c(3, rep(0, 10), rep(0.5, 10))
    gamma_ <- c(-1, rep(0.5, 10), rep(0, 10))

    # beta0 <- 0.1
    # beta_ <- c(beta0, rep(0, 10), rep(0.1, 10))
    # gamma0 <- 0.1
    # gamma_ <- c(-gamma0, rep(-0.05, 10), rep(0, 10))

    d <- length(beta_)
    X0 <- rep(1, N)
    #Xrest <- matrix(rnorm((d-1)*N), ncol=(d-1))
    # rbinom(100, 1, 0.5) -- bernoulli
    Xrest <- matrix(rbeta((d-1)*N, shape1=1, shape2=1), ncol=(d-1))
    X_design_matrix <- cbind(X0, Xrest)

    p <- length(beta_)
    Z0 <- rep(1, N)
    #Zrest <- matrix(rnorm((p-1)*N), ncol=(p-1))
    Zrest <- matrix(rbeta((p-1)*N, shape1=1, shape2=1), ncol=(p-1))
    Z_design_matrix <- cbind(Z0, Zrest)
  }
  y0 <- exp(X_design_matrix %*% beta_)
  print("y0")
  print(y0)
  mu <- Z_design_matrix %*% gamma_
  print("mu")
  print(mu)
  sigma_2 <- 1 ## NB

  # Transform parameters
  mu_IG <- y0/abs(mu)
  lambda_IG <- (y0/sigma_2)^2

  # Draw survival times and censoring times
  survival_times <- statmod::rinvgauss(N, mu_IG, lambda_IG)
  censoring_times <- statmod::rinvgauss(N, mu_IG, 100*lambda_IG)
  #censoring_times <- statmod::rinvgauss(N, 2*mu_IG, lambda_IG)

  #plot(survival_times)
  #points(censoring_times, col='red')

  observations <- censor_observations(survival_times, censoring_times)
  censored_survival_times <- observations$times
  observed <- observations$delta
  observations <- list(survival_times=censored_survival_times, delta=observed)
  true_parameters <- list(beta=beta_, gamma=gamma_)
  design_matrices <- list(X_design_matrix=X_design_matrix, Z_design_matrix=Z_design_matrix)
  return(list(
    observations=observations,
    true_parameters=true_parameters,
    design_matrices=design_matrices
  ))
}
