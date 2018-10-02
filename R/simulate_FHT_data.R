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

simulate_FHT_data <- function() {
  library(statmod)
  set.seed(2)
  N <- 1000


  # y0, beta, X
  beta_ <- c(1.5, 0.1, 0.2)
  X0 <- rep(1, N)
  X1 <- cbind(c(rep(1, 300), rep(2, 300), rep(-0.5, 400)))
  X1 <- scale(X1)
  X2 <- rnorm(N)
  X_design_matrix <- cbind(X0, X1, X2)
  X_design_matrix <- X_design_matrix + rnorm(prod(dim(X_design_matrix)), sd = 0.5)

  # mu, gamma, Z
  # with intercept and normalization
  gamma_ <- c(-1.0, -0.2, 0.1)
  Z0 <- rep(1, N)
  Z1 <- cbind(c(rep(1, 500), rep(-1, 500)))
  Z2 <- rnorm(N)
  Z_design_matrix <- cbind(Z0, Z1, Z2)
  # add noise
  Z_design_matrix <- Z_design_matrix + rnorm(prod(dim(Z_design_matrix)), sd = 0.5)

  # scale
  #Z_design_matrix <- scale(Z_design_matrix, center=T, scale=T)

  # add intercept
  #Z_design_matrix <- cbind(rep(1, N), Z_design_matrix)

  y0 <- exp(X_design_matrix %*% beta_)
  mu <- Z_design_matrix %*% gamma_
  sigma_2 <- 1 ## NB

  # Transform parameters
  mu_IG <- y0/abs(mu)
  lambda_IG <- (y0/sigma_2)^2

  # Draw survival times and censoring times
  survival_times <- statmod::rinvgauss(N, mu_IG, lambda_IG)
  censoring_times <- statmod::rinvgauss(N, mu_IG, lambda_IG*100)

  observations <- censor_observations(survival_times, censoring_times)
  censored_survival_times <- observations$t
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
