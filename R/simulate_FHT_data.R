#' Simulate FHT data
#'
#' Simulate FHT data
#'
#' @return \code{list} A list of lots of things
#' @return \code{list} A list of lots of things
#' among others this and this
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
  beta_ <- c(0.5, 1)
  gamma_ <- c(-0.2, -0.2)
  X1 <- cbind(c(rep(2, 500), rep(1, 500)))
  X2 <- cbind(c(rep(0, 500), rep(2, 500)))
  Z1 <- cbind(c(rep(1, 200), rep(4, 300), rep(2, 500)))
  Z2 <- cbind(c(rep(1, 200), rep(0, 300), rep(1, 500)))
  X_design_matrix <- cbind(X1, X2)
  X_design_matrix <- X_design_matrix + rnorm(prod(dim(X_design_matrix)), sd = 0.5)
  Z_design_matrix <- cbind(Z1, Z2)
  Z_design_matrix <- Z_design_matrix + rnorm(prod(dim(Z_design_matrix)), sd = 0.5)
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
