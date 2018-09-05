#' Generate test data
#'
#' Generate simulated test data: Survival times and observations
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

generate_test_data <- function() {
  library(statmod)
  set.seed(1)
  N <- 1000
  d <- 1
  p <- 1
  beta <- c(2, -1)
  gamma <- c(-2, -1)
  X <- cbind(c(rep(1, 500), rep(0, 500)))
  Z <- cbind(c(rep(1, 200), rep(4, 300), rep(10, 500)))
  y_intercepts <- rep(1, N)
  X_design_matrix <- cbind(y_intercepts, X)
  Z_design_matrix <- cbind(y_intercepts, Z)
  y0 <- exp(X_design_matrix %*% beta)
  mu <- Z_design_matrix %*% gamma
  sigma_2 <- 1 ## NB
  mu_IG <- - y0/mu
  lambda_IG <- (y0/sigma_2)^2
  survival_times <- rinvgauss(N, mu_IG, lambda_IG)
  #censoring_times <- rinvgauss(N, mu_IG, 10*lambda_IG)
  censoring_times <- rep(0.7, N)
  censored_survival_times <- survival_times
  censored_survival_times[is.na(survival_times)] <- censoring_times # censor!
  observed <- ifelse(censored_survival_times < censoring_times, 1, 0)
  times <- pmin(survival_times, censoring_times)
  observations <- list(survival_times=censored_survival_times, delta=observed)
  true_parameters <- list(beta=beta, gamma=gamma)
  design_matrices <- list(X_design_matrix=X_design_matrix, Z_design_matrix=Z_design_matrix)
  return(list(
    observations=observations,
    true_parameters=true_parameters,
    design_matrices=design_matrices
  ))
}
