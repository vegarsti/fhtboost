#' Run optimization on test data
#'
#' Some description
#'
#' @keywords keywords
#'
#' @export
#'
#' @examples
#' none

run_optimization_on_test_data <- function() {
  result <- generate_test_data()
  observations <- result$observations
  times <- observations$survival_times
  delta <- observations$delta
  minus_loglikelihood <- get_minus_loglikelihood_given_observations(times, delta)
  initial_parameters <- runif(4, min=1, max=2)
  optimized_parameters <- nlm(minus_loglikelihood, initial_parameters)$estimate
  beta_hat <- optimized_parameters[1:2]
  gamma_hat <- optimized_parameters[3:4]
  beta <- result$true_parameters$beta
  gamma <- result$true_parameters$gamma
  print(beta_hat)
  print(beta)
  print(gamma_hat)
  print(gamma)
}
