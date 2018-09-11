#' Run optimization on test data
#'
#' Function doing the whole shebang.
#'
#'
#' @keywords keywords
#'
#' @export
#'
#' @examples
#' R code here showing how your function works

run_optimization_on_test_data <- function() {
  result <- generate_test_data()
  observations <- result$observations
  t <- observations$survival_times
  delta <- observations$delta
  X <- result$design_matrices$X
  Z <- result$design_matrices$Z
  optimizable_function <- data_to_optimizable_function(X, Z, t, delta)
  initial_parameters <- runif(4, min=1, max=2)
  optimized_parameters <- nlm(optimizable_function, initial_parameters)$estimate
  beta_hat <- optimized_parameters[1:2]
  gamma_hat <- optimized_parameters[3:4]
  beta <- result$true_parameters$beta
  gamma <- result$true_parameters$gamma
  print(beta_hat)
  print(beta)
  print(gamma_hat)
  print(gamma)
}
