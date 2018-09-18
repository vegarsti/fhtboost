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
  # Generate test data
  result <- generate_test_data()
  observations <- result$observations
  times <- observations$survival_times
  delta <- observations$delta
  X <- result$design_matrices$X
  Z <- result$design_matrices$Z
  minus_FHT_loglikelihood <- data_to_optimizable_function(X, Z, times, delta)

  p <- dim(X)[2]
  d <- dim(Z)[2]

  # Run optimization
  initial_parameters <- runif(p+d, min=1, max=2)

  # Approach: Functional
  #print(minus_FHT_loglikelihood(initial_parameters))
  nlm_result <- nlm(minus_FHT_loglikelihood, initial_parameters)

  # Approach: Put everything into one function
  #print(FHT_minus_loglikelihood_with_all_parameters(initial_parameters, X, Z, times, delta))
  #nlm_result <- nlm(FHT_minus_loglikelihood_with_all_parameters, initial_parameters, X, Z, times, delta)

  optimized_parameters <- nlm_result$estimate

  # Checking etc.
  parameter_list <- parameter_vector_to_list(optimized_parameters, p, d)
  beta_hat <- parameter_list$beta
  gamma_hat <- parameter_list$gamma
  beta_ <- result$true_parameters$beta
  gamma_ <- result$true_parameters$gamma
  return(list(
    t=times,
    delta=delta,
    initial_parameters=initial_parameters,
    beta_hat=beta_hat,
    beta=beta_,
    gamma_hat=gamma_hat,
    gamma=gamma_
  ))
}
