FHT_minus_loglikelihood_with_all_parameters <- function(optimization_parameters, X, Z, times, delta) {
  loglikelihood_vector <- FHT_componentwise_minus_loglikelihood_with_parameters(parameter_list, X, Z, times, delta)
  loglikelihood_value <- sum(loglikelihood_vector)
  return(-loglikelihood_value)
}
