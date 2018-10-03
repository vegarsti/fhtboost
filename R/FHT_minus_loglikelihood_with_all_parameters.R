FHT_minus_loglikelihood_with_all_parameters <- function(beta_, gamma_, X, Z, times, delta) {
  loglikelihood_vector <- FHT_componentwise_minus_loglikelihood_with_parameters(beta_, gamma_, X, Z, times, delta)
  loglikelihood_value <- sum(loglikelihood_vector)
  return(-loglikelihood_value)
}
