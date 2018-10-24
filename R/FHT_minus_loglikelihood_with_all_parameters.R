FHT_minus_loglikelihood_with_all_parameters <- function(beta_, gamma_, X, Z, times, delta) {
  minus_loglikelihood_vector <- FHT_componentwise_minus_loglikelihood_with_parameters(beta_, gamma_, X, Z, times, delta)
  minus_loglikelihood_value <- sum(minus_loglikelihood_vector)
  return(minus_loglikelihood_value)
}
