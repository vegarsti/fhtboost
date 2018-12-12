maximum_likelihood_intercepts <- function(times, delta) {
  return(nlm(FHT_only_intercepts, c(0.1, 0.1), times, delta)$estimate)
}
