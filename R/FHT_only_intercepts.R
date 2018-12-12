FHT_only_intercepts <- function(parameters, times, delta) {
  y0_baseline <- exp(parameters[1])
  mu_baseline <- parameters[2]
  return(-sum(FHT_loglikelihood_with_y0_mu(y0_baseline, mu_baseline, times, delta)))
}
