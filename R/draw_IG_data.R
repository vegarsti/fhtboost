draw_IG_data <- function(y0, mu, N) {
  set.seed(2)
  sigma_2 <- 1
  mu_IG <- y0/(-mu)
  lambda_IG <- y0^2/sigma_2
  # Draw survival times and censoring times
  survival_times <- statmod::rinvgauss(N, mean=mu_IG, shape=lambda_IG)
  censoring_times <- statmod::rinvgauss(N, mean=mu_IG, shape=100*lambda_IG)
  observations <- censor_observations(survival_times, censoring_times)
  censored_survival_times <- observations$times
  observed <- observations$delta
  observations <- list(survival_times=censored_survival_times, delta=observed)
  return(observations)
}
