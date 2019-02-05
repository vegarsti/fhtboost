brier_score_on_uncensored_data <- function(times, delta, estimated_y0, estimated_mu) {
  max_time <- max(times)
  time_points <- 1000
  brier_times <- seq(from=0.1, to=max_time, length.out=time_points)
  observed_time_indexes <- delta == 1
  brier_scores <- sapply(brier_times, function(time) {
    brier_score_no_censoring(time, times[observed_time_indexes], estimated_y0[observed_time_indexes], estimated_mu[observed_time_indexes])
  })
  return(brier_scores)
}
