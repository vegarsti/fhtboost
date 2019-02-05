brier_score_on_censored_data <- function(times, delta, estimated_y0s, estimated_mus) {
  time_points <- 1000
  brier_times <- seq(from=0.1, to=max(times), length.out=time_points)
  brier_scores <- sapply(brier_times, function(time) {
    brier_score_with_censoring(time, times, delta, estimated_y0s, estimated_mus)
  })
  return(list(brier_times=brier_times, brier_scores=brier_scores))
}
