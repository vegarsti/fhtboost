#' @export

brier_r2 <- function(times, delta, estimated_y0s, estimated_mus, null_y0, null_mu, number_of_time_points=100) {
  non_null <- brier_score_on_censored_data(times, delta, estimated_y0s, estimated_mus, number_of_time_points=number_of_time_points)
  null <- brier_score_on_censored_data(times, delta, null_y0, null_mu, number_of_time_points=number_of_time_points)
  r2 <- 1 - non_null$brier_score/null$brier_score
  return(data.frame(times=non_null$brier_times, r2=r2, brier_score=non_null$brier_score, brier_null_score=null$brier_score))
}
