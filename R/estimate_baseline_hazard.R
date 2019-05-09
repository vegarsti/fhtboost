#' @export

estimate_baseline_hazard <- function(times, delta, linear_predictors) {
  # times, delta and linear_predictors must be sorted in correct time order
  N <- length(times)
  jumps <- rep(0, N)
  num_events <- rep(0, N)
  denominator <- rep(0, N)
  exp_lp <- exp(linear_predictors)
  for (i in 1:N) {
    current_time <- times[i]
    at_risk_indicator <- current_time <= times
    denominator[i] <- sum(at_risk_indicator * exp_lp)
    is_event <- delta[i]
    jumps[i] <- is_event/denominator[i]
  }
  A0 <- cumsum(jumps)
  return(A0)
}
