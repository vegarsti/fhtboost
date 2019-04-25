#' @export

kaplan_meier_estimate_of_censoring <- function(time, times, delta) {
  indexes_to_look_at <- which(times <= time)
  if (sum(indexes_to_look_at) == 0) {
    return(1)
  } else {
    times_to_look_at <- times[indexes_to_look_at]
    deltas_to_look_at <- delta[indexes_to_look_at]
    at_risk <- calculate_at_risk(times_to_look_at, times)
    terms <- 1 - (1-deltas_to_look_at)/at_risk
    return(prod(terms))
  }
}
