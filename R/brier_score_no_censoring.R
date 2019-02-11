#' @export

brier_score_no_censoring <- function(current_time, times, y0s, mus) {
  estimated_probabilities <- FHT_parametric_survival(current_time, mus, y0s)
  earlier <- times[times < current_time]
  after <- times[times >= current_time]
  before_indicator <- times < current_time
  after_indicator <- 1 - before_indicator
  return(mean(estimated_probabilities^2*before_indicator + (1-estimated_probabilities)^2*after_indicator))
}
