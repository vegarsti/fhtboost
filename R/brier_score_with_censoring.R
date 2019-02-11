#' @export

brier_score_with_censoring <- function(current_time, times, delta, y0s, mus) {
  estimated_probabilities <- FHT_parametric_survival(current_time, mus, y0s)
  before_indicator <- times < current_time
  after_indicator <- 1 - before_indicator
  estimates1 <- sapply(times, function(time) {
    kaplan_meier_estimate_of_censoring(time, times, delta)
  })
  term1 <- (estimated_probabilities^2*(before_indicator & delta))/estimates1
  estimate2 <- kaplan_meier_estimate_of_censoring(current_time, times, delta)
  term2 <- ((1-estimated_probabilities)^2*after_indicator)/estimate2
  terms <- term1 + term2
  return(mean(terms))
}
