#' @export

brier_score_with_censoring_on_times_with_probabilities <- function(times, delta, estimated_probabilities_matrix) {
  # estimated_probabilities_matrix is a matrix, with row i being probabilities of individuals i, and columns the times
  # estimated_probabilities[i, j] contains the estimated probability based
  # on the model from individual j, that individual i will die


  censoring_estimates <- sapply(times, function(time) {
    kaplan_meier_estimate_of_censoring(time, times, delta)
  })
  # must assume sorted !!!!!!

  # Remove last observation
  N <- length(times)
  times <- times[-N]
  delta <- delta[-N]
  censoring_estimates <- censoring_estimates[-N]
  estimated_probabilities_matrix <- estimated_probabilities_matrix[-N, -N]

  N <- length(times)
  final_scores <- rep(NA, N)
  for (i in 1:N) {
    current_time <- times[i]
    estimated_probabilities <- estimated_probabilities_matrix[, i]
    before_or_at_indicator <- times <= current_time
    after_indicator <- times > current_time
    term1 <- (estimated_probabilities^2*(before_or_at_indicator & delta))/censoring_estimates
    term2 <- ((1-estimated_probabilities)^2*after_indicator)/censoring_estimates[i]
    terms <- term1 + term2
    final_scores[i] <- mean(terms)
  }
  return(
    data.frame(times=times, brier_scores=final_scores)
  )
}

