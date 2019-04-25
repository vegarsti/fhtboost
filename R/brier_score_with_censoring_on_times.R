#' @export

brier_score_with_censoring_on_times <- function(times, delta, y0s, mus) {

  # Make them plain vectors so as to not have indexes messing up
  # y0s <- as.numeric(y0s)
  # mus <- as.numeric(mus)
  #
  # # Sort incoming
  # order_times <- order(times)
  # delta <- delta[order_times]
  # y0s <- y0s[order_times]
  # mus <- mus[order_times]
  # times <- sort(times)

  censoring_estimates <- sapply(times, function(time) {
    kaplan_meier_estimate_of_censoring(time, times, delta)
  })
  # must assume sorted !!!!!! v
  N <- length(times)
  times <- times[-N]
  delta <- delta[-N]
  y0s <- y0s[-N]
  mus <- mus[-N]
  censoring_estimates <- censoring_estimates[-N]

  estimated_probabilities_ <- matrix(nrow=90, ncol=90)
  N <- length(times)
  final_scores <- rep(NA, N)
  for (i in 1:N) {
    current_time <- times[i]
    estimated_probabilities <- FHT_parametric_survival(current_time, mus, y0s)
    estimated_probabilities_[, i] <- estimated_probabilities
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

