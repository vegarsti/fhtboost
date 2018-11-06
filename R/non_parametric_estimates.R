#' Run optimization on simulated data
#'
#' Function doing the whole shebang.
#'
#' @keywords keywords
#'
#' @export
#'
#' @examples
#' R code here showing how your function works

non_parametric_estimates <- function(survival_times, observed, continuous=FALSE) {
  sorted_t_vector <- sort(survival_times)
  diffs <- diff(sorted_t_vector)
  step_size <- 0.9*min(diffs[diffs > 0])
  if (!continuous) {
    times_sequence <- seq(0, max(survival_times), by=step_size)
  } else {
    times_sequence <- sorted_t_vector
  }
  at_risk <- calculate_at_risk(times_sequence, sorted_t_vector)
  sorted_observed <- observed[order(survival_times)]
  N <- length(times_sequence)
  jumps <- rep(0, N)
  jump_prods <- rep(1, N)
  sigma_2_hat_jumps <- rep(0, N)
  for (i in 2:N) {
    num_events <- sum((times_sequence[i-1] < sorted_t_vector & sorted_t_vector <= times_sequence[i]) * sorted_observed)
    jumps[i] <- num_events/at_risk[i]
    jump_prods[i] <- (1 - num_events/at_risk[i])
    if (num_events > 0) {
      z <- sum(1 / (at_risk[i] - (0:(num_events-1)))^2)
    }
    else {
      z <- 0
    }
    #sigma_2_hat_jumps[i] <- z
    sigma_2_hat_jumps[i] <- num_events/(at_risk[i]^2)
  }
  nelson_aalens <- cumsum(jumps)
  kaplan_meiers <- cumprod(jump_prods)
  sigma_2_hat <- cumsum(sigma_2_hat_jumps)
  tau_hat <- kaplan_meiers * sqrt(sigma_2_hat)
  alpha2 <- 0.5
  kaplan_meier_lower <- kaplan_meiers^exp(-alpha2*tau_hat/(kaplan_meiers * log(kaplan_meiers)))
  kaplan_meier_upper <- kaplan_meiers^exp(alpha2*tau_hat/(kaplan_meiers * log(kaplan_meiers)))
  return(list(
    nelson_aalens=nelson_aalens,
    kaplan_meiers=kaplan_meiers,
    kaplan_meier_lower=kaplan_meier_lower,
    kaplan_meier_upper=kaplan_meier_upper,
    times_sequence=times_sequence,
    sigma_2_hat=sigma_2_hat,
    tau_hat=tau_hat,
    at_risk=at_risk,
    sigma_2_hat_jumps=sigma_2_hat_jumps
  ))
}
