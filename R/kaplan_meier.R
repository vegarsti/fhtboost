#' Kaplan-Meier
#'
#' Non-parametric (Kaplan-Meier) estimates of the survival curve
#'
#' @param times Survival times (possibly censored)
#' @param delta Corresponding observation vector (0 is censored, 1 is actually observed)
#'
#' @return \code{kaplan_meier} Vector of Kaplan-Meier estimates
#'
#' @keywords keywords
#'
#' @export
#'
#' @examples
#' R code here showing how your function works

kaplan_meier <- function(times, delta) {
  sorted_times <- times[order(times)]
  sorted_delta <- delta[order(times)]
  at_risk <- calculate_at_risk(sorted_times, sorted_times)
  N <- length(times)
  kaplan_meier_vector <- rep(1, N)
  kaplan_meier_vector[2:N] <- kaplan_meier_vector[1:(N-1)] * (1 - sorted_delta[1:(N-1)]/at_risk[1:(N-1)])
  return(kaplan_meier_vector)
}
