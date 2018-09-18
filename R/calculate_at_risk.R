#' At risk
#'
#' Given survival times and a sequence of (possibly more granular) times, calculate how
#' many is at risk at any given time. An individual is at risk at time t if it has survival time >= t.
#'
#' @param times Vector of time points
#' @param survival_times Survival times (possibly censored)
#'
#' @return \code{at_risk} Vector of number of individuals at risk at each time in the original times.
#'
#' @keywords keywords
#'
#' @export
#'
#' @examples
#' R code here showing how your function works

calculate_at_risk <- function(times, survival_times) {
  at_risk <- sapply(times, function(ti) { sum(survival_times>=ti) })
  return(at_risk)
}
