#' Censor observations
#'
#' Censor survival times.
#'
#' @param uncensored_survival_times Vector of uncensored survival times
#' @param censoring_times Vector of times to censor
#'
#' @return \code{t} A vector of survival times, censored according to censoring times
#' @return \code{delta} Censoring observation vector
#'
#' @keywords keywords
#'
#' @export
#'
#' @examples
#' Code here.

censor_observations <- function(uncensored_survival_times, censoring_times) {
  censored_survival_times <- uncensored_survival_times
  censored_survival_times[is.na(uncensored_survival_times)] <- censoring_times # fill is.na with censoring time
  censored_survival_times <- ifelse(censored_survival_times < censoring_times, censored_survival_times, censoring_times)
  observed <- ifelse(censored_survival_times < censoring_times, 1, 0)
  times <- pmin(censored_survival_times, censoring_times)
  return(list(t=times, delta=observed))
}