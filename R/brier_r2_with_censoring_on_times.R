#' @export

brier_r2_with_censoring_on_times <- function(times, delta, y0s, mus, y0_null, mu_null) {
  brier_model <- brier_score_with_censoring_on_times(times, delta, y0s, mus)
  brier_null <- brier_score_with_censoring_on_times(times, delta, rep(y0_null, length(y0s)), rep(mu_null, length(mus)))
  return(1 - brier_model/brier_null)
}
