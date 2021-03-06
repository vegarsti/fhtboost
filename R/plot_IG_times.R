#' @export

plot_IG_times <- function(y0, mu) {
  obs <- draw_IG_data(y0, mu, N=1000)
  times <- obs$survival_times
  delta <- obs$delta
  non_para <- non_parametric_estimates(times, delta, continuous = TRUE)
  parametric_times <- seq(0.8, max(times), by=0.01)
  parametric_S <- FHT_parametric_survival(parametric_times, mean(mu), mean(y0))
  #plot(non_para$times_sequence, non_para$kaplan_meiers, typ='s', ylim=c(0, 1))
  #lines(parametric_times, parametric_S, col='red')
  #plot(parametric_times, parametric_S, col='red')


  times <- seq(from=0.01, to=500, by=0.01)

  mu <- -1
  y0s <- c(50, 20, 10, 5, 2.5)
  hazards <- c()
  for (y0 in y0s) {
    hazards <- cbind(hazards, FHT_hazard(times, mu, y0))
  }

  max_hazards = apply(hazards, 2, max)
  max_index <- which.max(max_hazards)

  plot(times, hazards[, max_index], typ='l')
  for (i in 1:length(y0s)) {
    if (i != max_index) {
      lines(times, hazards[, i], typ='l')
    }
  }
}
