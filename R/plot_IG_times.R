plot_IG_times <- function() {
  y0 <- 100
  mu <- -1
  obs <- draw_IG_data(y0, mu, N=1000)
  times <- obs$survival_times
  delta <- obs$delta
  non_para <- non_parametric_estimates(times, delta, continuous = TRUE)
  parametric_times <- seq(0.8, max(times), by=0.01)
  parametric_S <- FHT_parametric_survival(parametric_times, mean(mu), mean(y0))
  plot(non_para$times_sequence, non_para$kaplan_meiers, typ='s', ylim=c(0, 1))
  lines(parametric_times, parametric_S, col='red')
  #plot(parametric_times, parametric_S, col='red')
}
