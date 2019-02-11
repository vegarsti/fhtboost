#' @export

kaplan_meier_plot <- function(times, delta) {
  non_para <- non_parametric_estimates(times, delta, continuous = TRUE)
  plot(non_para$times_sequence, non_para$kaplan_meiers, typ='s', xlab="Time", ylab="Kaplan-Meier estimated survival probability")
}
