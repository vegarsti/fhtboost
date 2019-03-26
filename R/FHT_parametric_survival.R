#' @export

FHT_parametric_survival <- function(time, mu, y0) {
  return(pnorm((mu*time + y0)/sqrt(time)) - exp(-2*y0*mu)*pnorm((mu*time - y0))/sqrt(time))
}
