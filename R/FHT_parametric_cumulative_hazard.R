#' @export

FHT_parametric_cumulative_hazard <- function(times, mu, y0) {
  return(-log(pnorm((mu*times + y0)/sqrt(times)) - exp(-2*y0*mu)*pnorm((mu*times - y0))/sqrt(times)))
}
