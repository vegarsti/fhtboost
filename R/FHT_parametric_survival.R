#' @export

FHT_parametric_survival <- function(times, mu, y0) {
  return(pnorm((mu*times + y0)/sqrt(times)) - exp(-2*y0*mu)*pnorm((mu*times - y0))/sqrt(times))
}
