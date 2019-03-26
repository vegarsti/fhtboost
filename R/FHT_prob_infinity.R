#' @export
FHT_prob_infinity <- function(y0, mu) {
  return(1-exp(-2*y0*mu))
}
