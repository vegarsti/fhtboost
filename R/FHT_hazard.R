FHT_hazard <- function(times, mu, y0) {
  return(FHT_parametric_density(times, mu, y0)/FHT_parametric_survival(times, mu, y0))
}
