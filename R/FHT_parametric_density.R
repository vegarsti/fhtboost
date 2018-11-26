FHT_parametric_density <- function(times, mu, y0) {
  sigma2 <- 1
  return((y0/sqrt(2*pi*sigma2*times^3)*exp(-((y0+mu*times)^2/(2*sigma2*times)))))
}
