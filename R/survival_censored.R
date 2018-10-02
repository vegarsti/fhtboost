survival_censored <- function(y0, mu, sigma2, times) {
  return(pnorm((mu*times + y0)/sqrt(sigma2*times)) - exp(-2*y0*mu/sigma2)*pnorm((mu*times - y0)/sqrt(sigma2*times)))
}
