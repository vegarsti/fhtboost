#' @export

FHT_loglikelihood_with_y0_mu <- function(y0, mu, times, delta) {
  log_f <- log(y0) - 0.5*log(2*pi*times^3) - ((y0 + mu*times)^2)/(2*times)
  S <- pnorm((mu*times + y0)/sqrt(times)) - exp(-2*y0*mu)*pnorm((mu*times - y0))/sqrt(times)
  log_S <- log(S)
  #loglikelihood_vector <- delta*log_f + (1-delta)*log_S
  loglikelihood_vector <- ifelse(delta, log_f, log_S)
  return(loglikelihood_vector)
}
