#' Loss function for \code{mu}
#'
#' Derivative of minus log likelihood with respect to
#' one parameter: \code{mu}
#'
#' @param y0 Initial level of underlying Wiener process \code{y0}
#' @param mu Drift of underlying Wiener process \code{mu}
#'
#' @return gradient The gradient
#'
#' @keywords keywords
#'
#' @export
#'
#' @examples
#' R code here showing how your function works

loss_function_derivative_mu <- function(y0, mu, sigma2, times, delta) {
  # t is vector of times, delta is vector of observed = 1, or not = 0
  log_f <- - (y0 + mu*times)/sigma2
  log_S1 <- (times/sqrt(sigma2*times)) * dnorm((mu*times + y0)/(sqrt(sigma2*times)))
  log_S2 <- (2*y0/sigma2) * exp(-2*y0*mu/sigma2) * pnorm((mu*times-y0)/sqrt(sigma2*times))
  log_S3 <- -(times/sqrt(sigma2*times)) * exp(-2*y0*mu/sigma2)*dnorm((mu*times - y0)/sqrt(sigma2*times))
  S <- survival_censored(y0, mu, sigma2, times)
  log_S <- (log_S1 + log_S2 + log_S3)/S
  total <- ifelse(delta, log_f, log_S)
  return(total)
}
