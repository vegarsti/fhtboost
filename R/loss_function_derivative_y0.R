#' Loss function for \code{y0}
#'
#' Derivative of minus log likelihood with respect to
#' one parameter: \code{y0}
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

loss_function_derivative_y0 <- function(y0, mu, sigma2, times, delta) {
  # t is vector of times, delta is vector of observed = 1, or not = 0
  log_f <- 1/y0 - (y0 + mu*times)/(sigma2)#*times)
  log_S1 <- (1/sqrt(sigma2*times)) * dnorm((mu*times + y0)/(sqrt(sigma2*times)))
  log_S2 <- (2*mu/sigma2) * exp(-2*y0*mu/sigma2) * pnorm((mu*times-y0)/sqrt(sigma2*times))
  log_S3 <- 1/sqrt(sigma2*times) * exp(-2*y0*mu/sigma2)*dnorm((mu*times - y0)/sqrt(sigma2*times))
  S <- survival_censored(y0, mu, sigma2, times)
  log_S <- (log_S1 + log_S2 + log_S3)/S
  total <- ifelse(delta, log_f, log_S)
  #total <- exp(y0) * total
  return(total)
}
