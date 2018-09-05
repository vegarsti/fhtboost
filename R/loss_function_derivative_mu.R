#' Loss function \code{rho} differentiated wrt. \code{mu}
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

loss_function_derivative_mu <- function(y0, mu, sigma2, t, delta) {
  # t is vector of times, delta is vector of observed = 1, or not = 0
  log_f <- - (y0 + mu*t)/sigma2
  log_S1 <- (t/sqrt(sigma2*t)) * dnorm((mu*t + y0)/(sqrt(sigma2*t)))
  log_S2 <- (2*y0/sigma2) * exp(-2*y0*mu/sigma2) * pnorm((mu*t-y0)/sqrt(sigma2*t))
  log_S3 <- -t/sqrt(sigma2*t) * exp(-2*y0*mu/sigma2)*dnorm((mu*t - y0)/sqrt(sigma2*t))
  log_S <- log_S1 + log_S2 + log_S3
  total <- delta*log_f + (1-delta)*log_S
  return(-sum(total))
}
