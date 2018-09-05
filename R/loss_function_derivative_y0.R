#' Loss function \code{rho} differentiated wrt. \code{y0}
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

rho_y0 <- function(y0, mu, sigma2, t, delta) {
  # t is vector of times, delta is vector of observed = 1, or not = 0
  log_f <- 1/y0 - (y0 + mu*t)/sigma2
  log_S1 <- (1/sqrt(sigma2*t)) * dnorm((mu*t + y0)/(sqrt(sigma2*t)))
  log_S2 <- (2*mu/sigma2) * exp(-2*y0*mu/sigma2) * pnorm((mu*t-y0)/sqrt(sigma2*t))
  log_S3 <- 1/sqrt(sigma2*t) * exp(-2*y0*mu/sigma2)*dnorm((mu*t - y0)/sqrt(sigma2*t))
  log_S <- log_S1 + log_S2 + log_S3
  total <- delta*log_f + (1-delta)*log_S
  return(-sum(total))
}
