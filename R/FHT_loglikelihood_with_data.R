#' Get function for optimization
#'
#' Given observations (\code{t}, \code{delta}), return
#' a function evaluating the log likelihood for the FHT model
#'
#' @param t Survival times (possibly censored)
#' @param delta Corresponding observation indicator (1 = observed, 0 = censored)
#'
#' @return \code{FHT_loglikelihood} A function taking parameters \code{y0}, \code{mu}
#'
#' @keywords keywords
#'
#' @export
#'
#' @examples
#' Code examples.

FHT_loglikelihood_with_data <- function(t, delta) {
  FHT_loglikelihood_function <- function(y0, mu) {
    log_f <- log(y0) - 0.5*log(2*pi*t^3) - (y0 + mu*t)^2/(2*t)
    log_S <- log(pnorm((mu*t + y0)/sqrt(t)) - exp(-2*y0*mu)*pnorm((mu*t - y0))/sqrt(t))
    loglikelihood_value <- sum(delta*log_f + (1-delta)*log_S)
    return(loglikelihood_value)
  }
  return(FHT_loglikelihood_function)
}
