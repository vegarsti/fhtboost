#' Get function for optimization
#'
#' Given observations (\code{t}, \code{delta}), return
#' a function evaluating the log likelihood for the FHT model
#'
#' @param times Survival times (possibly censored)
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

FHT_loglikelihood_with_data <- function(times, delta) {
  FHT_loglikelihood_function <- function(y0, mu) {
    log_f <- log(y0) - 0.5*log(2*pi*times^3) - (y0 + mu*times)^2/(2*times)
    log_S <- log(pnorm((mu*times + y0)/sqrt(times)) - exp(-2*y0*mu)*pnorm((mu*times - y0))/sqrt(times))
    loglikelihood_value <- sum(delta*log_f) + sum((1-delta)*log_S)
    return(loglikelihood_value)
  }
  return(FHT_loglikelihood_function)
}
