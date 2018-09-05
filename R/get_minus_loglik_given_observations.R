#' Get function for optimization
#'
#' Given observations (\code{t}, \code{delta}), return
#' a function to be optimized, which evaluates the loss function
#' (the negative log likelihood)
#'
#' @param t Survival times (possibly censored)
#' @param delta Corresponding observation indicator (1 = observed, 0 = censored)
#'
#' @return \code{minusloglik} A function taking one argument -- a vector \code{parameters} -- which is to be optimized.
#'
#' @keywords keywords
#'
#' @export
#'
#' @examples
#' t <- c(0.1, 1.0)
#' delta <- c(0, 1)
#' minus_loglik <- get_minus_loglik(t, delta)
#' initial_parameters <- c(1, 1)
#' optimized_parameters <- nlm(minus_loglik, initial_parameters)$estimate

get_minus_loglikelihood_given_observations <- function(t, delta) {
  FHT_loglik <- function(y0, mu) {
    log_f <- log(y0) - 0.5*log(2*pi*t^3) - (y0 + mu*t)^2/(2*t)
    log_S <- log(pnorm((mu*t + y0)/sqrt(t)) - exp(-2*y0*mu)*pnorm((mu*t - y0))/sqrt(t))
    sum(delta*log_f + (1-delta)*log_S)
  }
  minusFHT_loglik <- function(parameters) {
    y0 <- parameters[1]
    mu <- parameters[2]
    -FHT_loglik(y0, mu)
  }
  return(minusFHT_loglik)
}
