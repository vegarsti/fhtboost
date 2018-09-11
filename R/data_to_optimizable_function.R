#' Function taking design matrices \code{X}, \code{Z} and data \code{t}, {delta}, and returns a function to be optimized in nlm.
#'
#' Function taking design matrices \code{X}, \code{Z} and data \code{t}, {delta}, and returns a function to be optimized in nlm.
#'
#' @param X Design matrix for \code{ln y0 = beta^T * X}
#' @param Z Design matrix for \code{mu = gamma^T * Z}
#' @param t Survival times (possibly censored)
#' @param delta Vector of same length as \code{t}, indicating if observations are actually observed (1) or censored (0)
#'
#' @return \code{optimizable_function} Function taking parameters as a vector (containing \code{beta}, \code{gamma}), to be optimized over.
#'
#' @keywords keywords
#'
#' @export
#'
#' @examples
#' R code here showing how your function works

data_to_optimizable_function <- function(X, Z, t, delta) {
  link_function <- FHT_link_function(X, Z)
  FHT_loglikelihood <- FHT_loglikelihood_with_data(t, delta)
  p <- dim(X)[2]
  d <- dim(Z)[2]
  optimizable_function <- function(optimization_parameters) {
    beta <- optimization_parameters[1:p]
    gamma <- optimization_parameters[p+1:p+d]
    parameters <- link_function(beta, gamma)
    y0 <- parameters$y0
    mu <- parameters$mu
    loglikelihood_value <- FHT_loglikelihood(y0, mu)
    return(-loglikelihood_value)
  }
  return(optimizable_function)
}
