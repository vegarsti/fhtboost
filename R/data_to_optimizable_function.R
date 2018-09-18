#' Data to optimizable function
#'
#' Takes design matrices \code{X}, \code{Z} and data \code{t}, {delta}, and returns a function to be optimized.
#'
#' @param X Design matrix for \code{ln y0 = beta^T * X}
#' @param Z Design matrix for \code{   mu = gamma^T * Z}
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

  minus_FHT_loglikelihood <- function(optimization_parameters) {
    parameter_list <- parameter_vector_to_list(optimization_parameters, p, d)
    FHT_parameters <- link_function(parameter_list$beta, parameter_list$gamma)
    minus_loglikelihood_value <- - FHT_loglikelihood(FHT_parameters$y0, FHT_parameters$mu)
    return(minus_loglikelihood_value)
  }

  return(minus_FHT_loglikelihood)
}
