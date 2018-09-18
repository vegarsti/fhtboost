#' Parameter vector to list
#'
#' Takes a parameter vector (from \code{nlm}) and splits it into correct parameter vectors.
#'
#' @param vec The vector in question
#' @param p Second dimension of design matrix \code{X}, i.e., size of \code{beta}.
#' @param d Second dimension of design matrix \code{Z}, i.e., size of \code{gamma}.
#'
#' @return \code{beta} Parameter vector
#' @return \code{gamma} Parameter vector
#'
#' @keywords keywords
#'
#' @export
#'
#' @examples
#' R code here showing how your function works

parameter_vector_to_list <- function(vec, p, d) {
  beta_vec <- vec[1:p]
  gamma_vec <- vec[(p+1):(p+d)]
  return(list(beta=beta_vec, gamma=gamma_vec))
}
