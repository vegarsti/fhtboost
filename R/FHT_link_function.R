#' Function taking design matrices \code{X}, \code{Z} and returns link function taking \code{beta}, \code{gamma}
#'
#' Function taking design matrices \code{X}, \code{Z} and returns link function taking \code{beta}, \code{gamma}
#'
#' @param X Design matrix for \code{ln y0 = beta^T * X}
#' @param Z Design matrix for \code{mu = gamma^T * Z}
#'
#' @return \code{link_function} Function taking \code{beta}, \code{gamma} and returning \code{y0}, \code{mu}.
#'
#' @keywords keywords
#'
#' @export
#'
#' @examples
#' R code here showing how your function works

FHT_link_function <- function(X, Z) {
  link_function <- function(beta, gamma) {
    y0 <- exp(X %*% beta)
    mu <- Z %*% gamma
    return(list(y0=y0, mu=mu))
  }
  return(link_function)
}
