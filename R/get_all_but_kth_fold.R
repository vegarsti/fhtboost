#' @export

get_all_but_kth_fold <- function(folds, k, K) {
  result <- c()
  for (i in 1:K) {
    if (i != k) {
      result <- c(result, folds[[i]])
    }
  }
  return(result)
}
