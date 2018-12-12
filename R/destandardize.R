destandardize <- function(matrix_to_scale, column_wise_weights, means) {
  p <- dim(matrix_to_scale)[2]
  if (is.null(dim(matrix_to_scale))) {
    matrix_to_scale <- matrix_to_scale + means*column_wise_weights
  }
  else {
    if (length(column_wise_weights) != p) {
      stop('Dimensions of matrix and vector of weights do not match!')
    }
    if (length(means) != p) {
      stop('Dimensions of matrix and vector of means do not match!')
    }
    for (j in 1:p) {
      matrix_to_scale[, j] <- (matrix_to_scale[, j] + means[j])*column_wise_weights[j]
    }
  }
  return(matrix_to_scale)
}
