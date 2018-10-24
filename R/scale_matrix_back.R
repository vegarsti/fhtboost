scale_matrix_back <- function(matrix_to_scale, column_wise_weights) {
  t(t(matrix_to_scale) * column_wise_weights)
}
