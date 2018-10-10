create_folds <- function(N, K) {
  indices <- sample(1:N)
  folds <- matrix(indices, nrow = K)
  return(folds)
}
