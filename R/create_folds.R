create_folds <- function(to_divide, K) {
  N <- length(to_divide)
  indices <- sample(to_divide)
  whole <- N %/% K
  rest <- N %% K
  folds <- list()
  end <- 0
  for (k in 1:K) {
    start <- end + 1
    entries <- whole + ifelse(k <= rest, 1, 0)
    end <- start + (entries-1)
    folds[[k]] <- indices[start:end]
  }
  return(folds)
}
