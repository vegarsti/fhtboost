create_folds_stratified <- function(delta, K) {
  observed <- which(delta == 1)
  censored <- which(delta == 0)
  observed_folds <- create_folds(observed, K)
  censored_folds <- create_folds(censored, K)
  folds <- list()
  for (i in 1:K) {
    folds[[i]] <- c(observed_folds[[i]], censored_folds[[i]])
  }
  return(folds)
}
