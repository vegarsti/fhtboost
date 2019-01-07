best_least_squares_update <- function(design_matrix, negative_gradient, number_of_predictors, predictors_to_consider, run_in_parallel=FALSE) { #, is_y0, y0) {
  parameter_estimates <- rep(NA, number_of_predictors)
  LARGE_NUMBER <- 90000000
  rss <- rep(LARGE_NUMBER, number_of_predictors)
  N <- dim(design_matrix)[1]
  negative_gradient_estimates <- matrix(NA, nrow=number_of_predictors, ncol=N)
  #print(sum(abs(negative_gradient)))
  # output_matrix <- matrix(negative_gradient, nrow=10)
  # print(output_matrix)
  # pause(20)

  # DO THIS IN PARALLEL
  if (run_in_parallel) {
    no_cores <- detectCores() - 1
    cluster <- makeCluster(no_cores, type="FORK")
    registerDoParallel(cluster)
    result <- foreach(j = predictors_to_consider, .combine = rbind)  %dopar% {
      design_matrix_j <- design_matrix[,j]
      parameter_estimate_j <- solve(t(design_matrix_j) %*% design_matrix_j) %*% t(design_matrix_j) %*% negative_gradient
      negative_gradient_estimates[j, ] <- design_matrix_j * parameter_estimate_j
      rss_j <- sum((negative_gradient_estimates[j, ] - negative_gradient)^2) / N
      return(c(parameter_estimate_j, rss_j))
    }
    stopCluster(cluster)
    rss <- result[, 2]
    parameter_estimates <- result[, 1]
    best_rss <- min(rss)
    best_predictor <- predictors_to_consider[which.min(rss)]
    parameter_estimates_best_predictor <- parameter_estimates[best_predictor]
  } else {
    for (j in predictors_to_consider) { # rewrite this to an apply call?
      design_matrix_j <- design_matrix[,j]
      parameter_estimates[j] <- solve(t(design_matrix_j) %*% design_matrix_j) %*% t(design_matrix_j) %*% negative_gradient
      negative_gradient_estimates[j, ] <- design_matrix_j * parameter_estimates[j]
      rss[j] <- sum((negative_gradient_estimates[j, ] - negative_gradient)^2) / N
    }
    best_predictor <- which.min(rss)
    best_rss <- min(rss)
    parameter_estimates_best_predictor <- parameter_estimates[best_predictor]
  }
  parameter_updates <- rep(0, number_of_predictors)
  parameter_updates[best_predictor] <- parameter_estimates_best_predictor
  return(list(parameter_updates=parameter_updates, rss=best_rss))
}
