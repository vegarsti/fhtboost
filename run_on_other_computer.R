library(devtools)
library(foreach)
library(readr)
library(doParallel)
install_github("vegarsti/fhtboost")
library(fhtboost)

#seeds <- c(4, 10, 13, 14, 32, 36, 69, 71, 76, 80, 88, 91, 98, 104, 107, 139, 151, 164, 179, 187, 214, 230, 236, 241, 256, 259, 261, 280, 284, 292, 296, 297, 298, 302, 306, 308, 309, 313, 323, 333, 347, 372, 388, 392, 403, 410, 412, 416, 418, 420, 422, 474, 495, 496)
seeds <- 1:150
no_cores <- detectCores() - 1
registerDoParallel(cores=no_cores)
cl <- makeCluster(no_cores)
directory <- "./correlated-no-intercept/"
N <- 500
N_test <- 1000
setup_type <- 'correlated'
add_noise <- FALSE
K <- 5
K_fold_repetitions <- 5
M <- 150
TEST_SEED <- 9000
boost_intercepts_continually <- FALSE

foreach(seed=seeds) %dopar% {
  run_CV_and_write_to_file(N, setup_type, add_noise, seed, K, K_fold_repetitions, directory, M, boost_intercepts_continually)
  estimate_model_and_validate_and_write_to_file(N, N_test, setup_type, add_noise, seed, directory, boost_intercepts_continually)
}
stopCluster(cl)

