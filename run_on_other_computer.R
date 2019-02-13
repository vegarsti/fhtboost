library(devtools)
library(foreach)
library(readr)
library(doParallel)
library(fhtboost)

# Configuration
#seeds <- 167:200
seeds <- c(1, 13, 26, 27, 41, 45, 49, 85, 88, 94, 99, 121, 150, 159, 162, 165, 166, 172, 174, 195, 196, 208, 212, 224, 227, 240, 253, 260, 261, 263, 269, 271, 289, 291, 292, 293, 303, 310, 320, 322, 357, 367, 371, 375, 379, 399, 402, 405, 406, 407, 420, 433, 443, 445, 449, 457, 466, 495)

no_cores <- detectCores() - 1
registerDoParallel(cores=no_cores)
cl <- makeCluster(no_cores)
directory <- "./simulations-no-intercept/"
N <- 500
N_test <- 1000
setup_type <- 'huge_clinical'
add_noise <- FALSE
K <- 5
K_fold_repetitions <- 5
M <- 130
TEST_SEED <- 9000
boost_intercepts_continually <- FALSE

foreach(seed=seeds) %dopar% {
  run_CV_and_write_to_file(N, setup_type, add_noise, seed, K, K_fold_repetitions, directory, M, boost_intercepts_continually)
  estimate_model_and_validate_and_write_to_file(N, setup_type, add_noise, seed, directory, boost_intercepts_continually)
}
stopCluster(cl)

