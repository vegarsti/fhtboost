rm(list=ls())

library(devtools)
library(foreach)
library(readr)
library(doParallel)
load_all() # or library(fhtboost)

# Configuration
number_of_data_sets <- 100

#seeds <- 1:number_of_data_sets
seeds <- 161:180

seed <- 433

no_cores <- detectCores() - 1
registerDoParallel(cores=no_cores)
cl <- makeCluster(no_cores)
directory <- "../dataset/non-cyclic-no-intercept/"
N <- 500
N_test <- 1000
setup_type <- 'huge_clinical'
add_noise <- FALSE
K <- 5
K_fold_repetitions <- 5
M <- 200
TEST_SEED <- 9000
boost_intercepts_continually <- FALSE

foreach(seed=seeds) %dopar% {
  run_CV_and_write_to_file(N, setup_type, add_noise, seed, K, K_fold_repetitions, directory, M, boost_intercepts_continually)
  estimate_model_and_validate_and_write_to_file(N, setup_type, add_noise, seed, directory, boost_intercepts_continually)
}
stopCluster(cl)
