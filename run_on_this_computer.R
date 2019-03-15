library(devtools)
library(foreach)
library(readr)
library(doParallel)
load_all()

seeds <- 1:500
# no_cores <- detectCores() - 1
# registerDoParallel(cores=no_cores)
# cl <- makeCluster(no_cores)

#### ONLY CHANGE HERE!!
scenario <- 'correlated'
boost_intercepts_continually <- FALSE
#### ONLY CHANGE HERE!!!

base_input_directory <- "../dataset/scenario-"
if (scenario == 'correlated') {
  not_seeds <- c(1, 3, 4, 6, 11, 23, 25, 26, 33, 35, 36, 37, 40, 41, 42, 45, 52, 54, 57, 61, 62, 63, 64, 65, 66, 68, 70, 72, 80, 81, 84, 87, 92, 95, 101, 107, 109, 110, 111, 114, 120, 121, 122, 127, 130, 133, 136, 149, 152, 155, 156, 158, 160, 168, 174, 180, 182, 193, 195, 196, 210, 211, 212, 217, 221, 223, 225, 233, 234, 242, 247, 252, 253, 257, 258, 259, 277, 278, 282, 285, 288, 293, 302, 313, 315, 316, 317, 319, 320, 334, 337, 338, 339, 340, 343, 347, 354, 355, 356, 357, 358, 361, 362, 365, 372, 381, 384, 389, 392, 393, 411, 412, 415, 421, 423, 428, 430, 433, 435, 439, 440, 442, 444, 445, 447, 456, 462, 466, 468, 472, 477, 479, 483, 485, 488, 493)
  scenario_directory <- paste0(base_input_directory, "correlated/")
} else if (scenario == 'huge_clinical') {
  not_seeds <- c(15, 17, 18, 24, 30, 57, 58, 63, 68, 70, 78, 79, 90, 93, 95, 103, 110, 113, 123, 140, 141, 144, 149, 160, 170, 171, 179, 181, 183, 184, 189, 193, 200, 201, 206, 213, 217, 218, 236, 238, 241, 244, 247, 248, 249, 250, 251, 256, 276, 278, 308, 314, 316, 317, 325, 337, 338, 344, 349, 361, 366, 376, 381, 382, 391, 394, 396, 404, 421, 427, 428, 435, 459, 461, 463, 464, 467, 472, 490, 491, 492, 493, 494, 500)
  scenario_directory <- paste0(base_input_directory, "non-correlated/")
} else {
  stop('Nonexisting scenario!')
}

seeds <- seeds[-not_seeds]
if (boost_intercepts_continually) {
  algorithm <- 'non-cyclic-intercept'
} else {
  algorithm <- 'non-cyclic-no-intercept'
}
input_directory <- paste0(scenario_directory, algorithm, '/')
directory <- input_directory

N <- 500
N_test <- 1000
setup_type <- scenario
add_noise <- FALSE
K <- 5
K_fold_repetitions <- 5
M <- 150
TEST_SEED <- 9000

#seeds <- run_again <- c(201, 203, 204, 206, 226, 235, 236, 237, 240, 241, 245, 254, 261, 262, 263, 264, 265, 266, 268, 270, 272, 280, 281, 284, 287, 292, 295, 301, 307, 309, 310, 311, 314, 321, 322, 327, 330, 333, 336, 349, 352, 360, 368, 374, 380, 382, 395, 396, 410, 417, 425, 434, 452, 453, 457, 458, 459, 478, 482)

foreach(seed=seeds) %dopar% {
  #run_CV_and_write_to_file(N, setup_type, add_noise, seed, K, K_fold_repetitions, directory, M, boost_intercepts_continually)
  estimate_model_and_validate_and_write_to_file(N, N_test, setup_type, add_noise, seed, directory, boost_intercepts_continually)
}
stopCluster(cl)

