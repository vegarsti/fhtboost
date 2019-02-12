rm(list=ls())
library(devtools)
load_all()

# Configuration
seeds <- 1:500
#not_seeds <- c(15, 17, 24, 30, 57, 58, 63, 68, 70, 78, 79, 90, 93, 95, 103, 110, 113, 123, 140, 141, 144, 149, 160)
not_seeds <- c(15, 17, 24, 30, 57, 58, 63, 68, 70, 78, 79, 90, 93, 95, 103, 110, 113, 123, 140, 141, 144, 149, 160, 170, 171, 179, 181, 183, 184, 189, 193, 200, 201, 206, 213, 217, 218, 236, 248, 249, 251, 256, 276, 278, 308, 314, 316, 317, 337, 338, 344, 349, 361, 366, 376, 381, 382, 391, 394, 396, 404, 421, 427, 428, 435, 459, 461, 463, 464, 467, 472, 490, 491, 492, 494, 500)
seeds <- seeds[-not_seeds]

B <- length(seeds)

directory <- "../dataset/"

# initialize
initial_seed <- seeds[1]
filename <- make_filename(directory, "beta", initial_seed)
beta <- read.csv(filename)
total_beta <- beta$x

filename <- make_filename(directory, "gamma", initial_seed)
gamma <- read.csv(filename)
total_gamma <- gamma$x

filename <- make_filename(directory, "summary", initial_seed)
summary <- read.csv(filename)
total_summary <- summary


# run on datasets
for (seed in seeds[-1]) {
  filename <- make_filename(directory, "beta", seed)
  beta <- read.csv(filename)
  filename <- make_filename(directory, "gamma", seed)
  gamma <- read.csv(filename)
  filename <- make_filename(directory, "summary", seed)
  summary <- read.csv(filename)
  total_beta <- cbind(total_beta, beta$x)
  total_gamma <- cbind(total_gamma, gamma$x)
  total_summary <- rbind(total_summary, summary)
}
total_beta <- total_beta[-1, ]
total_gamma <- total_gamma[-1, ]


# BETA
informative_p <- 35
p <- 10000
non_informative_p <- p - informative_p


result_beta <- t(apply(total_beta, 2, function(column) {
  non_null <- column != 0
  chosen_parameters <- sum(non_null)
  informative_chosen <- sum(non_null[1:informative_p])
  non_informative_chosen <- sum(non_null[(informative_p+1):p])
  c(chosen_parameters, informative_chosen, non_informative_chosen)
}))
result_beta <- data.frame("n"=result_beta[, 1], "informative"=result_beta[, 2], "non_informative"=result_beta[, 3])
result_beta$"informative_ratio" <- result_beta$informative / result_beta$n
result_beta$"non_informative_ratio" <- result_beta$non_informative / result_beta$n



non_null_parameter_count <- apply(total_beta, 1, function(row) length(non_null_parameters(row))) / B
number_of_chosen_parameters <- sum(non_null_parameter_count)

# parameters with signal
signal_parameters <- non_null_parameter_count[1:informative_p]
# sum(signal_parameters) / informative_p
non_signal_parameters <- non_null_parameter_count[(informative_p+1):p]
# sum(non_signal_parameters) / non_informative_p

TP_rate <- mean(signal_parameters)
FN_rate <- 1 - TP_rate
FP_rate <- mean(non_signal_parameters) # mean?
TN_rate <- 1 - FP_rate

confusion_matrix <- matrix(c(TP_rate, FN_rate, FP_rate, TN_rate), nrow=2)


# GAMMA
informative_d <- 5
d <- 15
non_informative_d <- d - informative_d



result_gamma <- t(apply(total_gamma, 2, function(column) {
  non_null <- column != 0
  chosen_parameters <- sum(non_null)
  informative_chosen <- sum(non_null[1:informative_d])
  non_informative_chosen <- sum(non_null[(informative_d+1):d])
  c(chosen_parameters, informative_chosen, non_informative_chosen)
}))
result_gamma <- data.frame("n"=result_gamma[, 1], "informative"=result_gamma[, 2], "non_informative"=result_gamma[, 3])
result_gamma$"informative_ratio" <- result_gamma$informative / result_gamma$n
result_gamma$"non_informative_ratio" <- result_gamma$non_informative / result_gamma$n
