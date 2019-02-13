rm(list=ls())
library(devtools)
load_all()

# Configuration
seeds <- 1:500
not_seeds <- c(15, 17, 24, 30, 57, 58, 63, 68, 70, 78, 79, 90, 93, 95, 103, 110, 113, 123, 140, 141, 144, 149, 160, 170, 171, 179, 181, 183, 184, 189, 193, 200, 201, 206, 213, 217, 218, 236, 238, 241, 244, 247, 248, 249, 250, 251, 256, 276, 278, 308, 314, 316, 317, 325, 337, 338, 344, 349, 361, 366, 376, 381, 382, 391, 394, 396, 404, 421, 427, 428, 435, 459, 461, 463, 464, 467, 472, 490, 491, 492, 494, 500)
seeds <- seeds[-not_seeds]

B <- length(seeds)

base_input_directory <- "../dataset/"
algorithm <- 'non-cyclic-intercept'
input_directory <- paste(base_input_directory, algorithm, sep='')
output_directory <- '../simulation-results/'
scenario <- 'uncorrelated'

informative_p <- 35
p <- 10000
non_informative_p <- p - informative_p
informative_d <- 5
d <- 15
non_informative_d <- d - informative_d



### READING RESULTS

# initialize
initial_seed <- seeds[1]
filename <- make_filename(input_directory, "beta", initial_seed)
beta <- read.csv(filename)
total_beta <- beta$x

filename <- make_filename(input_directory, "gamma", initial_seed)
gamma <- read.csv(filename)
total_gamma <- gamma$x

filename <- make_filename(input_directory, "summary", initial_seed)
summary <- read.csv(filename)
total_summary <- summary

# run on datasets
for (seed in seeds[-1]) {
  filename <- make_filename(input_directory, "beta", seed)
  beta <- read.csv(filename)
  filename <- make_filename(input_directory, "gamma", seed)
  gamma <- read.csv(filename)
  filename <- make_filename(input_directory, "summary", seed)
  summary <- read.csv(filename)
  total_beta <- cbind(total_beta, beta$x)
  total_gamma <- cbind(total_gamma, gamma$x)
  total_summary <- rbind(total_summary, summary)
}
total_beta <- total_beta[-1, ]
total_gamma <- total_gamma[-1, ]



### PROCESSING RESULTS ####
## BETA
summarize_variable_selection <- function(all_parameter_results, p, informative_p) {
  non_informative_p <- p - informative_p
  result <- t(apply(all_parameter_results, 2, function(column) {
    non_null <- column != 0
    chosen_parameters <- sum(non_null)
    informative_chosen <- sum(non_null[1:informative_p])
    non_informative_chosen <- sum(non_null[(informative_p+1):p])
    c(chosen_parameters, informative_chosen, non_informative_chosen)
  }))
  P <- informative_p
  N <- non_informative_p
  result <- data.frame(TP=result[, 2], FP=result[, 3])
  result$FN <- informative_p - result$TP
  result$TN <- non_informative_p - result$FP
  result$sensitivity <- result$TP/P # TPR
  result$specificity <- result$TN/N # TNR
  result$accuracy <- (result$TN+result$TP)/(N+P) # TNR
  result$FPR <- result$FP/N # false positive rate
  result$FNR <- result$FN/P # false negative rate
  #result_beta$FDR <- result_beta$FP/(result_beta$FP+result_beta$TP) # false discovery rate
  result$FOR <- result$FN/(result$FN+result$TN) # false omission rate
  result$NPV <- result$TN/(result$FN+result$TN) # false omission rate

  result_mean <- as.data.frame(colMeans(result))
  result_sd <- as.data.frame(apply(result, 2, sd))
  result_summary <- as.data.frame(cbind(result_mean, result_sd))
  colnames(result_summary) <- c("mean", "sd")
  result_summary$rownames <- c("TP", "FP", "FN", "TN", "sensitivity", "specificity", "accuracy", "FPR", "FNR", "FOR", "NPV")
  result_summary <- result_summary[-c(1:4), ]
  return(result_summary)
}

result_beta_summary <- summarize_variable_selection(total_beta, p, informative_p)
filename <- paste(algorithm, scenario, 'beta', sep='_')
full_filename <- paste(output_directory, filename, '.csv', sep='')
write.csv(result_beta_summary[, -3], full_filename, row.names=result_beta_summary$rownames)


# GAMMA
result_gamma_summary <- summarize_variable_selection(total_gamma, d, informative_d)
filename <- paste(algorithm, scenario, 'gamma', sep='_')
full_filename <- paste(output_directory, filename, '.csv', sep='')
write.csv(result_gamma_summary[, -3], full_filename, row.names=result_gamma_summary$rownames)


## ALL (results from each run)
summary_means <- colMeans(total_summary)
summary_sd <- apply(total_summary, 2, sd)
summary_min <- apply(total_summary, 2, min)
summary_max <- apply(total_summary, 2, max)
summary <- cbind(summary_means, summary_sd, summary_min, summary_max)
colnames(summary) <- c("mean", "sd", "min", "max")
summary <- as.data.frame(summary)

filename <- paste(algorithm, scenario, 'summary', sep='_')
full_filename <- paste(output_directory, filename, '.csv', sep='')
write.csv(summary, full_filename)
