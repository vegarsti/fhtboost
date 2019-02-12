rm(list=ls())
library(devtools)
load_all()

# Configuration
seeds <- 1:500
not_seeds <- c(15, 17, 24, 30, 57, 58, 63, 68, 70, 78, 79, 90, 93, 95, 103, 110, 113, 123, 140, 141, 144, 149, 160, 170, 171, 179, 181, 183, 184, 189, 193, 200, 201, 206, 213, 217, 218, 236, 248, 249, 251, 256, 276, 278, 308, 314, 316, 317, 337, 338, 344, 349, 361, 366, 376, 381, 382, 391, 394, 396, 404, 421, 427, 428, 435, 459, 461, 463, 464, 467, 472, 490, 491, 492, 494, 500)
seeds <- seeds[-not_seeds]

B <- length(seeds)

input_directory <- "../dataset/"
output_directory <- '../simulation-results/'
algorithm <- 'non-cyclic'
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
result_beta <- t(apply(total_beta, 2, function(column) {
  non_null <- column != 0
  chosen_parameters <- sum(non_null)
  informative_chosen <- sum(non_null[1:informative_p])
  non_informative_chosen <- sum(non_null[(informative_p+1):p])
  c(chosen_parameters, informative_chosen, non_informative_chosen)
}))
result_beta <- data.frame(beta_n=result_beta[, 1], beta_informative=result_beta[, 2], beta_non_informative=result_beta[, 3])

non_null_parameter_count_beta <- apply(total_beta, 1, function(row) length(non_null_parameters(row)))# / B
number_of_chosen_parameters_beta <- sum(non_null_parameter_count_beta)
total_number_of_parameters_beta <- p * B
total_number_of_informative_parameters_beta <- informative_p * B
total_number_of_non_informative_parameters_beta <- non_informative_p * B


informative_parameters <- non_null_parameter_count_beta[1:informative_p]
non_informative_parameters <- non_null_parameter_count_beta[(informative_p+1):p]

number_of_chosen_informative_parameters_beta <- sum(informative_parameters)
number_of_chosen_non_informative_parameters_beta <- sum(non_informative_parameters)

N <- total_number_of_non_informative_parameters_beta
P <- total_number_of_informative_parameters_beta
TP <- number_of_chosen_informative_parameters_beta
FN <- total_number_of_informative_parameters_beta - number_of_chosen_informative_parameters_beta
FP <- number_of_chosen_non_informative_parameters_beta
TN <- total_number_of_non_informative_parameters_beta - number_of_chosen_non_informative_parameters_beta

FPR <- FP/N # false positive rate
FNR <- FN/P # false negative rate
FDR <- FP/(FP+TP)
FOR <- FN/(FN+TN)
precision <- PPV <- TP/(TP+FP) # positive predictive value
sensitivity <- TPR <- TP/P # true positive rate
specificity <- TNR <- TN/N # true negative rate
accuracy <- (TP+TN)/(P+N)
confusion_matrix_beta <- matrix(c(TP, FN, FP, TN), nrow=2)

beta_results_df <- data.frame(
  N=N,P=P,TP=TP,FN=FN,FP=FP,TN=TN,
  false_positive_rate=FPR,false_negative_rate=FNR,false_discovery_rate=FDR,false_omission_rate=FOR,
  precision=precision,sensitivity=sensitivity,specificity=specificity,accuracy=accuracy
)


filename <- paste(algorithm, scenario, 'beta', sep='_')
full_filename <- paste(output_directory, filename, '.csv', sep='')
write.csv(beta_results_df, full_filename)



# GAMMA

result_gamma <- t(apply(total_gamma, 2, function(column) {
  non_null <- column != 0
  chosen_parameters <- sum(non_null)
  informative_chosen <- sum(non_null[1:informative_d])
  non_informative_chosen <- sum(non_null[(informative_d+1):d])
  c(chosen_parameters, informative_chosen, non_informative_chosen)
}))
result_gamma <- data.frame(gamma_n=result_gamma[, 1], gamma_informative=result_gamma[, 2], gamma_non_informative=result_gamma[, 3])


non_null_parameter_count_gamma <- apply(total_gamma, 1, function(row) length(non_null_parameters(row)))
number_of_chosen_parameters_gamma <- sum(non_null_parameter_count_gamma)
total_number_of_parameters_gamma <- d * B
total_number_of_informative_parameters_gamma <- informative_d * B
total_number_of_non_informative_parameters_gamma <- non_informative_d * B


informative_parameters_gamma <- non_null_parameter_count_gamma[1:informative_d]
non_informative_parameters_gamma <- non_null_parameter_count_gamma[(informative_d+1):d]

number_of_chosen_informative_parameters_gamma <- sum(informative_parameters_gamma)
number_of_chosen_non_informative_parameters_gamma <- sum(non_informative_parameters_gamma)

N <- total_number_of_non_informative_parameters_gamma
P <- total_number_of_informative_parameters_gamma
TP <- number_of_chosen_informative_parameters_gamma
FN <- total_number_of_informative_parameters_gamma - number_of_chosen_informative_parameters_gamma
FP <- number_of_chosen_non_informative_parameters_gamma
TN <- total_number_of_non_informative_parameters_gamma - number_of_chosen_non_informative_parameters_gamma

FPR <- FP/N # false positive rate
FNR <- FN/P # false negative rate
FDR <- FP/(FP+TP)
FOR <- FN/(FN+TN)
precision <- PPV <- TP/(TP+FP) # positive predictive value
sensitivity <- TPR <- TP/P # true positive rate
specificity <- TNR <- TN/N # true negative rate
accuracy <- (TP+TN)/(P+N)
confusion_matrix_gamma <- matrix(c(TP, FN, FP, TN), nrow=2)

gamma_results_df <- data.frame(
  N=N,P=P,TP=TP,FN=FN,FP=FP,TN=TN,
  false_positive_rate=FPR,false_negative_rate=FNR,false_discovery_rate=FDR,false_omission_rate=FOR,
  precision=precision,sensitivity=sensitivity,specificity=specificity,accuracy=accuracy
)


filename <- paste(algorithm, scenario, 'gamma', sep='_')
full_filename <- paste(output_directory, filename, '.csv', sep='')
write.csv(gamma_results_df, full_filename)



## ALL (results from each run)
all_results <- cbind(result_beta, result_gamma, total_summary)
filename <- paste(algorithm, scenario, 'all', sep='_')
full_filename <- paste(output_directory, filename, '.csv', sep='')
write.csv(all_results, full_filename)
