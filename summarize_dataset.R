rm(list=ls())
library(devtools)
load_all()

# Configuration
seeds <- 1:500
not_seeds <- c(15, 17, 18, 24, 30, 57, 58, 63, 68, 70, 78, 79, 90, 93, 95, 103, 110, 113, 123, 140, 141, 144, 149, 160, 170, 171, 179, 181, 183, 184, 189, 193, 200, 201, 206, 213, 217, 218, 236, 238, 241, 244, 247, 248, 249, 250, 251, 256, 276, 278, 308, 314, 316, 317, 325, 337, 338, 344, 349, 361, 366, 376, 381, 382, 391, 394, 396, 404, 421, 427, 428, 435, 459, 461, 463, 464, 467, 472, 490, 491, 492, 493, 494, 500)
seeds <- seeds[-not_seeds]

B <- length(seeds)

base_input_directory <- "../dataset/scenario-"
algorithm <- 'non-cyclic-no-intercept'
algorithm_string <- "without intercept"
scenario <- 'non-correlated'
input_directory <- paste0(base_input_directory, scenario, '/', algorithm, '/')
output_directory <- paste0('../dataset/scenario-', scenario, '/', algorithm, '-results/')

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
  result <- tryCatch({
    filename <- make_filename(input_directory, "beta", seed)
    beta <- read.csv(filename)
    filename <- make_filename(input_directory, "gamma", seed)
    gamma <- read.csv(filename)
    filename <- make_filename(input_directory, "summary", seed)
    summary <- read.csv(filename)
    total_beta <- cbind(total_beta, beta$x)
    total_gamma <- cbind(total_gamma, gamma$x)
    total_summary <- rbind(total_summary, summary)
  }, warning = function(w) {}, error = function(e) {
    cat(seed, '\n')
  }, finally = { })
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
  result$FDR <- result$FP/(result$FP+result$TP) # false discovery rate
  result$FOR <- result$FN/(result$FN+result$TN) # false omission rate
  result$NPV <- result$TN/(result$FN+result$TN) # false omission rate

  nan_indices <- which(is.nan(result$FDR))
  if (length(nan_indices) != 0) {
    result <- result[-nan_indices,]
  }

  result_mean <- as.data.frame(colMeans(result))
  result_median <- as.data.frame(apply(result, 2, median))
  result_sd <- as.data.frame(apply(result, 2, sd))
  result_summary <- as.data.frame(cbind(result_mean, result_median, result_sd))
  colnames(result_summary) <- c("mean", "median", "sd")
  result_summary$rownames <- c("TP", "FP", "FN", "TN", "sensitivity", "specificity", "accuracy", "FPR", "FNR", "FDR", "FOR", "NPV")
  result_summary <- result_summary[-c(1:4), ]
  return(list(summary=result_summary, whole=result))
}

result_beta_summary <- summarize_variable_selection(total_beta, p, informative_p)
result_beta_whole <- result_beta_summary$whole
result_beta_summary <- result_beta_summary$summary
filename <- paste(algorithm, scenario, 'beta', sep='_')
full_filename <- paste(output_directory, filename, '.csv', sep='')
write.csv(result_beta_summary[, -4], full_filename, row.names=result_beta_summary$rownames)
filename <- paste(algorithm, scenario, 'beta_whole', sep='_')
full_filename <- paste(output_directory, filename, '.csv', sep='')
write.csv(result_beta_whole, full_filename)


# GAMMA
result_gamma_summary <- summarize_variable_selection(total_gamma, d, informative_d)
result_gamma_whole <- result_gamma_summary$whole
result_gamma_summary <- result_gamma_summary$summary
filename <- paste(algorithm, scenario, 'gamma', sep='_')
full_filename <- paste(output_directory, filename, '.csv', sep='')
write.csv(result_gamma_summary[, -4], full_filename, row.names=result_gamma_summary$rownames)
filename <- paste(algorithm, scenario, 'gamma_whole', sep='_')
full_filename <- paste(output_directory, filename, '.csv', sep='')
write.csv(result_gamma_whole, full_filename)

filename <- "beta_variable_selection_boxplot.pdf"
full_filename <- paste0(output_directory, filename)
pdf(full_filename, width=12, height=6)
par(oma=c(0,1,0,0))
boxplot(as.numeric(result_beta_whole$FDR), as.numeric(result_beta_whole$sensitivity),
        as.numeric(result_beta_whole$specificity), horizontal=TRUE, axes=FALSE,
        main=paste0("Variable selection metrics for beta (uncorrelated scenario, ", algorithm_string, ")"))
axis(1)
axis(2, labels=c("FDR", "Sensitivity", "Specificity"), at=1:3, las=2)
box()
par(oma=c(0,0,0,0))
dev.off()


filename <- "gamma_variable_selection_boxplot.pdf"
full_filename <- paste0(output_directory, filename)
pdf(full_filename, width=12, height=6)
par(oma=c(0,1,0,0))
boxplot(as.numeric(result_gamma_whole$FDR), as.numeric(result_gamma_whole$sensitivity),
        as.numeric(result_gamma_whole$specificity), horizontal=TRUE, axes=FALSE,
        main=paste0("Variable selection metrics for gamma (uncorrelated scenario, ", algorithm_string, ")"))
axis(1)
axis(2, labels=c("FDR", "Sensitivity", "Specificity"), at=1:3, las=2)
box()
par(oma=c(0,0,0,0))
dev.off()






## ALL (results from each run)
summary_means <- colMeans(total_summary)
summary_median <- apply(total_summary, 2, median)
summary_sd <- apply(total_summary, 2, sd)
summary_min <- apply(total_summary, 2, min)
summary_max <- apply(total_summary, 2, max)
summary <- cbind(summary_means, summary_sd, summary_min, summary_max, summary_median)
colnames(summary) <- c("mean", "sd", "min", "max", "median")
summary <- as.data.frame(summary)

filename <- paste(algorithm, scenario, 'summary', sep='_')
full_filename <- paste(output_directory, filename, '.csv', sep='')
write.csv(summary, full_filename)

filename <- paste(algorithm, scenario, 'total_summary', sep='_')
full_filename <- paste(output_directory, filename, '.csv', sep='')
write.csv(total_summary, full_filename)

loglik_ylim <- c(min(min(total_summary$null_loglik), min(total_summary$loglik)), max(max(total_summary$null_loglik), max(total_summary$loglik)))

# Plot log likelihood and null likelihoods
plot(total_summary$loglik, ylim=loglik_ylim)
points(total_summary$null_loglik, pch='+', col='red')

# Plot deviance
plot(total_summary$deviance, ylim=c(min(total_summary$deviance), max(total_summary$deviance)))

deviance_intercept <- total_summary$deviance
deviance_no_intercept <- total_summary$deviance


directory <- paste0(base_input_directory, scenario, '/')
full_filename <- paste0(directory, 'deviance_both_boxplot_no_title.pdf')
pdf(full_filename, width=12, height=6)
boxplot(-deviance_intercept, -deviance_no_intercept, xlab='Difference in deviance', horizontal=TRUE,
        main="Difference in deviance on non-correlated data, with/without boosting intercept")
axis(2, labels=c("With intercept", "Without intercept"), at=1:2, las=2)
abline(v=0, lty=3)
dev.off()

