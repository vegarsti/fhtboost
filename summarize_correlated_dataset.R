rm(list=ls())
library(devtools)
load_all()

# Configuration
seeds <- 1:500
not_seeds <- c(1, 3, 4, 6, 11, 23, 25, 26, 33, 35, 36, 37, 40, 41, 42, 45, 52, 54, 57, 61, 62, 63, 64, 65, 66, 68, 70, 72, 80, 81, 84, 87, 92, 95, 101, 107, 109, 110, 111, 114, 120, 121, 122, 127, 130, 133, 136, 149, 152, 155, 156, 158, 160, 168, 174, 180, 182, 193, 195, 196, 210, 211, 212, 217, 221, 223, 225, 233, 234, 242, 247, 252, 253, 257, 258, 259, 277, 278, 282, 285, 288, 293, 302, 313, 315, 316, 317, 319, 320, 334, 337, 338, 339, 340, 343, 347, 354, 355, 356, 357, 358, 361, 362, 365, 372, 381, 384, 389, 392, 393, 411, 412, 415, 421, 423, 428, 430, 433, 435, 439, 440, 442, 444, 445, 447, 456, 462, 466, 468, 472, 477, 479, 483, 485, 488, 493)
seeds <- seeds[-not_seeds]

B <- length(seeds)

base_input_directory <- "../dataset/scenario-correlated/"
algorithm <- 'non-cyclic-no-intercept'
scenario <- 'correlated'
input_directory <- paste0(base_input_directory, algorithm, '/')
output_directory <- paste0(base_input_directory, algorithm, '-results/')

informative_p <- 35
p <- 10000
non_informative_p <- p - informative_p
informative_d <- 10
d <- 25
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
summarize_variable_selection <- function(all_parameter_results, p, informative_p_vector) {
  informative_p <- length(informative_p_vector)
  non_informative_p <- p - informative_p
  result <- t(apply(all_parameter_results, 2, function(column) {
    non_null <- column != 0
    chosen_parameters <- sum(non_null)
    informative_chosen <- sum(non_null[informative_p_vector])
    non_informative_chosen <- sum(non_null[-informative_p_vector])
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

informative_p_vector <- c(2, 3, 4, 12, 13, 14, 15, 22, 23, 24, 32, 33, 34, 35, 42, 43, 44, 52,
    53, 54, 55, 62, 63, 64, 72, 73, 74, 75, 82, 83, 84, 92, 93, 94, 95) - 1 # to adjust for removal of the intercept

result_beta_summary <- summarize_variable_selection(total_beta, p, informative_p_vector)
filename <- paste(algorithm, scenario, 'beta', sep='_')
full_filename <- paste(output_directory, filename, '.csv', sep='')
write.csv(result_beta_summary[, -3], full_filename, row.names=result_beta_summary$rownames)


# GAMMA
informative_d_vector <- c(2, 4, 6, 8, 10, 12, 15, 18, 21, 24) - 1 # to adjust for the intercept
result_gamma_summary <- summarize_variable_selection(total_gamma, d, informative_d_vector)
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

plot_file <- function(seed) {
  seed_string <- formatC(seed, width=4, flag="0")
  filename <- paste0('../dataset/scenario-correlated/non-cyclic-no-intercept/', seed_string, '_cv_loglik.csv')
  df <- read.csv(filename)
  plot(rowMeans(df), typ='l', ylim=c(min(apply(df, 2, min)), max(apply(df, 2, max))))
  for (k in 1:5) { lines(df[, k], lty=3) }
  return(rowMeans(df))
}

lol <- FALSE
if (lol) {
  base_input_directory <- "../dataset/scenario-correlated/"
  algorithm_no <- 'non-cyclic-no-intercept'
  scenario <- 'correlated'
  input_directory_no <- paste0(base_input_directory, algorithm_no, '/')
  algorithm_yes <- 'non-cyclic-intercept'
  input_directory_yes <- paste0(base_input_directory, algorithm_yes, '/')

  errs_loglik_no <- c()
  errs_dev_no <- c()
  errs_loglik_yes <- c()
  errs_dev_yes <- c()
  for (seed in seeds) {
    seed_string <- formatC(seed, width=4, flag="0")

    # no intercept
    filename <- paste0(input_directory_no, seed_string, '_cv_loglik.csv')
    df <- read.csv(filename)
    errs_loglik_no <- cbind(rowMeans(df), errs_loglik_no)
    filename <- paste0(input_directory_no, seed_string, '_cv_deviance.csv')
    df <- read.csv(filename)
    errs_dev_no <- cbind(rowMeans(df), errs_dev_no)

    # intercept
    filename <- paste0(input_directory_yes, seed_string, '_cv_loglik.csv')
    df <- read.csv(filename)
    errs_loglik_yes <- cbind(rowMeans(df), errs_loglik_yes)
    filename <- paste0(input_directory_yes, seed_string, '_cv_deviance.csv')
    df <- read.csv(filename)
    errs_dev_yes <- cbind(rowMeans(df), errs_dev_yes)
  }
  #pdf(paste0(paste(scenario, "cv", "loglik", sep='_'), '.pdf'), width=12, height=6)
  errs_loglik_yes_sd <- apply(errs_loglik_yes, 1, sd)
  errs_loglik_yes <- rowMeans(errs_loglik_yes)
  errs_loglik_no_sd <- apply(errs_loglik_no, 1, sd)
  errs_loglik_no <- rowMeans(errs_loglik_no)
  errs_dev_yes_sd <- apply(errs_dev_yes, 1, sd)
  errs_dev_yes <- rowMeans(errs_dev_yes)
  errs_dev_no_sd <- apply(errs_dev_no, 1, sd)
  errs_dev_no <- rowMeans(errs_dev_no)
  ylim_loglik <- c(min(min(errs_loglik_yes), min(errs_loglik_no)), max(max(errs_loglik_yes), max(errs_loglik_no)))
  ylim_loglik <- c(1100, 1300)
  plot(errs_loglik_yes, typ='l', ylab="Log-likelihood", xlab="Iteration", ylim=ylim_loglik, xlim=c(2, 120))
  lines(errs_loglik_no, typ='l', col='red')
  #dev.off()
  #pdf(paste0(paste(scenario, "cv", "deviance", sep='_'), '.pdf'), width=12, height=6)
  ylim_dev <- c(min(min(errs_dev_yes), min(errs_dev_no)), max(max(errs_dev_yes), max(errs_dev_no)))
  ylim_dev <- c(-200, 200)
  plot(errs_dev_yes, typ='l', ylab="Deviance", xlab="Iteration", ylim=ylim_dev, xlim=c(2, 120))
  lines(errs_dev_yes - 1.96*errs_dev_yes_sd)
  lines(errs_dev_yes + 1.96*errs_dev_yes_sd)
  lines(errs_dev_no, typ='l', col='red')
  lines(errs_dev_no - 1.96*errs_dev_no_sd, col='red')
  lines(errs_dev_no + 1.96*errs_dev_no_sd, col='red')
  #dev.off()
}
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

total_summary_no <- read.csv('../dataset/scenario-correlated/non-cyclic-no-intercept-results/non-cyclic-no-intercept_correlated_total_summary.csv')
total_summary_yes <- read.csv('../dataset/scenario-correlated/non-cyclic-intercept-results/non-cyclic-intercept_correlated_total_summary.csv')

# Plot deviance
#ylim_deviance <- c(min(total_summary$deviance), max(total_summary$deviance))
filename <- paste0(paste(scenario, algorithm, "cv", "deviance", sep='_'), '.pdf')
filename <- paste0('../../text/figures/', filename)
#pdf(filename, width=12, height=6)
#plot(total_summary$deviance, ylim=ylim_deviance, xlab="Seed", ylab="Deviance")
#hist(total_summary_yes$deviance, ylim=c(0, 40), col=rgb(red=1,green=0,blue=0, alpha=0.2), breaks=seq(-100, 250, 10), xlab="Deviance", main="Deviances on test set")
#hist(total_summary_no$deviance, ylim=c(0, 40), col=rgb(red=0,green=0,blue=1, alpha=0.2), breaks=seq(-100, 250, 10), xlab="Deviance", main="Deviances on test set", add=T)
abline(v=0, col='red', lwd=3)
#dev.off()
