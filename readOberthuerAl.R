#################
# read the data #
#################
directory <- 'neuroblastoma/'
filename <- paste0(directory, 'webclinical.txt')
tmp <- read.table(filename, header = TRUE, sep = '\t') # data collected from the web
tmp <- tmp[seq(1, 751, by = 2), ]
tmp <- tmp[ , -c(2, 3, 4, 6, 7)]
tmp[ , 1] <- substr(tmp[ , 1], 1, 4)
tmp[ , 3] <- substr(tmp[ , 3], 1, 4)
tmp[ , 4] <- substr(tmp[ , 4], 1, 4)
tmp[ , 5] <- substr(tmp[ , 5], 1, 5)
for (i in 1:376)
{
  if (is.na(suppressWarnings(as.numeric(tmp[i, 1])))) tmp[i, 1] <- substr(tmp[i, 1], 1, 3)
  if (is.na(suppressWarnings(as.numeric(tmp[i, 1])))) tmp[i, 1] <- substr(tmp[i, 1], 1, 2)
  if (is.na(suppressWarnings(as.numeric(tmp[i, 1])))) tmp[i, 1] <- substr(tmp[i, 1], 1, 1)
  if (is.na(suppressWarnings(as.numeric(tmp[i, 3])))) tmp[i, 3] <- substr(tmp[i, 3], 1, 3)
  if (is.na(suppressWarnings(as.numeric(tmp[i, 3])))) tmp[i, 3] <- substr(tmp[i, 3], 1, 2)
  if (is.na(suppressWarnings(as.numeric(tmp[i, 3])))) tmp[i, 3] <- substr(tmp[i, 3], 1, 1)
  if (is.na(suppressWarnings(as.numeric(tmp[i, 4])))) tmp[i, 4] <- substr(tmp[i, 4], 1, 3)
  if (is.na(suppressWarnings(as.numeric(tmp[i, 4])))) tmp[i, 4] <- substr(tmp[i, 4], 1, 2)
  if (is.na(suppressWarnings(as.numeric(tmp[i, 4])))) tmp[i, 4] <- substr(tmp[i, 4], 1, 1)
  ifelse(tmp[i, 5] == 'alive', tmp[i, 5] <- 0, tmp[i, 5] <- 1)
}
tmp[ , 1] <- as.numeric(tmp[ , 1])
tmp[ , 3] <- as.numeric(tmp[ , 3])
tmp[ , 4] <- as.numeric(tmp[ , 4])
medianAge <- median(tmp$Age)

filename <- paste0(directory, 'NB2004_clin_and_genes_neuroblastoma_noTies.txt')
boev <- read.table(filename, header = TRUE, sep = '\t')
# data provided by H. M. Boevelstad
clinicalData <- cbind(boev[ , 1:3], rep(NA, dim(boev)[1]))
supp <- tmp[c(288:329, 331, 332, 335, 338, 339, 341, 343:353, 355:369, 371:376), 4]
clinicalData[283:362, 4] <- as.numeric(supp > medianAge) # the last observations are in order
for (i in 1:282)
  for (j in 1:287)
    if ((trunc(boev[i, 1] * 365) == tmp[j, 4]) & boev[i, 2] == tmp[j, 5])
      clinicalData[i, 4] <- as.numeric(tmp[j, 1] > medianAge)
# survival time and censoring status are unique key. Only possible issue: observation 247.
colnames(clinicalData) <- c('time', 'status', 'risk', 'age')

molecularData <- boev[, -c(1:3)]
colnames(molecularData) <- paste('X', 1:9978, sep = '')

rm(tmp, boev, i, j, medianAge, supp)

oberthur_filename <- 'preproc_Oberthur_data.Rdata'
save(molecularData, clinicalData, file=oberthur_filename)

### HERE MY WORK STARTS (data prep)
library(devtools)
install_github("vegarsti/fhtboost")
library(fhtboost) # load_all()
oberthur_filename <- 'preproc_Oberthur_data.Rdata'
load(oberthur_filename)
has_age_observations <- which(!is.na(clinicalData[, 4]))

X <- as.matrix(scale(molecularData[has_age_observations, ]))
Z <- as.matrix(scale(clinicalData[has_age_observations, c(3, 4)]))
times <- clinicalData$time[has_age_observations]
delta <- clinicalData$status[has_age_observations]

# Divide into train and test
seed <- 2
set.seed(seed)
# Divide into test and train. test approx 1/3
K <- 3
folds <- create_folds_stratified(delta, K)
test_indices <- sort(folds[[1]])
train_indices <- sort(c(folds[[2]], folds[[3]]))

## TRAIN
ones_train <- rep(1, length(train_indices))
times_train <- times[train_indices]
delta_train <- delta[train_indices]
X_train_rest <- X[train_indices, ]
X_train <- as.matrix(cbind(ones_train, X_train_rest))
Z_train_rest <- Z[train_indices, ]
Z_train <- as.matrix(cbind(ones_train, Z_train_rest))

## TEST
ones_test <- rep(1, length(test_indices))
times_test <- times[test_indices]
delta_test <- delta[test_indices]
X_test_rest <- X[test_indices, ]
X_test <- as.matrix(cbind(ones_test, X_test_rest))
Z_test_rest <- Z[test_indices, ]
Z_test <- as.matrix(cbind(ones_test, Z_test_rest))


# remove "NA" age observations
# na_age_observations <- which(is.na(clinicalData[, 4]))
# indices_to_remove_train <- c()
# indices_to_remove_test <- c()
# for (i in 1:length(train_indices)) {
#   if (any(train_indices[i] == na_age_observations)) {
#     indices_to_remove_train <- c(indices_to_remove_train, i)
#   }
# }
# for (i in 1:length(test_indices)) {
#   if (any(test_indices[i] == na_age_observations)) {
#     indices_to_remove_test <- c(indices_to_remove_test, i)
#   }
# }
# train_indices <- train_indices[-indices_to_remove_train]
# test_indices <- test_indices[-indices_to_remove_test]

# ones <- rep(1, length(train_indices))
#
# times_train <- clinicalData$time[train_indices]
# delta_train <- clinicalData$status[train_indices]
# X_train_rest <- scale(molecularData[train_indices, ])
# X_train <- cbind(ones, X_train_rest)
# Z_train_rest <- as.matrix(scale(clinicalData[, c(3, 4)])[train_indices, ])
# Z_train <- as.matrix(cbind(ones, Z_train_rest))

### INITIAL ANALYSIS

non_para <- non_parametric_estimates(times_train, delta_train, continuous = TRUE)
plot(non_para$times_sequence, non_para$kaplan_meiers, typ='s', xlab="Time", ylab="Kaplan-Meier estimated survival probability", ylim=c(0, 1))

# find intercepts to start with

M <- m_stop <- 100 # ??
K_fold_repetitions <- 10
K <- 5
boost_intercepts_continually <- FALSE

## RUN CV
CV_result <- run_CV(
  M, K_fold_repetitions, K, X_train, Z_train, times_train, delta_train,
  boost_intercepts_continually=boost_intercepts_continually
)
plot(rowMeans(CV_result$CV_errors_K_loglik), typ='l')


### WRITE CV RESULT TO FILE
seed_string <- formatC(seed, width=2, flag="0")
directory <- "../dataset/oberthuer/"
full_filename <- paste0(directory, seed_string, "_", "loglik.csv")
write.csv(CV_result$CV_errors_K_loglik, file=full_filename, row.names=FALSE)
full_filename <- paste0(directory, seed_string, "_", "deviance.csv")
write.csv(CV_result$CV_errors_K_deviance, file=full_filename, row.names=FALSE)


## READ CV RESULT FROM FILE
seed_string <- formatC(seed, width=2, flag="0")
directory <- "../dataset/oberthuer/"
full_filename <- paste0(directory, seed_string, "_", "loglik.csv")
logliks <- read.csv(full_filename)
full_filename <- paste0(directory, seed_string, "_", "deviance.csv")
deviances <- read.csv(full_filename)

### POST PROCESSING AND PLOTTING
logliks <- CV_result$CV_errors_K_loglik
ylims <- c(min(apply(logliks, 2, min)), max(apply(logliks, 2, max)))
m_stop_from_CV <- which.min(rowMeans(logliks))
#m_stop_from_CV <- which.max(rowMeans(deviances))
directory <- "../dataset/oberthuer/"
full_filename <- paste0(directory, seed_string, "_", "loglik.pdf")
pdf(full_filename, width=12, height=6)
plot(rowMeans(logliks), typ='l', ylim=ylims, ylab="Negative log-likelihood", xlab="Boosting iteration")
for (k in 1:K_fold_repetitions) {
  lines(logliks[, k], lty=3, col=rgb(0, 0, 0, alpha = 0.5))
}
abline(v=m_stop_from_CV, lwd=2, col='red')
legend(
  'topright',
  legend=c("Sum of log-lik. on test set in 5-fold CV", "Mean of log-lik. sums in 5-fold CV", "Iteration number which minimizes mean"),
  col=c(rgb(0, 0, 0, alpha = 0.5), 'black', 'red'),
  lty=c(3, 1, 1),
  lwd=c(1, 1, 2)
)
dev.off()

result <- boosting_run(
  times=times_train,
  delta=delta_train,
  X=X_train,
  Z=Z_train,
  m_stop=m_stop_from_CV,
  boost_intercepts_continually=FALSE,
  should_print=FALSE
)

beta_hat <- result$final_parameters$beta_hat_final
gamma_hat <- result$final_parameters$gamma_hat_final
y0_hat <- exp(X_train %*% beta_hat)
mu_hat <- Z_train %*% gamma_hat

betas <- data.frame(cbind(non_null_parameters(beta_hat) - 1, beta_hat[non_null_parameters(beta_hat)]))
names(betas) <- c("j", "beta_j")
full_filename <- paste0(directory, seed_string, "_", "beta.csv")
write.csv(betas, file=full_filename, row.names=FALSE)
gammas <- data.frame(cbind((1:length(gamma_hat)) - 1, gamma_hat))
names(gammas) <- c("j", "gamma_j")
full_filename <- paste0(directory, seed_string, "_", "gamma.csv")
write.csv(gammas, file=full_filename, row.names=FALSE)



# Run on test set
beta_hat_null <- rep(0, dim(X_test)[2])
beta_hat_null[1] <- beta_hat[1]
gamma_hat_null <- rep(0, dim(Z_test)[2])
gamma_hat_null[1] <- gamma_hat[1]
test_null_loglikelihood <- FHT_minus_loglikelihood_with_all_parameters(
  beta_hat_null, gamma_hat_null, X_test, Z_test, times_test, delta_test
)
test_loglikelihood <- FHT_minus_loglikelihood_with_all_parameters(
  beta_hat, gamma_hat, X_test, Z_test, times_test, delta_test
)
test_difference_of_deviance <- 2*(test_loglikelihood - test_null_loglikelihood)
loglikelihood_df <- data.frame(
  null_loglikelihood=test_null_loglikelihood,
  loglikelihood=test_loglikelihood,
  deviance=test_difference_of_deviance
)
full_filename <- paste0(directory, seed_string, "_", "test_result.csv")
write.csv(loglikelihood_df, file=full_filename, row.names=FALSE)




# i1 <- 211
# y01 <- y0_hat[i1]
# mu1 <- mu_hat[i1]
#
# i2 <- 174
# y02 <- y0_hat[i2]
# mu2 <- mu_hat[i2]
#
#
# times_to_plot <- seq(0.01, 13, by=0.01)
# #survivals <- FHT_parametric_survival(times_to_plot, mu, y0)
# #cum_hazards <- FHT_parametric_cumulative_hazard(times_to_plot, mu, y0)
# hazard1 <- FHT_hazard(times_to_plot, mu1, y01)
# hazard2 <- FHT_hazard(times_to_plot, mu2, y02)
# hazard_base <- FHT_hazard(times_to_plot, mu1, exp(result$final_parameters$beta_hat_final[1]))
# plot(times_to_plot, hazard1, typ='l')
# lines(times_to_plot, hazard2, col='red')
# lines(times_to_plot, hazard_base, col='blue')
#
# # survivals <- rep(0, length(times_to_plot))
# # for (i in 1:length(times_to_plot)) {
# #   survivals[i] <- FHT_parametric_survival(times_to_plot[i], mu, y0)
# # }
# plot(times_to_plot, survivals, typ='l', ylim=c(0, 1))
# plot(times_to_plot, cum_hazards, typ='l', ylim=c(0, 1))
# plot(times_to_plot, hazard, typ='l')
#
#
#
# ### CYCLIC
# cyclic_result_3 <- cyclic_boosting_run(
#   times=times_train,
#   delta=delta_train,
#   X=X_train,
#   Z=Z_train,
#   m_stop_y0=30,
#   m_stop_mu=3,
#   boost_intercepts_continually=FALSE,
#   should_print=FALSE
# )
#
# cyclic_result_2 <- cyclic_boosting_run(
#   times=times_train,
#   delta=delta_train,
#   X=X_train,
#   Z=Z_train,
#   m_stop_y0=30,
#   m_stop_mu=2,
#   boost_intercepts_continually=FALSE,
#   should_print=FALSE
# )
#
# cyclic_result_1 <- cyclic_boosting_run(
#   times=times_train,
#   delta=delta_train,
#   X=X_train,
#   Z=Z_train,
#   m_stop_y0=30,
#   m_stop_mu=1,
#   boost_intercepts_continually=FALSE,
#   should_print=FALSE
# )
#
# cyclic_result_0 <- cyclic_boosting_run(
#   times=times_train,
#   delta=delta_train,
#   X=X_train,
#   Z=Z_train,
#   m_stop_y0=30,
#   m_stop_mu=0,
#   boost_intercepts_continually=FALSE,
#   should_print=FALSE
# )
#
# plot(cyclic_result_3$loss, ylim=c(150, 200), typ='l')
# lines(cyclic_result_2$loss, typ='l', col='red')
# lines(cyclic_result_1$loss, typ='l', col='blue')
#
#
#
# ### COX
# library(mboost)
# library(survival)
# mstop_CV <- 100 # ?? no. of iterations
# cox_cv <- function(seed, boosting_model) {
#   set.seed(seed)
#   return(apply(cvrisk(
#     boosting_model,
#     folds=cv(model.weights(boosting_model),
#              type='kfold',
#              B=10 # B is number of folds
#     ),
#     papply = lapply
#   ),
#   2,
#   mean
#   ))
# }
#
# cox_model <- glmboost(
#   y=Surv(times_train, delta_train),
#   x=as.matrix(cbind(Z_train_rest, X_train_rest)),
#   family=CoxPH()
# )
# par1 <- which.min(apply(sapply(1:mstop_CV, cox_cv, boosting_model=cox_model), 1, mean))
# # "cvrisk" is broken for coxPH...!!!!
#
# cox_df <- data.frame(cbind(X_train_rest, Z_train_rest))
# cox_result <- coxph(Surv(times_train, delta_train) ~ ., data=cox_df, iter.max=0, x=TRUE)
#
# ##### IG FHT on test data
#
# times_test <- clinicalData$time[test_indices]
# delta_test <- clinicalData$status[test_indices]
# X_test_rest <- scale(molecularData[test_indices, ])
# Z_test_rest <- as.matrix(clinicalData[, c(3, 4)][test_indices, ])
# ones_test <- rep(1, length(times_test))
# X_test <- cbind(ones_test, X_test_rest)
# Z_test <- as.matrix(cbind(ones_test, Z_test_rest))
#
# best_intercepts <- maximum_likelihood_intercepts(times_test, delta_test) # ???!!!!
# null_y0 <- best_intercepts[1]
# null_mu <- best_intercepts[2]
#
# beta_null <- c(null_y0, rep(0, (dim(X_test)[2]-1)))
# gamma_null <- c(null_mu, rep(0, (dim(Z_test)[2]-1)))
# null_model_loglikelihood <- - FHT_minus_loglikelihood_with_all_parameters(beta_null, gamma_null, X_test, Z_test, times_test, delta_test)
# full_model_loglikelihood <- - FHT_minus_loglikelihood_with_all_parameters(result$final_parameters$beta_hat_final, result$final_parameters$gamma_hat_final, X_test, Z_test, times_test, delta_test)
# deviance <- 2 * (full_model_loglikelihood - null_model_loglikelihood)
#
#
# #### COX
#
# #install.packages("riskRegression")
# library(pec)
# test_data <- data.frame(cbind(Z_test_rest, X_test_rest))
# prederr <- pec(cox_result, formula=Surv(times_test, delta_test) ~ 1, data=test_data)
# prederr2 <- pec(cox_result, formula=Surv(times_test, delta_test) ~ ., data=test_data)
