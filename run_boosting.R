rm(list=ls())

library(devtools)
library(profvis) # pause
library(foreach)
library(readr)
load_all()

seed <- 10
#for (seed in 3:5) {
simulated_data <- simulate_FHT_data(N=500, setup_type='correlated', add_noise=FALSE, seed=seed)
times <- simulated_data$observations$survival_times
delta <- simulated_data$observations$delta
X <- simulated_data$design_matrices$X
Z <- simulated_data$design_matrices$Z
beta_true <- simulated_data$true_parameters$beta
gamma_true <- simulated_data$true_parameters$gamma

# Kaplan-Meier plot
non_para <- non_parametric_estimates(times, delta, continuous = TRUE)
plot(non_para$times_sequence, non_para$kaplan_meiers, typ='s', xlab="Time", ylab="Kaplan-Meier estimated survival probability")

# Cross validation
do_CV <- TRUE
boost_intercepts_continually <- FALSE
if (do_CV) {
  # DIVIDE INTO K FOLDS
  K <- 5
  K_fold_repetitions <- 2 # or 10
  M <- 100 # should be guaranteed in over fitting space
  CV_result <- run_CV(M, K_fold_repetitions, K, X, Z, times, delta, boost_intercepts_continually)
  CV_errors_K <- CV_result$CV_errors_K_deviance
  CV_errors <- rowMeans(CV_errors_K)
  y_max <- max(apply(CV_errors_K, 1, max))
  plot(CV_errors, typ='l', ylim=c(0, y_max))
  # plot CV results
  Ks <- dim(CV_errors_K)[2]
  for (k in 1:Ks) {
    lines(CV_errors_K[, k], lty=3)
  }
  m_stop <- which.max(CV_errors)
} else {
  m_stop <- 60 # or some other value; needs to be the minimizer (somewhat)
}
print(m_stop)
#}

# Dimensions
d <- dim(X)[2]
p <- dim(Z)[2]
N <- dim(X)[1]

# Loglikelihood of the null model:
best_intercepts <- maximum_likelihood_intercepts(times, delta)
y0 <- best_intercepts[1]
mu <- best_intercepts[2]
null_y0 <- y0
null_mu <- mu
null_model_loglikelihood <- - sum(FHT_loglikelihood_with_y0_mu(y0, mu, times, delta))

find_joint_maximum <- TRUE
if (find_joint_maximum) {
  ### FIND MAX ###
  minus_FHT_loglikelihood <- data_to_optimizable_function(X, Z, times, delta)

  p <- dim(X)[2]
  d <- dim(Z)[2]

  # Run optimization
  initial_parameters <- runif(p+d, min=-0.5, max=0.5) # may need to adjust to get a viable initial value
  initial_parameters <- rep(0, p+d)
  nlm_result <- nlm(minus_FHT_loglikelihood, initial_parameters)
  optimized_parameters <- nlm_result$estimate
  maximum_likelihood <- nlm_result$minimum
} else {
  maximum_likelihood <- 2000
}
# DO BOOSTING

result <- boosting_run(times, delta, X, Z, m_stop, boost_intercepts_continually=FALSE, should_print=FALSE)

estimated_y0s <- exp(X %*% result$final_parameters$beta_hat_final)
estimated_mus <- Z %*% result$final_parameters$gamma_hat_final
#estimated_survival <- FHT_parametric_survival(times, mu, y0)

# Plot Brier scores
# tt <- Sys.time()
# brier_result <- brier_score_on_censored_data(times, delta, estimated_y0s, estimated_mus)
# cat('time: ', Sys.time() - tt)
# plot(brier_result$brier_times, brier_result$brier_scores, typ='l')
#
# # Plot Brier R2
# tt <- Sys.time()
# r2_result <- brier_r2(times, delta, estimated_y0s, estimated_mus, null_y0 = y0, null_mu = mu)
# cat('time: ', Sys.time() - tt)
# plot(r2_result$r2_times, r2_result$r2, typ='l', ylim=c(0, 1))


### PLOTTING ###

# settings / meta data
ylim_vector <- c(min(result$loss, maximum_likelihood - 100), result$loss[1] + 100)
plot_title <- ''
xlabel <- 'Iteration'
ylabel <- 'Negative log likelihood'

plot(result$loss, typ='l', lty=1, main=plot_title, xlab=xlabel, ylab=ylabel)#, ylim=ylim_vector)
#lines(result_no_intercept_boosting$loss, lty=3)
abline(h=maximum_likelihood, col='red')

y0_post <- exp(result$final_parameters$gamma_hat_final[1])
mu_post <- result$final_parameters$beta_hat_final[1]
#null_model_loglikelihood_post <- - sum(FHT_loglikelihood_with_y0_mu(y0_post, mu, times, delta))
# doesn't work because NaNs produced?

colors <- c('black', 'black', 'red', 'red')
ltypes <- c(2, 3, 1, 1)
lwd <- c(1, 1, 1, 3)
ylim_vector <- c(min(result_more_steps$loss, maximum_likelihood - 100), result$loss[1] + 100)
plot(result_more_steps$loss, typ='l', lty=2, main=plot_title, xlab=xlabel, ylab=ylabel, ylim=ylim_vector)
lines(result_no_intercept_boosting$loss, lty=3)
abline(h=maximum_likelihood, col='red')
abline(v=m_stop, col='red', lwd=3)
legend(
  x='topright',
  legend=c('Boosting with changing the intercept', 'Boosting without changing the intercept', 'Maximum likelihood value', 'm_stop'),
  col=colors,
  lty = ltypes,
  lwd = lwd
)


# Plot parameters
plot(abs(result$final_parameters$gamma_hat_final), ylim=c(0, 0.2))
abline(h=0.1)
abline(v=11)

plot(abs(result$final_parameters$beta_hat_final), ylim=c(0, 0.2))
abline(h=0.1)
abline(v=11)
