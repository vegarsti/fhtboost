rm(list=ls())

library(devtools)
library(profvis) # pause
library(foreach)
#library(doParallel)
load_all()

# Get simulated data
# setup_type is one of 'small_dense', 'small_sparse', 'huge', 'huge_clinical', 'correlated'
simulated_data <- simulate_FHT_data(N=500, setup_type='correlated', add_noise=FALSE)
times <- simulated_data$observations$survival_times
delta <- simulated_data$observations$delta
X <- simulated_data$design_matrices$X
Z <- simulated_data$design_matrices$Z
beta_true <- simulated_data$true_parameters$beta
gamma_true <- simulated_data$true_parameters$gamma

find_joint_maximum <- FALSE

# Kaplan-Meier plot
non_para <- non_parametric_estimates(times, delta, continuous = TRUE)
plot(non_para$times_sequence, non_para$kaplan_meiers, typ='s', xlab="Time", ylab="Kaplan-Meier estimated survival probability")

# Dimensions
d <- dim(X)[2]
p <- dim(Z)[2]
N <- dim(X)[1]

# Cross validation
do_CV <- FALSE

if (do_CV) {
  # DIVIDE INTO K FOLDS
  K <- 5
  K_fold_repetitions <- 1 # or 10
  M <- 100 # should be guaranteed in over fitting space
  CV_result <- run_CV(M, K_fold_repetitions, K, X, Z, times, delta)
  CV_errors <- CV_result$CV_errors
  CV_errors_k <- CV_result$CV_errors_k
  plot(CV_errors, typ='l')
  m_stop <- which.min(CV_errors)
  # plot CV results
  plot(CV_errors, typ='l', lty=3)
  for (k in 1:K_fold_repetitions) {
    lines(CV_errors_k[, k], lty=3)
  }
} else {
  m_stop <- 60 # or some other value; needs to be the minimizer (somewhat)
}

# Loglikelihood of the null model:
best_intercepts <- maximum_likelihood_intercepts(times, delta)
y0 <- best_intercepts[1]
mu <- best_intercepts[2]
null_model_loglikelihood <- - sum(FHT_loglikelihood_with_y0_mu(y0, mu, times, delta))

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

tt <- Sys.time()
result <- boosting_run(times, delta, X, Z, m_stop, boost_intercepts_continually=TRUE, should_print=FALSE, run_in_parallel=FALSE)
cat('time: ', Sys.time() - tt)

#result_no_intercept_boosting <- boosting_run(times, delta, X, Z, m_stop, boost_intercepts_continually=FALSE, should_print=FALSE)
#result_more_steps <- boosting_run(times, delta, X, Z, m_stop+100, boost_intercepts_continually=TRUE)
#result_no_intercept_boosting <- boosting_run(times, delta, X, Z, m_stop+50, boost_intercepts_continually=FALSE)

### PLOTTING ###

# settings / meta data
ylim_vector <- c(min(result$loss, maximum_likelihood - 100), result$loss[1] + 100)
plot_title <- ''
xlabel <- 'Iteration'
ylabel <- 'Negative log likelihood'

plot(result$loss, typ='l', lty=1, main=plot_title, xlab=xlabel, ylab=ylabel)#, ylim=ylim_vector)
#lines(result_no_intercept_boosting$loss, lty=3)
#abline(h=maximum_likelihood, col='red')

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
