rm(list=ls())

library(devtools)
load_all()

# Get simulated data
simulated_data <- simulate_FHT_data(dense=TRUE) # kind = 'normal' or 'uniform'
times <- simulated_data$observations$survival_times
delta <- simulated_data$observations$delta
X <- simulated_data$design_matrices$X
Z <- simulated_data$design_matrices$Z
beta_true <- simulated_data$true_parameters$beta
gamma_true <- simulated_data$true_parameters$gamma

# Kaplan-Meier plot
non_para <- non_parametric_estimates(times, delta, continuous = TRUE)
plot(non_para$times_sequence, non_para$kaplan_meiers, typ='s')

# Dimensions
d <- dim(X)[2]
p <- dim(Z)[2]
N <- dim(X)[1]

# Cross validation
do_CV <- TRUE

if (do_CV) {
  # DIVIDE INTO K FOLDS
  K <- 10
  K_fold_repetitions <- 5 # or 10
  M <- 100 # should be guaranteed in over fitting space
  CV_result <- run_CV(M, K_fold_repetitions, K, X, Z, times, delta)
  CV_errors <- CV_result$CV_errors
  CV_errors_k <- CV_result$CV_errors_k
  plot(CV_errors, typ='l')
  m_stop <- which.min(CV_errors)
} else {
  m_stop <- 21 # or some other value; needs to be the minimizer (somewhat)
}

# plot CV results
plot(CV_errors, typ='l', lty=3)
for (k in 1:K_fold_repetitions) {
  lines(CV_errors_k[, k], lty=3)
}

# Loglikelihood of the null model:
best_intercepts <- maximum_likelihood_intercepts(times, delta)
y0 <- best_intercepts[1]
mu <- best_intercepts[2]
null_model_loglikelihood <- - sum(FHT_loglikelihood_with_y0_mu(y0, mu, times, delta))

### FIND MAX ###
minus_FHT_loglikelihood <- data_to_optimizable_function(X, Z, times, delta)

p <- dim(X)[2]
d <- dim(Z)[2]

# Run optimization
#initial_parameters <- runif(p+d, min=-0.5, max=0.5) # may need to adjust to get a viable initial value
initial_parameters <- rep(0.1, p+d)
nlm_result <- nlm(minus_FHT_loglikelihood, initial_parameters)
optimized_parameters <- nlm_result$estimate
maximum_likelihood <- nlm_result$minimum

# DO BOOSTING
result <- boosting_run(times, delta, X, Z, m_stop, boost_intercepts_continually=TRUE)
result_more_steps <- boosting_run(times, delta, X, Z, m_stop+50, boost_intercepts_continually=TRUE)
result_no_intercept_boosting <- boosting_run(times, delta, X, Z, m_stop+50, boost_intercepts_continually=FALSE)

### PLOTTING ###

# settings / meta data
ylim_vector <- c(1400, 2000)
plot_title <- 'Different loss function results'
xlabel <- 'Iteration'
ylabel <- 'Loss function'

plot(result$loss, typ='l', lty=2, main=plot_title, xlab=xlabel, ylab=ylabel, ylim=ylim_vector)
abline(h=maximum_likelihood, col='red')
plot(result_more_steps$loss, typ='l', lty=2, main=plot_title, xlab=xlabel, ylab=ylabel, ylim=ylim_vector)
lines(result_no_intercept_boosting$loss, lty=3)
abline(v=m_stop, col='red', lwd=2)
abline(h=maximum_likelihood, col='red')
