rm(list=ls())
library(devtools)
load_all()

# num_seeds <- 100
# seeds <- 1:num_seeds
#
# maximum_likelihoods <- rep(0, num_seeds)
# loss_intercept <- c()
# loss_no_intercept <- c()

#for (seed in seeds) {
seed <- 10
simulated_data <- simulate_FHT_data(setup_type='small_dense', add_noise=FALSE, seed=seed)
times <- simulated_data$observations$survival_times
delta <- simulated_data$observations$delta
X <- simulated_data$design_matrices$X
Z <- simulated_data$design_matrices$Z
beta_true <- simulated_data$true_parameters$beta
gamma_true <- simulated_data$true_parameters$gamma

# Kaplan-Meier plot
non_para <- non_parametric_estimates(times, delta, continuous = TRUE)
plot(non_para$times_sequence, non_para$kaplan_meiers, typ='s', xlab="Time", ylab="Kaplan-Meier estimated survival probability")

m_stop <- 35

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
result_continually <- boosting_run(times, delta, X, Z, m_stop, boost_intercepts_continually=TRUE, should_print=FALSE)

loss_intercept <- result_continually$loss
loss_no_intercept <- result$loss

# settings / meta data
ylim_vector <- c(min(result$loss, maximum_likelihood - 100), result$loss[1] + 100)
plot_title <- ''
xlabel <- 'Iteration'
ylabel <- 'Negative log likelihood'
ylim_vector <- c(maximum_likelihood-40, max(result$loss))


filename <- 'case1_fixed_only.pdf'
filename <- paste0('../../text/figures/', filename)
pdf(filename, width=12, height=6)
plot(result$loss, typ='l', lty=1, main=plot_title, xlab=xlabel, ylab=ylabel, ylim=ylim_vector)
abline(h=maximum_likelihood, col='blue', lwd=2)
legend(
  'topright',
  legend = c("Boosting with fixed intercept", "Max log-likelihood"),
  col = c('black', 'blue'),
  lty = c(1, 2, 1),
  lwd = c(1, 1, 2)
)
dev.off()

#}
###
filename <- 'small_example.pdf'
filename <- paste0('../../text/figures/', filename)
#pdf(filename, width=12, height=6)
plot_title <- "Boosting with changing intercept recovers the maximum likelihood estimates"
plot(result$loss, typ='l', lty=1, main=plot_title, xlab=xlabel, ylab=ylabel, ylim=ylim_vector)
lines(result_continually$loss, lty=2, col='red')
abline(h=maximum_likelihood, col='blue', lwd=2)
legend(
  'topright',
  legend = c("Boosting with fixed intercept", "Boosting with changing of intercept", "Max log-likelihood"),
  col = c('black', 'red', 'blue'),
  lty = c(1, 2, 1),
  lwd = c(1, 1, 2)
)
#dev.off()
###









### BRIER
beta_hat <- result_continually$final_parameters$beta_hat_final
gamma_hat <- result_continually$final_parameters$gamma_hat_final
y0_hat <- exp(X %*% beta_hat)
mu_hat <- Z %*% gamma_hat
briers <- sapply(times, function(time) { brier_score_with_censoring(time, times, delta, y0s=y0_hat, mus=mu_hat) })
brier_score_with_censoring_on_times
