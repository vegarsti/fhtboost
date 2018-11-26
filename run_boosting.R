# Get simulated data
simulated_data <- simulate_FHT_data()
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
do_CV <- FALSE

if (do_CV) {
  # DIVIDE INTO K FOLDS
  K <- 10
  K_fold_repetitions <- 10
  M <- 50 ### M STOP
  CV_errors <- run_CV(M, K_fold_repetitions, K, X, Z, times, delta)
  loss <- rowSums(loss_K)
  plot(CV_errors, typ='l')
  m_stop <- which.min(CV_errors)
} else {
  m_stop <- 15
}

# Loglikelihood of the null model:
null_model_loglikelihood <- FHT_minus_loglikelihood_with_all_parameters(beta_=rep(0, d), gamma_=rep(0, p), X, Z, times, delta)

result <- boosting_run(times, delta, X, Z, m_stop)
result_more_steps <- boosting_run(times, delta, X, Z, m_stop+50)

plot(result$loss, typ='l', lty=2, main='Different loss function results', ylim=c(1200, 1500), xlab='Iteration', ylab='Loss function')
plot(result_more_steps$loss, typ='l', lty=2, main='Different loss function results', ylim=c(1200, 1500), xlab='Iteration', ylab='Loss function')
abline(v=m_stop, col='red', lwd=2)
