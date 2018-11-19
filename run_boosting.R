simulated_data <- simulate_FHT_data()
times <- simulated_data$observations$survival_times
delta <- simulated_data$observations$delta
X <- simulated_data$design_matrices$X
Z <- simulated_data$design_matrices$Z
beta_true <- simulated_data$true_parameters$beta
gamma_true <- simulated_data$true_parameters$gamma

# non_para <- non_parametric_estimates(times, delta, continuous = TRUE)
# plot(non_para$times_sequence, non_para$kaplan_meiers, typ='s')

# dimensions
d <- dim(X)[2]
p <- dim(Z)[2]
N <- dim(X)[1]

### SANITY CHECK -> does nlm recover the parameters? ###
minus_FHT_loglikelihood_nlm <- data_to_optimizable_function(X, Z, times, delta)
# # Run optimization
initial_parameters <- runif(p+d, min=0.1, max=0.5)
nlm_result <- nlm(minus_FHT_loglikelihood_nlm, initial_parameters)
nlm_parameter_list <- parameter_vector_to_list(nlm_result$estimate, d, p)
beta_from_nlm <- nlm_parameter_list$beta
beta_0_from_nlm <- beta_from_nlm[1]
gamma_from_nlm <- nlm_parameter_list$gamma
gamma_0_from_nlm <- gamma_from_nlm[1]

# DIVIDE INTO K FOLDS
K <- 10
K_fold_repetitions <- 10
M <- 50 ### M STOP
CV_errors <- run_CV(M, K_fold_repetitions, K, X, Z, times, delta)
loss <- rowSums(loss_K)
plot(CV_errors, typ='l')
m_stop <- which.min(CV_errors)

result <- boosting_run(times, delta, X, Z, m_stop)
result_more_steps <- boosting_run(times, delta, X, Z, m_stop+30)

# Plot loss functions
plot(result$loss, typ='l', lty=2, main='Different loss function results', ylim=c(1200, 1500), xlab='Iteration', ylab='Loss function')
plot(result_more_steps$loss, typ='l', lty=2, main='Different loss function results', ylim=c(1200, 1500), xlab='Iteration', ylab='Loss function')
abline(h=nlm_loss, col='blue')
abline(v=m_stop, col='blue', lty=3)
legend("topright", legend=c("With", "Without", "ML intercept, no changing", "ML (nlm)", "m_stop"), col=c("black", "red", "black", "blue", "blue"), lty=c(2, 1, 3, 1, 3), cex=1.5)

# Plot beta0
joint_optimized <- nlm_result$estimate[1]
beta0s <- result_wo_nlm$parameters$beta_hats[, 1]
plot(beta0s, typ='l', ylim=c(0.5, 2), lty=2, main='Beta intercept (in y0)')
abline(h=joint_optimized, col='red')
abline(h=beta0s[1], col='blue')
abline(v=m_stop_mu, col='blue', lty=3)
legend("topleft", legend=c("Path", "Initial", "ML", "m_stop"), col=c("black", "red", "blue", "blue"), lty=c(2, 1, 1, 3), cex=1.5)

# Plot gamma0
joint_optimized <- nlm_result$estimate[d+1]
gamma0s <- result_wo_nlm$parameters$gamma_hats[, 1]
plot(gamma0s, typ='l', ylim=c(-1, 0), lty=2, main='Gamma intercept (in mu)')
abline(h=joint_optimized, col='red')
abline(h=gamma0s[1], col='blue')
abline(v=m_stop_mu, col='blue', lty=3)
legend("topleft", legend=c("Path", "Initial", "ML", "m_stop"), col=c("black", "red", "blue", "blue"), lty=c(2, 1, 1, 3), cex=1.5)
