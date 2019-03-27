library(devtools)
# install_github("vegarsti/fhtboost")
library(foreach)
library(doParallel)
library(fhtboost)
# load_all()
oberthur_filename <- 'preproc_Oberthur_data.Rdata'
load(oberthur_filename)
has_age_observations <- which(!is.na(clinicalData[, 4]))

directory <- "oberthur_clinical/"

X <- as.matrix(scale(molecularData[has_age_observations, ]))
Z <- as.matrix(scale(clinicalData[has_age_observations, c(3, 4)]))
times <- clinicalData$time[has_age_observations]
delta <- clinicalData$status[has_age_observations]

no_cores <- detectCores() - 1
registerDoParallel(cores=no_cores)
cl <- makeCluster(no_cores)

# OPTIONS
K_fold_repetitions <- 10
K <- 5
boost_intercepts_continually <- FALSE
run_clinical <- FALSE
run_both <- FALSE
run_genetic <- TRUE

seeds <- 1:100

foreach(seed=seeds) %dopar% {
  set.seed(seed)
  seed_string <- formatC(seed, width=2, flag="0")

  # Divide into test and train. test approx 1/3
  num_folds <- 3
  folds <- create_folds_stratified(delta, num_folds)
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

  ### INITIAL ANALYSIS
  non_para <- non_parametric_estimates(times_train, delta_train, continuous = TRUE)
  full_filename <- paste0(directory, seed_string, "_", "kaplan_meier.pdf")
  pdf(full_filename, width=12, height=6)
  plot(
    non_para$times_sequence, non_para$kaplan_meiers, typ='s', xlab="Time",
    ylab="Kaplan-Meier estimated survival probability", ylim=c(0, 1)
  )
  dev.off()


  # Both
  if (run_both) {
    M <- m_stop <- 100 # ??
    CV_result <- run_CV(
      M, K_fold_repetitions, K, X_train, Z_train, times_train, delta_train,
      boost_intercepts_continually=boost_intercepts_continually
    )

    ### WRITE CV RESULT TO FILE
    full_filename <- paste0(directory, seed_string, "_", "loglik.csv")
    write.csv(CV_result$CV_errors_K_loglik, file=full_filename, row.names=FALSE)
    full_filename <- paste0(directory, seed_string, "_", "deviance.csv")
    write.csv(CV_result$CV_errors_K_deviance, file=full_filename, row.names=FALSE)


    ## READ CV RESULT FROM FILE
    seed_string <- formatC(seed, width=2, flag="0")
    full_filename <- paste0(directory, seed_string, "_", "loglik.csv")
    logliks <- read.csv(full_filename)
    full_filename <- paste0(directory, seed_string, "_", "deviance.csv")
    deviances <- read.csv(full_filename)

    ### POST PROCESSING AND PLOTTING
    logliks <- CV_result$CV_errors_K_loglik
    ylims <- c(min(apply(logliks, 2, min)), max(apply(logliks, 2, max)))
    m_stop_from_CV <- which.min(rowMeans(logliks))
    #m_stop_from_CV <- which.max(rowMeans(deviances))
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
  }









  if (run_clinical) {
    # Clinical
    M_clinical <- 10
    CV_result_clinical <- run_CV_clinical(
      M_clinical, K_fold_repetitions, K, X_train, Z_train, times_train, delta_train,
      boost_intercepts_continually=boost_intercepts_continually
    )

    boosting_type <- "clinical"

    ### WRITE CV RESULT TO FILE
    full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "loglik.csv")
    write.csv(CV_result_clinical$CV_errors_K_loglik, file=full_filename, row.names=FALSE)
    full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "deviance.csv")
    write.csv(CV_result_clinical$CV_errors_K_deviance, file=full_filename, row.names=FALSE)


    ## READ CV RESULT FROM FILE
    seed_string <- formatC(seed, width=2, flag="0")
    full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "loglik.csv")
    logliks <- read.csv(full_filename)
    full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "deviance.csv")
    deviances <- read.csv(full_filename)

    ### POST PROCESSING AND PLOTTING
    logliks <- CV_result_clinical$CV_errors_K_loglik
    ylims <- c(min(apply(logliks, 2, min)), max(apply(logliks, 2, max)))
    m_stop_from_CV <- which.min(rowMeans(logliks))
    #m_stop_from_CV <- which.max(rowMeans(deviances))
    full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "loglik.pdf")
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

    result <- cyclic_boosting_run(
      times=times_train,
      delta=delta_train,
      X=X_train,
      Z=Z_train,
      m_stop_y0=1,
      m_stop_mu=m_stop_from_CV,
      boost_intercepts_continually=boost_intercepts_continually,
      should_print=FALSE
    )

    beta_hat <- result$final_parameters$beta_hat_final
    gamma_hat <- result$final_parameters$gamma_hat_final
    y0_hat <- exp(X_train %*% beta_hat)
    mu_hat <- Z_train %*% gamma_hat

    betas <- data.frame(cbind(non_null_parameters(beta_hat) - 1, beta_hat[non_null_parameters(beta_hat)]))
    names(betas) <- c("j", "beta_j")
    full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "beta.csv")
    write.csv(betas, file=full_filename, row.names=FALSE)
    gammas <- data.frame(cbind((1:length(gamma_hat)) - 1, gamma_hat))
    names(gammas) <- c("j", "gamma_j")
    full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "gamma.csv")
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
    full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "test_result.csv")
    write.csv(loglikelihood_df, file=full_filename, row.names=FALSE)
  }


  if (run_genetic) {
    # Genetic
    M_genetic <- 50
    CV_result_genetic <- run_CV_genetic(
      M_genetic, K_fold_repetitions, K, X_train, Z_train, times_train, delta_train,
      boost_intercepts_continually=boost_intercepts_continually
    )

    boosting_type <- "genetic"


    ### WRITE CV RESULT TO FILE
    full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "loglik.csv")
    write.csv(CV_result_genetic$CV_errors_K_loglik, file=full_filename, row.names=FALSE)
    full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "deviance.csv")
    write.csv(CV_result_genetic$CV_errors_K_deviance, file=full_filename, row.names=FALSE)


    ## READ CV RESULT FROM FILE
    seed_string <- formatC(seed, width=2, flag="0")
    full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "loglik.csv")
    logliks <- read.csv(full_filename)
    full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "deviance.csv")
    deviances <- read.csv(full_filename)

    ### POST PROCESSING AND PLOTTING
    logliks <- CV_result_genetic$CV_errors_K_loglik
    ylims <- c(min(apply(logliks, 2, min)), max(apply(logliks, 2, max)))
    m_stop_from_CV <- which.min(rowMeans(logliks))
    #m_stop_from_CV <- which.max(rowMeans(deviances))
    full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "loglik.pdf")
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

    result <- cyclic_boosting_run(
      times=times_train,
      delta=delta_train,
      X=X_train,
      Z=Z_train,
      m_stop_y0=m_stop_from_CV,
      m_stop_mu=1,
      boost_intercepts_continually=boost_intercepts_continually,
      should_print=FALSE
    )

    beta_hat <- result$final_parameters$beta_hat_final
    gamma_hat <- result$final_parameters$gamma_hat_final
    y0_hat <- exp(X_train %*% beta_hat)
    mu_hat <- Z_train %*% gamma_hat

    betas <- data.frame(cbind(non_null_parameters(beta_hat) - 1, beta_hat[non_null_parameters(beta_hat)]))
    names(betas) <- c("j", "beta_j")
    full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "beta.csv")
    write.csv(betas, file=full_filename, row.names=FALSE)
    gammas <- data.frame(cbind((1:length(gamma_hat)) - 1, gamma_hat))
    names(gammas) <- c("j", "gamma_j")
    full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "gamma.csv")
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
    full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "test_result.csv")
    write.csv(loglikelihood_df, file=full_filename, row.names=FALSE)
  }


}
stopCluster(cl)
