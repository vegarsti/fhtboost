library(devtools)
# install_github("vegarsti/fhtboost")
library(foreach)
library(doParallel)
# library(fhtboost)
load_all()
oberthur_filename <- 'preproc_Oberthur_data.Rdata'
load(oberthur_filename)
has_age_observations <- which(!is.na(clinicalData[, 4]))

X <- as.matrix(scale(molecularData[has_age_observations, ]))
Z <- as.matrix(scale(clinicalData[has_age_observations, c(3, 4)]))
times <- clinicalData$time[has_age_observations]
delta <- clinicalData$status[has_age_observations]


###########################


rm(molecularData, clinicalData)

# OPTIONS
K_fold_repetitions <- 10
K <- 5
boost_intercepts_continually <- FALSE

# directory <- "oberthur/"
directory <- "../dataset/oberthuer/oberthur_all/"

# boost_intercepts_continually <- FALSE
# boost_intercepts_continually <- TRUE

# Set up parallel things
no_cores <- detectCores() - 1
registerDoParallel(cores=no_cores)
cl <- makeCluster(no_cores)

seeds <- 1:100
foreach(seed=seeds) %dopar% {
  set.seed(seed)
  seed_string <- formatC(seed, width=3, flag="0")

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

  ## TEST &&& SORT THESE
  times_test <- times[test_indices]
  delta_test <- delta[test_indices]
  X_test_rest <- X[test_indices, ]
  Z_test_rest <- Z[test_indices, ]
  order_times <- order(times_test)
  delta_test <- delta_test[order_times]
  X_test_rest <- X_test_rest[order_times, ]
  Z_test_rest <- Z_test_rest[order_times, ]
  times_test <- sort(times_test)

  ones_test <- rep(1, length(test_indices))
  X_test <- as.matrix(cbind(ones_test, X_test_rest))
  Z_test <- as.matrix(cbind(ones_test, Z_test_rest))

  # Both
  boosting_type <- "both"

    M <- m_stop <- 200
    CV_result <- run_CV(
      M, K_fold_repetitions, K, X_train, Z_train, times_train, delta_train,
      boost_intercepts_continually=boost_intercepts_continually
    )

    full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "loglik.csv")
    # write.csv(CV_result$CV_errors_K_loglik, file=full_filename, row.names=FALSE)
    #full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "deviance.csv")
    # write.csv(CV_result$CV_errors_K_deviance, file=full_filename, row.names=FALSE)
    # logliks <- CV_result$CV_errors_K_loglik

    logliks <- read.csv(full_filename)


    ### POST PROCESSING AND PLOTTING
    ylims <- c(min(apply(logliks, 2, min)), max(apply(logliks, 2, max)))
    m_stop_from_CV <- which.min(rowMeans(logliks))
    #m_stop_from_CV <- which.max(rowMeans(deviances))
    # full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "loglik.pdf")
    # full_filename <- paste0(tex_figures_directory, "example_cv_loglik.pdf")
    # pdf(full_filename, width=12, height=6)
    # plot(rowMeans(logliks), typ='l', ylim=c(80, 100), ylab="Negative log-likelihood", xlab="Boosting iteration", xlim=c(0, 100))
    # for (k in 1:K_fold_repetitions) {
    #   lines(logliks[, k], lty=3, col=rgb(0, 0, 0, alpha = 0.5))
    # }
    # abline(v=m_stop_from_CV, lwd=2, col='red')
    # legend(
    #   'topright',
    #   legend=c("Sum of log-lik. on test set in 5-fold CV", "Mean of log-lik. sums in 5-fold CV", "Iteration number which minimizes mean"),
    #   col=c(rgb(0, 0, 0, alpha = 0.5), 'black', 'red'),
    #   lty=c(3, 1, 1),
    #   lwd=c(1, 1, 2)
    # )
    # dev.off()

    result <- boosting_run(
      times=times_train,
      delta=delta_train,
      X=X_train,
      Z=Z_train,
      m_stop=m_stop_from_CV,
      boost_intercepts_continually=boost_intercepts_continually,
      should_print=FALSE
    )

    beta_hat <- result$final_parameters$beta_hat_final
    gamma_hat <- result$final_parameters$gamma_hat_final
    y0_hat <- exp(X_train %*% beta_hat)
    mu_hat <- Z_train %*% gamma_hat

    # betas <- data.frame(cbind(non_null_parameters(beta_hat) - 1, beta_hat[non_null_parameters(beta_hat)]))
    # names(betas) <- c("j", "beta_j")
    # full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "beta.csv")
    # write.csv(betas, file=full_filename, row.names=FALSE)
    # gammas <- data.frame(cbind((1:length(gamma_hat)) - 1, gamma_hat))
    # names(gammas) <- c("j", "gamma_j")
    # full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "gamma.csv")
    # write.csv(gammas, file=full_filename, row.names=FALSE)



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

    # Calculate Brier score
    y0_hat <- as.numeric(exp(X_test %*% beta_hat))
    mu_hat <- as.numeric(Z_test %*% gamma_hat)
    y0_null <- rep(exp(beta_hat_null)[1], length(y0_hat))
    mu_null <- rep(gamma_hat_null[1], length(mu_hat))

    estimated_probabilities <- sapply(times_test, function(current_time) {
      FHT_parametric_survival(current_time, mu_hat, y0_hat)
    })
    brier_score_df <- brier_score_with_censoring_on_times_with_probabilities(
      times=times_test, delta=delta_test,
      estimated_probabilities_matrix=estimated_probabilities
    )

    estimated_probabilities_null <- sapply(times_test, function(current_time) {
      FHT_parametric_survival(current_time, mu_null, y0_null)
    })
    brier_null <- brier_score_with_censoring_on_times_with_probabilities(
      times=times_test, delta=delta_test,
      estimated_probabilities_matrix=estimated_probabilities_null
    )
    brier_df_both <- data.frame(
      times=brier_score_df$times,
      brier_scores_null=brier_null$brier_scores,
      brier_scores_model=brier_score_df$brier_scores
    )
    #full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "brier_data.csv")
    #write.csv(brier_df_both, file=full_filename, row.names=FALSE)

    # Plot Brier
    # full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "brier_r2_s.pdf")
    # pdf(full_filename, width=12, height=6)
    # plot(
    #   times_test_non_null_both[order(times_test_non_null_both)],
    #   brier_r2s_non_null_both[order(times_test_non_null_both)],
    #   typ='s', ylim=c(-1, 1),
    #   ylab="Brier R2", xlab="Time"
    # )
    # rug(times_test_non_null_both)
    # abline(h=0, lty=3)
    # dev.off()
    #
    #
    # loglikelihood_df <- data.frame(
    #   null_loglikelihood=test_null_loglikelihood,
    #   loglikelihood=test_loglikelihood,
    #   deviance=test_difference_of_deviance
    # )
    # full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "test_result.csv")
    # write.csv(loglikelihood_df, file=full_filename, row.names=FALSE)










    boosting_type <- "clinical"
    # Clinical
    # M_clinical <- 10
    # CV_result_clinical <- run_CV_clinical(
    #   M_clinical, K_fold_repetitions, K, X_train, Z_train, times_train, delta_train,
    #   boost_intercepts_continually=boost_intercepts_continually
    # )
    #
    # ### WRITE CV RESULT TO FILE
    # full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "loglik.csv")
    # write.csv(CV_result_clinical$CV_errors_K_loglik, file=full_filename, row.names=FALSE)
    # full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "deviance.csv")
    # write.csv(CV_result_clinical$CV_errors_K_deviance, file=full_filename, row.names=FALSE)


    ## READ CV RESULT FROM FILE
    full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "loglik.csv")
    logliks <- read.csv(full_filename)
    # full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "deviance.csv")
    # deviances <- read.csv(full_filename)

    ### POST PROCESSING AND PLOTTING
    #ylims <- c(min(apply(logliks, 2, min)), max(apply(logliks, 2, max)))
    m_stop_from_CV <- which.min(rowMeans(logliks))
    #m_stop_from_CV <- which.max(rowMeans(deviances))
    # full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "loglik.pdf")
    # pdf(full_filename, width=12, height=6)
    # plot(rowMeans(logliks), typ='l', ylim=ylims, ylab="Negative log-likelihood", xlab="Boosting iteration")
    # for (k in 1:K_fold_repetitions) {
    #   lines(logliks[, k], lty=3, col=rgb(0, 0, 0, alpha = 0.5))
    # }
    # abline(v=m_stop_from_CV, lwd=2, col='red')
    # legend(
    #   'topright',
    #   legend=c("Sum of log-lik. on test set in 5-fold CV", "Mean of log-lik. sums in 5-fold CV", "Iteration number which minimizes mean"),
    #   col=c(rgb(0, 0, 0, alpha = 0.5), 'black', 'red'),
    #   lty=c(3, 1, 1),
    #   lwd=c(1, 1, 2)
    # )
    # dev.off()

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

    # betas <- data.frame(cbind(non_null_parameters(beta_hat) - 1, beta_hat[non_null_parameters(beta_hat)]))
    # names(betas) <- c("j", "beta_j")
    # full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "beta.csv")
    # write.csv(betas, file=full_filename, row.names=FALSE)
    # gammas <- data.frame(cbind((1:length(gamma_hat)) - 1, gamma_hat))
    # names(gammas) <- c("j", "gamma_j")
    # full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "gamma.csv")
    # write.csv(gammas, file=full_filename, row.names=FALSE)
    #


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


    # Calculate Brier score
    y0_hat <- as.numeric(exp(X_test %*% beta_hat))
    mu_hat <- as.numeric(Z_test %*% gamma_hat)
    y0_null <- rep(exp(beta_hat_null)[1], length(y0_hat))
    mu_null <- rep(gamma_hat_null[1], length(mu_hat))

    estimated_probabilities <- sapply(times_test, function(current_time) {
      FHT_parametric_survival(current_time, mu_hat, y0_hat)
    })
    brier_score_df <- brier_score_with_censoring_on_times_with_probabilities(
      times=times_test, delta=delta_test,
      estimated_probabilities_matrix=estimated_probabilities
    )

    estimated_probabilities_null <- sapply(times_test, function(current_time) {
      FHT_parametric_survival(current_time, mu_null, y0_null)
    })
    brier_null <- brier_score_with_censoring_on_times_with_probabilities(
      times=times_test, delta=delta_test,
      estimated_probabilities_matrix=estimated_probabilities_null
    )
    brier_df_clinical <- data.frame(
      times=brier_score_df$times,
      brier_scores_null=brier_null$brier_scores,
      brier_scores_model=brier_score_df$brier_scores
    )

    full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "brier_data.csv")
    write.csv(brier_df_clinical, file=full_filename, row.names=FALSE)


    # Plot Brier
    # full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "brier_r2_s.pdf")
    # pdf(full_filename, width=12, height=6)
    # plot(
    #   times_test_non_null_clinical[order(times_test_non_null_clinical)],
    #   brier_r2s_non_null_clinical[order(times_test_non_null_clinical)],
    #   typ='s', ylim=c(-1, 1),
    #   ylab="Brier R2", xlab="Time"
    # )
    # rug(times_test_non_null_clinical)
    # abline(h=0, lty=3)
    # dev.off()
    #
    #
    #
    # loglikelihood_df <- data.frame(
    #   null_loglikelihood=test_null_loglikelihood,
    #   loglikelihood=test_loglikelihood,
    #   deviance=test_difference_of_deviance
    # )
    # full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "test_result.csv")
    # write.csv(loglikelihood_df, file=full_filename, row.names=FALSE)



    boosting_type <- "genetic"
    # Genetic
    # M_genetic <- 50
    # CV_result_genetic <- run_CV_genetic(
    #   M_genetic, K_fold_repetitions, K, X_train, Z_train, times_train, delta_train,
    #   boost_intercepts_continually=boost_intercepts_continually
    # )
    #
    #
    # ### WRITE CV RESULT TO FILE
    # full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "loglik.csv")
    # write.csv(CV_result_genetic$CV_errors_K_loglik, file=full_filename, row.names=FALSE)
    # full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "deviance.csv")
    # write.csv(CV_result_genetic$CV_errors_K_deviance, file=full_filename, row.names=FALSE)


    ## READ CV RESULT FROM FILE
    full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "loglik.csv")
    logliks <- read.csv(full_filename)
    # full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "deviance.csv")
    # deviances <- read.csv(full_filename)

    ### POST PROCESSING AND PLOTTING
    #logliks <- CV_result_genetic$CV_errors_K_loglik
    #ylims <- c(min(apply(logliks, 2, min)), max(apply(logliks, 2, max)))
    m_stop_from_CV <- which.min(rowMeans(logliks))
    #m_stop_from_CV <- which.max(rowMeans(deviances))
    # full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "loglik.pdf")
    # pdf(full_filename, width=12, height=6)
    # plot(rowMeans(logliks), typ='l', ylim=ylims, ylab="Negative log-likelihood", xlab="Boosting iteration")
    # for (k in 1:K_fold_repetitions) {
    #   lines(logliks[, k], lty=3, col=rgb(0, 0, 0, alpha = 0.5))
    # }
    # abline(v=m_stop_from_CV, lwd=2, col='red')
    # legend(
    #   'topright',
    #   legend=c("Sum of log-lik. on test set in 5-fold CV", "Mean of log-lik. sums in 5-fold CV", "Iteration number which minimizes mean"),
    #   col=c(rgb(0, 0, 0, alpha = 0.5), 'black', 'red'),
    #   lty=c(3, 1, 1),
    #   lwd=c(1, 1, 2)
    # )
    # dev.off()

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

    # betas <- data.frame(cbind(non_null_parameters(beta_hat) - 1, beta_hat[non_null_parameters(beta_hat)]))
    # names(betas) <- c("j", "beta_j")
    # full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "beta.csv")
    # write.csv(betas, file=full_filename, row.names=FALSE)
    # gammas <- data.frame(cbind((1:length(gamma_hat)) - 1, gamma_hat))
    # names(gammas) <- c("j", "gamma_j")
    # full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "gamma.csv")
    # write.csv(gammas, file=full_filename, row.names=FALSE)



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



    # Calculate Brier score
    y0_hat <- as.numeric(exp(X_test %*% beta_hat))
    mu_hat <- as.numeric(Z_test %*% gamma_hat)
    y0_null <- rep(exp(beta_hat_null)[1], length(y0_hat))
    mu_null <- rep(gamma_hat_null[1], length(mu_hat))

    estimated_probabilities <- sapply(times_test, function(current_time) {
      FHT_parametric_survival(current_time, mu_hat, y0_hat)
    })
    brier_score_df <- brier_score_with_censoring_on_times_with_probabilities(
      times=times_test, delta=delta_test,
      estimated_probabilities_matrix=estimated_probabilities
    )

    estimated_probabilities_null <- sapply(times_test, function(current_time) {
      FHT_parametric_survival(current_time, mu_null, y0_null)
    })
    brier_null <- brier_score_with_censoring_on_times_with_probabilities(
      times=times_test, delta=delta_test,
      estimated_probabilities_matrix=estimated_probabilities_null
    )
    brier_df_genetic <- data.frame(
      times=brier_score_df$times,
      brier_scores_null=brier_null$brier_scores,
      brier_scores_model=brier_score_df$brier_scores
    )

    full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "brier_data.csv")
    write.csv(brier_df_genetic, file=full_filename, row.names=FALSE)

    # Plot Brier
    # full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "brier_r2_s.pdf")
    # pdf(full_filename, width=12, height=6)
    # plot(
    #   times_test_non_null_genetic[order(times_test_non_null_genetic)],
    #   brier_r2s_non_null_genetic[order(times_test_non_null_genetic)],
    #   typ='s', #ylim=c(-1, 1),
    #   ylab="Brier R2", xlab="Time"
    # )
    # rug(times_test_non_null_genetic)
    # abline(h=0, lty=3)
    # dev.off()
    #
    #
    #
    # loglikelihood_df <- data.frame(
    #   null_loglikelihood=test_null_loglikelihood,
    #   loglikelihood=test_loglikelihood,
    #   deviance=test_difference_of_deviance,
    #   average_brier_r2=average_brier_r2
    # )
    # full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "test_result.csv")
    # write.csv(loglikelihood_df, file=full_filename, row.names=FALSE)



    # y_min <- min(min(brier_r2s_non_null_both), min(brier_r2s_non_null_clinical), min(brier_r2s_non_null_genetic))
    # y_max <- max(max(brier_r2s_non_null_both), max(brier_r2s_non_null_clinical), max(brier_r2s_non_null_genetic))
    # ylim <- c(y_min, y_max)
    #
    # colors <- c("black", "red", "blue")

    # full_filename <- paste0(directory, seed_string, "_", "all_briers.pdf")
    # pdf(full_filename, width=12, height=6)
    # plot(
    #   times_test_non_null_both[order(times_test_non_null_both)],
    #   brier_r2s_non_null_both[order(times_test_non_null_both)],
    #   col=colors[1],
    #   ylim=ylim, xlab="Time", ylab="Brier R2", typ="s"
    # )
    # lines(
    #   times_test_non_null_clinical[order(times_test_non_null_clinical)],
    #   brier_r2s_non_null_clinical[order(times_test_non_null_clinical)],
    #   col=colors[2],
    #   typ='s'
    # )
    # lines(
    #   times_test_non_null_genetic[order(times_test_non_null_genetic)],
    #   brier_r2s_non_null_genetic[order(times_test_non_null_genetic)],
    #   col=colors[3],
    #   typ='s'
    # )
    # legend(
    #   'bottomright',
    #   legend=c("Full model", "Clinical", "Genetic"),
    #   col=colors,
    #   lty = 1
    # )
    # abline(h=0, lty=3)
    # dev.off()

    # brier_mean_df <- data.frame(
    #   full=mean_brier_r2_both,
    #   clinical=mean_brier_r2_clinical,
    #   genetic=mean_brier_r2_genetic
    # )
    # full_filename <- paste0(directory, seed_string, "_", "mean_brier_scores.csv")
    # write.csv(brier_mean_df, file=full_filename, row.names=FALSE)
    #
    # brier_median_df <- data.frame(
    #   full=median_brier_r2_both,
    #   clinical=median_brier_r2_clinical,
    #   genetic=median_brier_r2_genetic
    # )
    # full_filename <- paste0(directory, seed_string, "_", "median_brier_scores.csv")
    # write.csv(brier_median_df, file=full_filename, row.names=FALSE)





    library(CoxBoost)
    library(pec)

    K_fold_repetitions <- 10
    boosting_type <- "cox"
    K <- 5
    repeated_cross_validation_cox <- function(seed, maxstep, time, status, xx, penalty, K, unpen.index=NULL) {
      set.seed(seed)
      cv.CoxBoost(time, status, xx, penalty=penalty, K=K, maxstepno=maxstep, unpen.index=unpen.index)$mean.logplik
    }

    MAX_STEPS <- 300 # more?

    # Formula: CoxBoost estimator equal to mboost Cox ??
    nu <- 0.1
    N <- length(times_train)
    lambda <- N*(1 - nu)/nu

    design_matrix <- as.matrix(cbind(X_train[, -1], Z_train[, -1]))
    repeated_cv_result_cox <- sapply(
      1:K_fold_repetitions, repeated_cross_validation_cox, maxstep=MAX_STEPS, time=times_train, status=delta_train,
      xx=design_matrix, penalty=lambda, K=K
    )
    mstop <- which.max(apply(repeated_cv_result_cox, 1, mean))
    ymin <- min(apply(repeated_cv_result_cox, 2, min))
    ymax <- max(apply(repeated_cv_result_cox, 2, max))
    ylim <- c(ymin, ymax)
    # plot(rowMeans(repeated_cv_result_cox), typ='l', ylim=ylim)
    # for (i in 1:K_fold_repetitions) {
    #   lines(repeated_cv_result_cox[, i], lty=3)
    # }
    # abline(v=mstop, lwd=2, col='red')
    cox_model <- CoxBoost(
      time=times_train, status=delta_train, x=design_matrix, penalty=lambda, stepno=mstop
    )


    # Predict CoxBoost
    design_matrix_test <- data.frame(X_test[, -1], Z_test[, -1])
    linear_predictors <- as.numeric(predict(
      cox_model, newdata=design_matrix_test, newtime=times_test, newstatus=delta_test, at.step=mstop, type="lp"
    ))

    # current_time
    estimate_baseline_hazard <- function(times_test, delta_test, linear_predictors) {
      N <- length(times_test)
      jumps <- rep(0, N)
      num_events <- rep(0, N)
      denominator <- rep(0, N)
      # for each timepoint t_i
      order_times <- order(times_test)
      times_sorted <- sort(times_test)
      delta_sorted <- delta_test[order_times]
      exp_lp_sorted <- exp(linear_predictors)[order_times]
      for (i in 1:N) {
        current_time <- times_sorted[i]
        at_risk_indicator <- current_time <= times_sorted
        denominator[i] <- sum(at_risk_indicator * exp_lp_sorted)
        is_event <- delta_sorted[i]
        jumps[i] <- is_event/denominator[i]
      }
      A0 <- cumsum(jumps)
      return(A0)
    }
    A0 <- estimate_baseline_hazard(times_test, delta_test, linear_predictors)
    # baseline_hazard <- exp(-A0)
    # plot(times_sorted, baseline_hazard, typ='s', ylim=c(0, 1), lwd=3)
    # for (i in 1:5) {
    #   lines(times_sorted, exp(-A0*exp_lp_sorted[i]), typ='s', col='red', lty=3)
    # }


    # BRIER SCORES STUFF
    exp_lp_sorted <- exp(linear_predictors)
    N <- length(times_test)
    estimated_probabilities_matrix_cox <- t(sapply(1:N, function(i) exp(-A0*exp_lp_sorted[i])))

    brier_score_df <- brier_score_with_censoring_on_times_with_probabilities(
      times=times_test, delta=delta_test,
      estimated_probabilities_matrix=estimated_probabilities_matrix_cox
    )

    estimated_probabilities_null_matrix_cox <- t(matrix(rep(exp(-A0), N), nrow=N))

    brier_null <- brier_score_with_censoring_on_times_with_probabilities(
      times=times_test, delta=delta_test,
      estimated_probabilities_matrix=estimated_probabilities_null_matrix_cox
    )
    brier_df_cox <- data.frame(
      times=brier_score_df$times,
      brier_scores_null=brier_null$brier_scores,
      brier_scores_model=brier_score_df$brier_scores
    )

    full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "brier_data.csv")
    write.csv(brier_df_cox, file=full_filename, row.names=FALSE)

    # plot(
    #   times_test_non_null_both[order(times_test_non_null_both)],
    #   brier_r2s_non_null_both[order(times_test_non_null_both)],
    #   typ='s', ylim=c(-1, 1),
    #   ylab="Brier R2", xlab="Time"
    # )
    # abline(h=0, lty=3)
    # plot(times_test_non_null_cox, brier_r2s_non_null_cox, col='red', typ='s')
    # abline(h=0, lty=3)
    # lines(
    #   times_test_non_null_both[order(times_test_non_null_both)],
    #   brier_r2s_non_null_both[order(times_test_non_null_both)],
    #   typ='s'
    # )




    boosting_type <- "cox_mandatory"

    vector_of_clinical_indexes <- c(9979, 9980)
    repeated_cv_result_cox_mandatory <- sapply(
      1:K_fold_repetitions, repeated_cross_validation_cox, maxstep=MAX_STEPS, time=times_train, status=delta_train,
      xx=design_matrix, penalty=lambda, K=K, unpen.index=vector_of_clinical_indexes
    )
    mstop_mandatory <- which.max(apply(repeated_cv_result_cox_mandatory, 1, mean))

    # PLOTTING
    ymin <- min(apply(repeated_cv_result_cox_mandatory, 2, min))
    ymax <- max(apply(repeated_cv_result_cox_mandatory, 2, max))
    ylim <- c(ymin, ymax)
    # plot(rowMeans(repeated_cv_result_cox_mandatory), typ='l', ylim=ylim)
    # for (i in 1:K_fold_repetitions) {
    #   lines(repeated_cv_result_cox_mandatory[, i], lty=3)
    # }
    # abline(v=mstop_mandatory, lwd=2, col='red')

    cox_model_mandatory <- CoxBoost(
      time=times_train, status=delta_train, x=design_matrix, penalty=lambda,
      stepno=mstop_mandatory,
      unpen.index=vector_of_clinical_indexes
    )
    linear_predictors_mandatory <- as.numeric(predict(
      cox_model_mandatory, newdata=design_matrix_test, newtime=times_test, newstatus=delta_test,
      at.step=mstop_mandatory, type="lp"
    ))
    A0_mandatory <- estimate_baseline_hazard(times_test, delta_test, linear_predictors_mandatory)
    exp_lp_mandatory_sorted <- exp(linear_predictors_mandatory)

    estimated_probabilities_matrix_cox_mandatory <- t(sapply(1:N, function(i) exp(-A0*exp_lp_mandatory_sorted[i])))

    brier_score_df <- brier_score_with_censoring_on_times_with_probabilities(
      times=times_test, delta=delta_test,
      estimated_probabilities_matrix=estimated_probabilities_matrix_cox_mandatory
    )

    estimated_probabilities_null_matrix_cox_mandatory <- t(matrix(rep(exp(-A0_mandatory), N), nrow=N))
    brier_null <- brier_score_with_censoring_on_times_with_probabilities(
      times=times_test, delta=delta_test,
      estimated_probabilities_matrix=estimated_probabilities_null_matrix_cox_mandatory
    )
    brier_df_cox_mandatory <- data.frame(
      times=brier_score_df$times,
      brier_scores_null=brier_null$brier_scores,
      brier_scores_model=brier_score_df$brier_scores
    )


    full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "brier_data.csv")
    write.csv(brier_df_cox_mandatory, file=full_filename, row.names=FALSE)


    # plot(brier_df_cox$times, brier_df_cox$brier_scores_model, typ='s')
    # lines(brier_df_cox_mandatory$times, brier_df_cox_mandatory$brier_scores_model, typ='s', col='red')





}
stopCluster(cl)

#
# Analyzing Brier scores
mean_briers <- c()
median_briers <- c()
for (seed in 1:100) {
  seed_string <- formatC(seed, width=3, flag="0")
  full_filename <- paste0(directory, seed_string, "_", "mean_brier_scores.csv")
  mean_briers <- rbind(mean_briers, read.csv(full_filename))
  full_filename <- paste0(directory, seed_string, "_", "median_brier_scores.csv")
  median_briers <- rbind(median_briers, read.csv(full_filename))
}

xlab <- bquote(.("Mean Brier") ~ R^2)
boxplot(
  mean_briers$full, mean_briers$genetic, mean_briers$clinical,
  xlab=xlab,
  horizontal=TRUE
)
axis(2, labels=c("Full", "Genetic", "Clinical"), at=1:3, las=2)
abline(v=0, lty=3)

xlab <- bquote(.("Median Brier") ~ R^2)
boxplot(
  median_briers$full, median_briers$genetic, median_briers$clinical,
  xlab=xlab,
  horizontal=TRUE
)
axis(2, labels=c("Full", "Genetic", "Clinical"), at=1:3, las=2)
abline(v=0, lty=3)


seed <- 30



seed <- 30
seed_string <- formatC(seed, width=3, flag="0")
boosting_type <- "cox"
full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "brier_data.csv")
brier_df_cox <- read.csv(full_filename)
plot(brier_df_cox$times, brier_df_cox$brier_scores_model, typ='s')
boosting_type <- "both"
full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "brier_data.csv")
brier_df_both <- read.csv(full_filename)
lines(brier_df_both$times, brier_df_both$brier_scores_model, typ='s', col='red')
abline(h=0, lty=3)

boosting_type <- "cox_mandatory"
full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "brier_data.csv")
brier_df_cox_mandatory <- read.csv(full_filename)
lines(brier_df_cox_mandatory$times, brier_df_cox_mandatory$brier_scores_model, typ='s', col='blue')

boosting_type <- "genetic"
full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "brier_data.csv")
brier_df_cox_mandatory <- read.csv(full_filename)
lines(brier_df_cox_mandatory$times, brier_df_cox_mandatory$brier_scores_model, typ='s', col='black')
