# Create summary plots and shit
#####################################
tex_figures_directory <- "../../text/figures/"
full_filename <- paste0(tex_figures_directory, "age_scatter.pdf")
pdf(full_filename, width=12, height=6)
plot(clinicalData[has_age_observations, 4], xlab="", ylab="Age", xaxt='n')
dev.off()

full_filename <- paste0(tex_figures_directory, "age_hist_finer.pdf")
pdf(full_filename, width=12, height=6)
hist(clinicalData[has_age_observations, 4], xlab="Age", ylab="Frequency", breaks=seq(0, 10000, by=250),
     main="", xlim=c(0, 9000))
dev.off()


directory <- "../dataset/oberthuer/oberthur_all/"
full_filename <- paste0(directory, seed_string, "_", boosting_type, "_", "loglik.csv")
logliks <- read.csv(full_filename)
ylims <- c(min(apply(logliks, 2, min)), max(apply(logliks, 2, max)))
m_stop_from_CV <- which.min(rowMeans(logliks))
full_filename <- paste0(tex_figures_directory, "example_cv_loglik.pdf")
pdf(full_filename, width=12, height=6)
plot(rowMeans(logliks), typ='l', ylim=c(80, 100), ylab="Negative log-likelihood", xlab="Boosting iteration", xlim=c(0, 100))
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


# Wiener process
y0 <- 2.00
mu <- 0.077
wiener_seed <- 6
set.seed(wiener_seed)
N <- 10000
time_increment <- 0.01
time_end <- 12
times <- seq(0, time_end, by=time_increment)
increments <- matrix(
  rnorm(
    N*length(times), mean=0, sd=sqrt(time_increment)
  ), ncol=N)
noise_matrix <- apply(increments, 2, cumsum)
y0_matrix <- matrix(rep(y0, length(times)*N), ncol=N, byrow=T)
mu_matrix <- matrix(rep(mu*times, N), ncol=N)
w <- y0_matrix + mu_matrix + noise_matrix
below_0_2 <- apply(w, 2, min) < 0
below_0 <- apply(w, 2, min) < 0
dead <- sum(below_0)
red_colors <- which(below_0)
colors <- rep(rgb(0, 0, 0, alpha = 1), N)
colors[red_colors] <- rgb(1, 0, 0, alpha = 0.6)
full_filename <- paste0(tex_figures_directory, "example_wieners", wiener_seed, ".pdf")
pdf(full_filename, width=12, height=6)
matplot(times, w, type='l', ylab="Health Y(t)", xlab="Time t", ylim=c(0, 10), xlim=c(0, time_end), col=colors)
abline(h=0, lty=3)
legend(
  'topleft',
  legend=c("Processes still alive at t = 12", "Processes dead by t = 12"),
  col=c("Black", rgb(1, 0, 0, alpha = 0.6)),
  lty=c(1, 1)
)
dev.off()




### PLOT DEVIANCE DIFF
directory <- "../dataset/oberthuer/oberthur_all/"
boosting_types <- c("both", "clinical", "genetic", "cox", "cox_mandatory")
get_deviance <- function(boosting_type) {
  M <- 100
  deviances <- rep(NA, M)
  for (seed in 1:M) {
    seed_string <- formatC(seed, width=3, flag="0")
    df <- read.csv(paste0(directory, seed_string, '_', boosting_type, '_', "test_result.csv"))
    deviances[seed] <- df$deviance
  }
  return(deviances)
}
deviance_df <- data.frame(
  "both"=get_deviance("both"),
  "clinical"=get_deviance("clinical"),
  "genetic"=get_deviance("genetic")#,
  #"cox"=get_deviance("cox"),
  #"cox_mandatory"=get_deviance("cox_mandatory")
)
full_filename <- paste0(tex_figures_directory, "deviance_FHT.pdf")
pdf(full_filename, width=12, height=6)
boxplot(
  deviance_df,
  xlab="Difference of deviance",
  horizontal=TRUE,
  yaxt='n'
)
axis(2, labels=c("Full", "Clinical", "Genetic"), at=1:3, las=2)
abline(v=0, lty=3)
dev.off()


# Histograms
tex_figures_directory <- "../../text/figures/"
full_filename <- paste0(tex_figures_directory, "deviances_histogram.pdf")
pdf(full_filename, width=12, height=6)
hist(deviance_df$both, breaks=seq(-200, 50, by=5), col=rgb(0.7,0.7,0.7,1/3), main='', xlab="Difference of deviance")
hist(deviance_df$genetic, breaks=seq(-200, 50, by=5), col=rgb(0.1,0.1,0.1,1/3), add=TRUE)
abline(v=median(deviance_df$both), col='red', lwd=3)
abline(v=median(deviance_df$genetic), col='firebrick', lwd=3)
abline(v=0, lty=3)
legend(
  'topleft',
  legend=c("Full model", "Genomic model", "Median of full model", "Median of genomic model"),
  col=c(rgb(0.7,0.7,0.7,1/3), rgb(0.1,0.1,0.1,1/3), "Red", "Firebrick"),
  lty=c(1, 1, 1, 1),
  lwd=c(10, 10, 3, 3)
)
dev.off()


# Plot Brier scores
directory <- "../dataset/oberthuer/oberthur_all/"
boosting_types <- c("both", "clinical", "genetic", "cox", "cox_mandatory")
get_brier_scores <- function(seed, boosting_type) {
  seed_string <- formatC(seed, width=3, flag="0")
  brier_df <- read.csv(paste0(directory, seed_string, '_', boosting_type, '_', "brier_data.csv"))
  return(brier_df)
}

seed <- 30
both_brier <- get_brier_scores(seed, "both")
clinical_brier <- get_brier_scores(seed, "clinical")
genetic_brier <- get_brier_scores(seed, "genetic")
cox_brier <- get_brier_scores(seed, "cox")
cox_m_brier <- get_brier_scores(seed, "cox_mandatory")

# merge all brier to one df
# brier_df
brier_df <- data.frame(
  times=both_brier$times,
  both=both_brier$brier_scores_model,
  clinical=clinical_brier$brier_scores_model,
  genetic=genetic_brier$brier_scores_model,
  cox=cox_brier$brier_scores_model,
  cox_m=cox_m_brier$brier_scores_model,
  FHT_null=both_brier$brier_scores_null
)
plot(brier_df$times, brier_df$both, typ='s', ylim=c(0, 0.3))
abline(h=0, lty=3)
lines(brier_df$times, brier_df$FHT_null, col='red', typ='s')
lines(brier_df$times, brier_df$cox, col='blue', typ='s')

numerical_brier_integrated <- function(times, briers) {
  time_steps <- diff(times)
  N <- length(briers)
  brier_score_until_last <- briers[1:(N-1)]
  sum(time_steps * brier_score_until_last)
}

numerical_brier_integrated_div_time <- function(times, briers) {
  time_steps <- diff(times)
  N <- length(briers)
  brier_score_until_last <- briers[1:(N-1)]
  time_length <- times[length(times)] - times[1]
  sum(time_steps * brier_score_until_last) / time_length
}

all_integrated_brier <- c()
all_integrated_brier_div_time <- c()
for (seed in 1:100) {
  both_brier <- get_brier_scores(seed, "both")
  clinical_brier <- get_brier_scores(seed, "clinical")
  genetic_brier <- get_brier_scores(seed, "genetic")
  cox_brier <- get_brier_scores(seed, "cox")
  cox_m_brier <- get_brier_scores(seed, "cox_mandatory")

  # merge all brier to one df
  # brier_df
  brier_df <- data.frame(
    times=both_brier$times,
    both=both_brier$brier_scores_model,
    clinical=clinical_brier$brier_scores_model,
    genetic=genetic_brier$brier_scores_model,
    cox=cox_brier$brier_scores_model,
    cox_m=cox_m_brier$brier_scores_model,
    FHT_null=both_brier$brier_scores_null
  )

  # start_index <- min(which(both_brier$times > 1))
  # end_index <- max(which(both_brier$times < 8))
  # indices <- start_index:end_index
  indices <- 1:length(brier_df$times)

  integrated_brier_score_both <- numerical_brier_integrated(brier_df$times[indices], brier_df$both[indices])
  integrated_brier_score_clinical <- numerical_brier_integrated(brier_df$times[indices], brier_df$clinical[indices])
  integrated_brier_score_genetic <- numerical_brier_integrated(brier_df$times[indices], brier_df$genetic[indices])
  integrated_brier_score_cox <- numerical_brier_integrated(brier_df$times[indices], brier_df$cox[indices])
  integrated_brier_score_cox_mandatory <- numerical_brier_integrated(brier_df$times[indices], brier_df$cox_m[indices])
  integrated_brier_score_FHT_null <- numerical_brier_integrated(brier_df$times[indices], brier_df$FHT_null[indices])
  integrated_brier_scores <- data.frame(
    "both"=integrated_brier_score_both,
    "clinical"=integrated_brier_score_clinical,
    "genetic"=integrated_brier_score_genetic,
    "cox"=integrated_brier_score_cox,
    "cox_mandatory"=integrated_brier_score_cox_mandatory,
    "FHT_null"=integrated_brier_score_FHT_null
  )
  all_integrated_brier <- rbind(all_integrated_brier, integrated_brier_scores)

  # Write to file
  seed_string <- formatC(seed, width=2, flag="0")
  full_filename <- paste0(directory, seed_string, '_', "integrated_brier_scores.csv")
  write.csv(integrated_brier_scores, full_filename, row.names = FALSE)

  integrated_brier_score_both <- numerical_brier_integrated_div_time(brier_df$times[indices], brier_df$both[indices])
  integrated_brier_score_clinical <- numerical_brier_integrated_div_time(brier_df$times[indices], brier_df$clinical[indices])
  integrated_brier_score_genetic <- numerical_brier_integrated_div_time(brier_df$times[indices], brier_df$genetic[indices])
  integrated_brier_score_cox <- numerical_brier_integrated_div_time(brier_df$times[indices], brier_df$cox[indices])
  integrated_brier_score_cox_mandatory <- numerical_brier_integrated_div_time(brier_df$times[indices], brier_df$cox_m[indices])
  integrated_brier_score_FHT_null <- numerical_brier_integrated_div_time(brier_df$times[indices], brier_df$FHT_null[indices])
  integrated_brier_scores_div_time <- data.frame(
    "both"=integrated_brier_score_both,
    "clinical"=integrated_brier_score_clinical,
    "genetic"=integrated_brier_score_genetic,
    "cox"=integrated_brier_score_cox,
    "cox_mandatory"=integrated_brier_score_cox_mandatory,
    "FHT_null"=integrated_brier_score_FHT_null
  )
  all_integrated_brier_div_time <- rbind(all_integrated_brier_div_time, integrated_brier_scores_div_time)
  full_filename <- paste0(directory, seed_string, '_', "integrated_brier_scores_div_time.csv")
  write.csv(integrated_brier_scores_div_time, full_filename, row.names = FALSE)
}
full_filename <- paste0(directory, "all_integrated_brier_scores.csv")
write.csv(all_integrated_brier, full_filename, row.names=FALSE)
full_filename <- paste0(directory, "all_integrated_brier_scores_div_time.csv")
write.csv(all_integrated_brier_div_time, full_filename, row.names=FALSE)


full_filename <- paste0(tex_figures_directory, "integrated_brier_boxplot.pdf")
pdf(full_filename, width=12, height=6)
par(oma=c(0,4,0,0))
boxplot(
  all_integrated_brier_div_time,
  xlab="Integrated Brier scores",
  horizontal=TRUE,
  yaxt='n',
  ylim=c(0, 0.3)
)
labels <- c("FHT (Full)", "FHT (Clinical)", "FHT (Genomic)", "CoxBoost", "CoxBoost (mand.)", "FHT (Null)")
axis(2, labels=labels, at=1:length(labels), las=2)
abline(v=0, lty=3)
par(oma=c(0,0,0,0))
dev.off()




# BRIER: Both vs Cox
full_filename <- paste0(tex_figures_directory, "brier_cox_both.pdf")
pdf(full_filename, width=12, height=6)
par(oma=c(0,1.5,0,0))
plot(cox_brier$times, cox_brier$brier_scores_model, typ='s', xlab="Time", ylab="", col="black", ylim=c(0, 0.3), las=1, yaxt="n")
mtext(expression("Brier score"),side=2,las=1,line=1)
axis(2, labels=c("0.0", "0.1", "0.2", "0.3"), at=c(0.0, 0.1, 0.2, 0.3), las=2)
lines(both_brier$times, both_brier$brier_scores_model, typ='s', col='red', lty=5)
#lines(genetic_brier$times, genetic_brier$brier_scores_model, typ='s', col='blue')
abline(h=0, lty=3)
legend(
  'topleft',
  legend=c("CoxBoost model","Full FHT model"),
  lty=c(1, 5),
  col=c("Black", "Red")
)
dev.off()
par(oma=c(0,0,0,0))

# BRIER FHT
full_filename <- paste0(tex_figures_directory, "brier_FHT.pdf")
pdf(full_filename, width=12, height=6)
par(oma=c(0,1.5,0,0))
plot(both_brier$times, both_brier$brier_scores_model, typ='s', xlab="Time", ylab="", lty=5,
     ylim=c(0, 0.3), las=1, yaxt="n")
mtext(expression("Brier score"),side=2,las=1,line=1)
axis(2, labels=c("0.0", "0.1", "0.2", "0.3"), at=c(0.0, 0.1, 0.2, 0.3), las=2)
lines(clinical_brier$times, clinical_brier$brier_scores_model, typ='s', col='red', lty=1)
lines(genetic_brier$times, genetic_brier$brier_scores_model, typ='s', col='blue', lty=3)
abline(h=0, lty=3)
legend(
  'topleft',
  legend=c("Full FHT model", "Clinical FHT model", "Genetic FHT model"),
  lty=c(5, 1, 3),
  col=c("Black", "Red", "Blue")
)
dev.off()
par(oma=c(0,0,0,0))

# BRIER: COX AND GENETIC
full_filename <- paste0(tex_figures_directory, "brier_cox_genetic.pdf")
pdf(full_filename, width=12, height=6)
par(oma=c(0,1.5,0,0))
plot(cox_brier$times, cox_brier$brier_scores_model, typ='s', xlab="Time", ylab="", lty=1, las=1, ylim=c(0, 0.3), yaxt="n")
mtext(expression("Brier score"),side=2,las=1,line=1)
axis(2, labels=c("0.0", "0.1", "0.2", "0.3"), at=c(0.0, 0.1, 0.2, 0.3), las=2)
lines(genetic_brier$times, genetic_brier$brier_scores_model, typ='s', col='red', lty=5)
abline(h=0, lty=3)
legend(
  'topleft',
  legend=c("CoxBoost model", "Genetic FHT model"),
  lty=c(1, 5),
  col=c("Black", "Red")
)
par(oma=c(0,0,0,0))
dev.off()



# BRIER: COX MANDTORY
full_filename <- paste0(tex_figures_directory, "brier_cox_mandatory.pdf")
pdf(full_filename, width=12, height=6)
par(oma=c(0,1.5,0,0))
plot(cox_brier$times, cox_brier$brier_scores_model, typ='s', xlab="Time", ylab="", lty=1, ylim=c(0, 0.3), las=1, yaxt='n')
mtext(expression("Brier score"),side=2,las=1,line=1)
axis(2, labels=c("0.0", "0.1", "0.2", "0.3"), at=c(0.0, 0.1, 0.2, 0.3), las=2)
lines(cox_m_brier$times, cox_m_brier$brier_scores_model, typ='s', col="red", lty=5)
abline(h=0, lty=3)
legend(
  'topleft',
  legend=c("CoxBoost model", "CoxBoost model with mandatory clinical covariates"),
  lty=c(1, 5),
  col=c("Black", "Red")
)
par(oma=c(0,0,0,0))
dev.off()



numerical_brier_integrated <- function(times, briers) {
  time_steps <- diff(times)
  N <- length(briers)
  brier_score_until_last <- briers[1:(N-1)]
  sum(time_steps * brier_score_until_last)
}


# BRIER: ALL
full_filename <- paste0(tex_figures_directory, "brier_all.pdf")
pdf(full_filename, width=12, height=6)
plot(cox_brier$times, cox_brier$brier_scores_model, typ='s', xlab="Time", ylab="Brier score", lty=1)
lines(both_brier$times, both_brier$brier_scores_model, typ='s', lty=6, col="Black")
lines(genetic_brier$times, genetic_brier$brier_scores_model, typ='s', col='red', lty=5)
lines(clinical_brier$times, clinical_brier$brier_scores_model, typ='s', col='steelblue', lty=4)
lines(cox_m_brier$times, cox_m_brier$brier_scores_model, typ='s', col="black")
abline(h=0, lty=3)
legend(
  'topleft',
  legend=c("CoxBoost model", "Full FHT model", "Clinical FHT model", "Genetic FHT model"),
  lty=c(1, 6, 5, 4),
  col=c("Black", "Black", "Red", "Steelblue")
)
dev.off()




full_filename <- paste0(tex_figures_directory, "deviance_FHT.pdf")
pdf(full_filename, width=12, height=6)
boxplot(
  deviance_df,
  xlab="Difference of deviance",
  horizontal=TRUE,
  yaxt='n'
)
axis(2, labels=c("Full", "Clinical", "Genetic"), at=1:3, las=2)
abline(v=0, lty=3)
dev.off()


# get deviance df
get_deviance_on_seed <- function(seed, boosting_type) {
  seed_string <- formatC(seed, width=3, flag="0")
  df <- read.csv(paste0(directory, seed_string, '_', boosting_type, '_', "test_result.csv"))
  return(df$deviance)
}
seed <- 31
deviance_seed_df <- data.frame(
  "both"=get_deviance_on_seed(seed, "both"),
  "clinical"=get_deviance_on_seed(seed, "clinical"),
  "genetic"=get_deviance_on_seed(seed, "genetic")#,
  #"cox"=get_deviance_on_seed(30, "cox"),
  #"cox_mandatory"=get_deviance_on_seed(30, "cox_mandatory")
)















#########################
# NOT CORRELATED DEVIANCES
#
seeds <- 1:500
not_seeds <- c(15, 17, 18, 24, 30, 57, 58, 63, 68, 70, 78, 79, 90, 93, 95, 103, 110, 113, 123, 140, 141, 144, 149, 160, 170, 171, 179, 181, 183, 184, 189, 193, 200, 201, 206, 213, 217, 218, 236, 238, 241, 244, 247, 248, 249, 250, 251, 256, 276, 278, 308, 314, 316, 317, 325, 337, 338, 344, 349, 361, 366, 376, 381, 382, 391, 394, 396, 404, 421, 427, 428, 435, 459, 461, 463, 464, 467, 472, 490, 491, 492, 493, 494, 500)
seeds <- seeds[-not_seeds]
initial_seed <- seeds[1]
base_input_directory <- "../dataset/scenario-"
scenario <- 'non-correlated'
algorithm <- 'non-cyclic-no-intercept'
input_directory_no <- paste0(base_input_directory, scenario, '/', algorithm, '/')
algorithm <- 'non-cyclic-intercept'
input_directory_yes <- paste0(base_input_directory, scenario, '/', algorithm, '/')

deviances_no_intercept <- c()
deviances_yes_intercept <- c()
for (seed in seeds) {
  summary <- read.csv(make_filename(input_directory_no, "summary", seed))
  deviances_no_intercept <- c(deviances_no_intercept, -summary$deviance)
  summary <- read.csv(make_filename(input_directory_yes, "summary", seed))
  deviances_yes_intercept <- c(deviances_yes_intercept, -summary$deviance)
}

full_filename <- paste0(tex_figures_directory, "deviances_simulation_not_correlated.pdf")
pdf(full_filename, width=12, height=6)
boxplot(
  deviances_yes_intercept, deviances_no_intercept,
  xlab="Difference of deviance",
  horizontal=TRUE,
  yaxt='n'
)
axis(2, labels=c("With", "Without"), at=1:2, las=2)
abline(v=0, lty=3)
dev.off()




#########################
# CORRELATED DEVIANCES
#
seeds <- 1:500
not_seeds <- c(1, 3, 4, 6, 11, 23, 25, 26, 33, 35, 36, 37, 40, 41, 42, 45, 52, 54, 57, 61, 62, 63, 64, 65, 66, 68, 70, 72, 80, 81, 84, 87, 92, 95, 101, 107, 109, 110, 111, 114, 120, 121, 122, 127, 130, 133, 136, 149, 152, 155, 156, 158, 160, 168, 174, 180, 182, 193, 195, 196, 210, 211, 212, 217, 221, 223, 225, 233, 234, 242, 247, 252, 253, 257, 258, 259, 277, 278, 282, 285, 288, 293, 302, 313, 315, 316, 317, 319, 320, 334, 337, 338, 339, 340, 343, 347, 354, 355, 356, 357, 358, 361, 362, 365, 372, 381, 384, 389, 392, 393, 411, 412, 415, 421, 423, 428, 430, 433, 435, 439, 440, 442, 444, 445, 447, 456, 462, 466, 468, 472, 477, 479, 483, 485, 488, 493)
seeds <- seeds[-not_seeds]

initial_seed <- seeds[1]
base_input_directory <- "../dataset/scenario-correlated/"
algorithm <- 'non-cyclic-no-intercept'
input_directory_no <- paste0(base_input_directory, algorithm, '/')
algorithm <- 'non-cyclic-intercept'
input_directory_yes <- paste0(base_input_directory, algorithm, '/')
deviances_no_intercept <- c()
deviances_yes_intercept <- c()
for (seed in seeds) {
  summary <- read.csv(make_filename(input_directory_no, "summary", seed))
  deviances_no_intercept <- c(deviances_no_intercept, -summary$deviance)
  summary <- read.csv(make_filename(input_directory_yes, "summary", seed))
  deviances_yes_intercept <- c(deviances_yes_intercept, -summary$deviance)
}

full_filename <- paste0(tex_figures_directory, "deviances_simulation_correlated.pdf")
pdf(full_filename, width=12, height=6)
boxplot(
  deviances_yes_intercept, deviances_no_intercept,
  xlab="Difference of deviance",
  horizontal=TRUE,
  yaxt='n'
)
axis(2, labels=c("With", "Without"), at=1:2, las=2)
abline(v=0, lty=3)
dev.off()







directory <- "../dataset/oberthuer/oberthur_all/"
seeds <- 1:100
boosting_types <- c("both", "clinical", "genetic", "cox", "cox_mandatory")
for (seed in seeds) {
  seed_string <- formatC(seed, width=3, flag="0")
  boosting_filename <- function(boosting_type) {
    return(paste0(directory, seed_string, '_', boosting_type, '_', "brier_data.csv"))
  }
  both <- read.csv(boosting_filename("both"))
  clinical <- read.csv(boosting_filename("clinical"))
  genetic <- read.csv(boosting_filename("genetic"))
  cox <- read.csv(boosting_filename("cox"))
  cox_mandatory <- read.csv(boosting_filename("cox_mandatory"))

  lower_limit <- 5
  start_time <- max(
    min(both$times[(both$times > lower_limit)]),
    min(clinical$times[(clinical$times > lower_limit)]),
    min(genetic$times[(genetic$times > lower_limit)]),
    min(cox$times[(cox$times > lower_limit)]),
    min(cox_mandatory$times[(cox_mandatory$times > lower_limit)])
  )

  upper_limit <- 10
  end_time <- min(
    max(both$times[(both$times < upper_limit)]),
    max(clinical$times[(clinical$times < upper_limit)]),
    max(genetic$times[(genetic$times < upper_limit)]),
    max(cox$times[(cox$times < upper_limit)]),
    max(cox_mandatory$times[(cox_mandatory$times < upper_limit)])
  )

  start_index <- which(both$times == start_time)
  end_index <- which(both$times == end_time)
  integrated_brier_score_both <- numerical_brier_integrated(
    both$times[start_index:end_index],
    both$brier_scores_model[start_index:end_index]
  )

  start_index_clinical <- which(clinical$times == start_time)
  end_index_clinical <- which(clinical$times == end_time)
  integrated_brier_score_clinical <- numerical_brier_integrated(
    clinical$times[start_index:end_index],
    clinical$brier_scores_model[start_index:end_index]
  )

  start_index_genetic <- which(genetic$times == start_time)
  end_index_genetic <- which(genetic$times == end_time)
  integrated_brier_score_genetic <- numerical_brier_integrated(
    genetic$times[start_index:end_index],
    genetic$brier_scores_model[start_index:end_index]
  )


  start_index <- which(cox$times == start_time)
  end_index <- which(cox$times == end_time)
  integrated_brier_score_cox <- numerical_brier_integrated(
    cox$times[start_index:end_index],
    cox$brier_scores_model[start_index:end_index]
  )

  start_index <- which(cox_mandatory$times == start_time)
  end_index <- which(cox_mandatory$times == end_time)
  integrated_brier_score_cox_mandatory <- numerical_brier_integrated(
    cox_mandatory$times[start_index:end_index],
    cox_mandatory$brier_scores_model[start_index:end_index]
  )

  integrated_brier_scores <- data.frame(
    "both"=integrated_brier_score_both,
    "clinical"=integrated_brier_score_clinical,
    "genetic"=integrated_brier_score_genetic,
    "cox"=integrated_brier_score_cox,
    "cox_mandatory"=integrated_brier_score_cox_mandatory
  )

  full_filename <- paste0(directory, seed_string, '_', "integrated_brier_scores_5_10.csv")
  write.csv(integrated_brier_scores, full_filename, row.names = FALSE)
}


directory <- "../dataset/oberthuer/oberthur_all/"
seeds <- 1:100
integrated_brier_scores <- c()
for (seed in seeds) {
  seed_string <- formatC(seed, width=3, flag="0")
  full_filename <- paste0(directory, seed_string, '_', "integrated_brier_scores.csv")
  df <- read.csv(full_filename)
  integrated_brier_scores <- rbind(integrated_brier_scores, as.numeric(df))
}
integrated_brier_scores <- data.frame(integrated_brier_scores)
colnames(integrated_brier_scores) <- c("both","clinical","genetic","cox","cox_mandatory")

median_integrated_brier_scores <- data.frame(apply(integrated_brier_scores, 2, median))
full_filename <- paste0(directory, "median_integrated_brier_scores.csv")
write.csv(median_integrated_brier_scores, full_filename, row.names = c("Both", "Clinical", "Genetic", "Cox", "Cox mandatory"))

full_filename <- paste0(tex_figures_directory, "integrated_brier_boxplot.pdf")
#pdf(full_filename, width=12, height=6)
par(oma=c(0,4,0,0))
#boxplot(integrated_brier_scores, ylim=c(0, 3), horizontal=T, yaxt='n')
boxplot(integrated_brier_scores, ylim=c(0, max(apply(integrated_brier_scores, 2, max))), horizontal=T, yaxt='n')
abline(v=0, lty=2)
axis(2, labels=c("Full", "Clinical", "Genetic", "Cox", "Cox (mandatory)"), at=1:5, las=2)
par(oma=c(0,0,0,0))
#dev.off()
