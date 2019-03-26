#install.packages("BiocManager")
library(BiocManager)
#BiocManager::install("ArrayExpress")
library(ArrayExpress)
system('mkdir E-GEOD-39582')
tmp <- getAE('E-GEOD-39582', type = 'processed', path = 'E-GEOD-39582')
tmpProc <- getcolproc(tmp)
tmpE <- procset(tmp, tmpProc[2])

# extract molecular data
gene <- t(exprs(tmpE))
# extract clinical data
info <- pData(tmpE)
if(!identical(rownames(info), rownames(gene))) stop('Please check the data!!!')

time <- info[ , 20]
status <- info[ , 21]
clin <- info[ , c(4, 7, 10, 23)]
colnames(clin) <- c('sex', 'age', 'subtype', 'stage')
# separate training and test sets
label <- info[ , 13] == 'discovery'

# check for possible missing values and remove the observations with missing data
nona <- apply(is.na(cbind(time, status, clin)), 1, sum) == 0
gene <- gene[nona, ]
time <- time[nona]
status <- status[nona]
clin <- clin[nona, ]
label <- label[nona]
# in case of time equal to 0, add an arbitrarily small value to avoid computational issues
# time[time == 0] <- 1e-9
# time[time == 0] <- 0.1
# time[time == 0] <- 0.01 # causes differences in parameters!!
# # time[time == 0] <- 1
zero_time_indices <- which(time == 0)
time[zero_time_indices] <- 0.1
# transform sex in a dummy variable
clin[clin[ , 1] == 'female', 1] <- -1
clin[clin[ , 1] == 'male', 1] <- 1
clin[ ,1] <- as.numeric(clin[ , 1])

# codify sub with dummy variables
sub2 <- sub3 <- sub4 <- sub5 <- sub6 <- rep(-1, length(time))
sub2[clin[ , 3] == 'C2'] <- 1
sub3[clin[ , 3] == 'C3'] <- 1
sub4[clin[ , 3] == 'C4'] <- 1
sub5[clin[ , 3] == 'C5'] <- 1
sub6[clin[ , 3] == 'C6'] <- 1

# codify stage with dummy variables
stage2 <- stage3 <- stage4 <- rep(-1, length(time))
stage2[clin[ , 4] == '2'] <- 1
stage3[clin[ , 4] == '3'] <- 1
stage4[clin[ , 4] == '4'] <- 1

# create the matrix with the clinical variables
clin <- cbind(clin[ , 1:2], sub2, sub3, sub4, sub5, sub6, stage2, stage3, stage4)
# to avoid probelms with names starting with a number
colnames(gene) <- paste0('G', colnames(gene))

# remove the working quantities
rm(tmp, tmpE, tmpNames, tmpProc, info, nona, stage2, stage3, stage4,
   sub2, sub3, sub4, sub5, sub6, sdAge, sdGene)
system('rm -r E-GEOD-39582')
#marisa_filename <- 'preproc_Marisa_data.Rdata'
#save(time, status, clin, gene, label, file=marisa_filename)


marisa_filename <- 'preproc_Marisa_data.Rdata'
load(marisa_filename)


### HERE STARTS MY PART
test_indices <- which(!label)
train_indices <- which(label)

early_times <- which(time < 0.1)
# remove "NA" age observations
indices_to_remove_train <- c()
for (i in 1:length(train_indices)) {
  if (any(train_indices[i] == early_times)) {
    indices_to_remove_train <- c(indices_to_remove_train, i)
  }
}
train_indices <- train_indices[-indices_to_remove_train]


# TRAINING SET
ones <- rep(1, length(train_indices))
X_train <- cbind(ones, scale(gene)[train_indices, ])
Z_train <- cbind(ones, scale(clin)[train_indices, ])
times_train <- time[train_indices]
delta_train <- status[train_indices]

nlm_result <- maximum_likelihood_intercepts(times_train, delta_train)

# Just looking at observed
# X_train <- X_train[delta_train == 1, ]
# Z_train <- Z_train[delta_train == 1, ]
# times_train <- times_train[delta_train == 1]
# delta_train <- delta_train[delta_train == 1]


### INITIAL ANALYSIS
library(devtools)
load_all()
non_para <- non_parametric_estimates(times_train, delta_train, continuous = TRUE)
plot(non_para$times_sequence, non_para$kaplan_meiers, typ='s', xlab="Time", ylab="Kaplan-Meier estimated survival probability")

parametric_times <- seq(0.1, max(times_train), by=0.01)
nlm_result <- maximum_likelihood_intercepts(times_train, delta_train)
cat(nlm_result)
y0 <- exp(nlm_result[1])
mu <- nlm_result[2]
parametric_S <- FHT_parametric_survival(parametric_times, mu, y0)
lines(parametric_times, parametric_S, col='red')

beta_null <- 1
mus <- seq(-2, 0, by=0.01)
logliks <- rep(0, length(mus))

for (i in 1:length(mus)) {
  mu <- mus[i]
  parameters <- c(beta_null, mu)
  logliks[i] <- FHT_only_intercepts(parameters, times_train, delta_train)
}
plot(mus, logliks, typ='l', main=beta_null, ylim=c(0, 50000))
cat(beta_null, min(logliks))

### RUN FHTBOOST
load_all()

M <- 30 # ??
K_fold_repetitions <- 10
K <- 10
boost_intercepts_continually <- FALSE

result2 <- boosting_run(
  times=times_train,
  delta=delta_train,
  X=X_train,
  Z=Z_train,
  m_stop=M,
  boost_intercepts_continually=FALSE,
  should_print=TRUE
)
plot(result2$loss, typ='l')

#nlm_result <- maximum_likelihood_intercepts(times_train, delta_train)
beta_null <- 0.0
gamma_null <- -0.1
parameters <- c(beta_null, gamma_null)

y0 <- exp(beta_null)
mu <- gamma_null
sigma2 <- 1

vector <- FHT_loglikelihood_with_y0_mu(y0, mu, times_train, delta_train)

y0_deriv <- loss_function_derivative_y0(y0, mu, sigma2, times_train, delta_train)
mu_deriv <- loss_function_derivative_mu(y0, mu, sigma2, times_train, delta_train)






########### COXBOOST
nCV <- 5 # number of repetitions for cross validation
library("CoxBoost")
lambda<-(sum(label)-1)*9
rep_cv.lb <- function(seed, maxstep, time, status, xx, penalty, unpen.index=NULL) {
  set.seed(seed)
  return(cv.CoxBoost(time, status, xx, penalty=penalty, maxstepno=maxstep, unpen.index=unpen.index)$mean.logplik)
}

par2<-which.max(apply(sapply(1:nCV,rep_cv.lb,maxstep=100,time=time[label],status=status[label],
                               xx=as.matrix(cbind(clin,gene)[label,]),penalty=lambda),1,mean))
ctime$cb<-system.time(mod2<-CoxBoost(time[label],status[label],as.matrix(cbind(clin,gene)[label,]),penalty=lambda,stepno=par2))
data2<-data.frame(time,status,cbind(clin,gene)[,mod2$coefficients[par2+1,]!=0])[label,]
if(ncol(data2)==3) colnames(data2)[3]<-colnames(cbind(clin,gene))[mod2$coefficients[par2+1,]!=0]
models$lb.naive<-coxph(Surv(time,status)~.,data=data2,init=mod2$coefficients[par2+1,mod2$coefficients[par2+1,]!=0],iter.max=0)
