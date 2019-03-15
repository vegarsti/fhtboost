#################
# read the data #
#################
directory <- 'neuroblastoma/'
filename <- paste0(directory, 'webclinical.txt')
tmp <- read.table(filename, header = TRUE, sep = '\t') # data collected from the web
tmp <- tmp[seq(1, 751, by = 2), ]
tmp <- tmp[ , -c(2, 3, 4, 6, 7)]
tmp[ , 1] <- substr(tmp[ , 1], 1, 4)
tmp[ , 3] <- substr(tmp[ , 3], 1, 4)
tmp[ , 4] <- substr(tmp[ , 4], 1, 4)
tmp[ , 5] <- substr(tmp[ , 5], 1, 5)
for (i in 1:376)
{
  if (is.na(suppressWarnings(as.numeric(tmp[i, 1])))) tmp[i, 1] <- substr(tmp[i, 1], 1, 3)
  if (is.na(suppressWarnings(as.numeric(tmp[i, 1])))) tmp[i, 1] <- substr(tmp[i, 1], 1, 2)
  if (is.na(suppressWarnings(as.numeric(tmp[i, 1])))) tmp[i, 1] <- substr(tmp[i, 1], 1, 1)
  if (is.na(suppressWarnings(as.numeric(tmp[i, 3])))) tmp[i, 3] <- substr(tmp[i, 3], 1, 3)
  if (is.na(suppressWarnings(as.numeric(tmp[i, 3])))) tmp[i, 3] <- substr(tmp[i, 3], 1, 2)
  if (is.na(suppressWarnings(as.numeric(tmp[i, 3])))) tmp[i, 3] <- substr(tmp[i, 3], 1, 1)
  if (is.na(suppressWarnings(as.numeric(tmp[i, 4])))) tmp[i, 4] <- substr(tmp[i, 4], 1, 3)
  if (is.na(suppressWarnings(as.numeric(tmp[i, 4])))) tmp[i, 4] <- substr(tmp[i, 4], 1, 2)
  if (is.na(suppressWarnings(as.numeric(tmp[i, 4])))) tmp[i, 4] <- substr(tmp[i, 4], 1, 1)
  ifelse(tmp[i, 5] == 'alive', tmp[i, 5] <- 0, tmp[i, 5] <- 1)
}
tmp[ , 1] <- as.numeric(tmp[ , 1])
tmp[ , 3] <- as.numeric(tmp[ , 3])
tmp[ , 4] <- as.numeric(tmp[ , 4])
medianAge <- median(tmp$Age)

filename <- paste0(directory, 'NB2004_clin_and_genes_neuroblastoma_noTies.txt')
boev <- read.table(filename, header = TRUE, sep = '\t')
# data provided by H. M. Boevelstad
clinicalData <- cbind(boev[ , 1:3], rep(NA, dim(boev)[1]))
supp <- tmp[c(288:329, 331, 332, 335, 338, 339, 341, 343:353, 355:369, 371:376), 4]
clinicalData[283:362, 4] <- as.numeric(supp > medianAge) # the last observations are in order
for (i in 1:282)
  for (j in 1:287)
    if ((trunc(boev[i, 1] * 365) == tmp[j, 4]) & boev[i, 2] == tmp[j, 5])
      clinicalData[i, 4] <- as.numeric(tmp[j, 1] > medianAge)
# survival time and censoring status are unique key. Only possible issue: observation 247.
colnames(clinicalData) <- c('time', 'status', 'risk', 'age')

molecularData <- boev[, -c(1:3)]
colnames(molecularData) <- paste('X', 1:9978, sep = '')

rm(tmp, boev, i, j, medianAge, supp)


### HERE MY WORK STARTS (data prep)

X <- molecularData
Z <- clinicalData[, c(3, 4)]
times <- clinicalData$time
delta <- clinicalData$status

# Divide into train and test
set.seed(1)
# N is 362
N <- dim(molecularData)[1]
N_test <- 70 # slightly above 25%
test_indices <- sample(1:N, N_test)
times_test <- times[test_indices]
times_train <- times[-test_indices]
delta_test <- delta[test_indices]
delta_train <- delta[-test_indices]
X_test <- X[test_indices, ]
Z_test <- Z[test_indices, ]
X_train <- X[-test_indices, ]
Z_train <- Z[-test_indices, ]
