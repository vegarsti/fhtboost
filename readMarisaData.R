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
time[time == 0] <- 1e-9

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


### HERE STARTS MY PART
test_indices <- which(!label)
train_indices <- which(label)
n_train <- length(train_indices)

# TRAINING SET
ones <- rep(1, n_train)
gene_standardized <- scale(gene)
clin_standardized <- scale(clin)
X_train <- cbind(ones, gene_standardized[train_indices, ])
Z_train <- cbind(ones, clin_standardized[train_indices, ])
times_train <- time[train_indices]
delta_train <- status[train_indices]


### RUN FHTBOOST
library(devtools)
load_all()

M <- 10 # ??
K_fold_repetitions <- 10
K <- 10
boost_intercepts_continually <- FALSE

result <- boosting_run(
  times=times_train,
  delta=delta_train,
  X=X_train,
  Z=Z_train,
  m_stop=M,
  boost_intercepts_continually=FALSE,
  should_print=TRUE
)
