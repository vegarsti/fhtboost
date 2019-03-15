#' @export

generate_clinical <- function(n.obs=200,tot.genes=10000,n.groups=100,n.clin=NULL,n.gene=NULL,mean.n.gene=15,
                        mu.g=6,sigma.g=0.65,mu.c=1,sigma.c=0.5,rho.c=0,rho.b=0,rho.g=0,phi=0.1,nu=10,tau=20) {
  # n.obs (integer) = number of observations
  # tot.genes (integer) = number of molecular predictors
  # n.groups (integer) = number of pathways (molecular predictors correlated to each other)
  # n.clin (vector) = number of clinical predictors for each group (if NULL, no clinical predictors are generated)
  # n.gene (vector) = number of molecular predictors for each group (if NULL, all the group sizes are generated randomly, if its length is smaller than n.groups the unspecified size are generated randomly as well)
  # mean.n.gene (integer) = average sizes of molecular predictors for group. Relevant only if n.gene is NULL or length(n.gene)<n.groups
  # mu.c (integer) = mean of the log-normal distribution of the clinical variables
  # sigma.c (integer) = standard deviation of the log-normal distribution of the clinical variables
  # mu.g (integer) = mean of the log-normal distribution of the genes
  # sigma.g (integer) = standard deviation of the log-normal distribution of the genes
  # rho.g (vector) = correlation within each block of genes (equal for all of them). Default (for all or only the part of the groups for which it is not specified) is 0.
  # rho.c (vector) = vector containing the correlation between the clinical predictors in each pathway
  # rho.b (vector) = vector containing the correlation between the clinical and the molecular predictors in each pathway
  # phi (integer) = standard deviation of the normal distribution modeling the multiplicative noise to the signal
  # nu (integer) = mean of the normal distribution modeling the additive noise to the signal
  # tau (integer) = standard deviation of the distribution modeling the additive noise to the signal

  require(mvtnorm)
  require(corpcor)

  if(tot.genes<n.groups) stop('The number of genes must be bigger than the number of pathways\n')
  # if the groups sizes are not provided, generate them
  length.n.gene<-length(n.gene)
  if(tot.genes<sum(n.gene)) stop('The number of genes must be bigger than the total number of genes in the pathways\n')
  {if(is.null(n.gene)) length.n.gene<-0
    else if (length.n.gene<n.groups) warning(paste0('The length of n.gene is smaller than ',n.groups,'. The sizes of the remaining groups is generated randomly\n'))}
  if (length.n.gene<n.groups)
  {
    n.gene<-c(n.gene,round(rnorm(n.groups-length.n.gene,mean.n.gene,0.3*mean.n.gene)))
    n.gene[n.gene<1]<-1 # to have no empty group
  }
  # update the number of groups generated
  sumBs<-sum(n.gene)
  # if we generate more genes than those indicated in tot.gene, tell how many genes are actually generated
  {if(sumBs>tot.genes)
  {
    tot.genes<-sumBS
    warnings(paste0('Total number of simulated genes equal to ',sumBS,'\n'))
  }
    else n.gene[n.groups+1]<-tot.genes-sumBs} # the last block contains all the genes not belonging to the first n.groups blocks
  # if the correlation among genes is not specified (for all or only part of the groups), we set it to 0
  if(length(rho.g)<n.groups) rho.g<-c(rho.g,rep(0,n.groups-length(rho.g)))

  # check for the clinical structure and complete the lists
  ifelse(is.null(n.clin),length.n.clin<-0,length.n.clin<-length(n.clin))
  # check if the number of clinical groups is reasonable
  if(length.n.clin>n.groups) {
    n.clin<-n.clin[1:n.groups]
    warnings(paste0('Number of clinical groups too large, only the first ',length.n.gene,' are used.\n'))
  }
  # n.clin[(length.n.clin+1):n.groups]<-0
  # this caused an error for me. it also set the last ones to zero

  # if correlation among clinical predictors is not specified (for all or only part of the groups) set it to 0
  if(length(rho.c)<n.groups) rho.c<-c(rho.c,rep(0,n.groups-length(rho.c)))
  # if correlation between clinical and molecular predictors is not specified (for all or only part of the groups) set it to 0
  if(length(rho.b)<n.groups) rho.b<-c(rho.b,rep(0,n.groups-length(rho.b)))

  # generate the data
  Clin<-NULL
  Gene<-NULL
  for (i in 1:n.groups)
  {
    # generate the covariance matrix
    Sigma.h<-rbind(cbind(rep(sigma.c,n.clin[i])%*%t(rep(sigma.c,n.clin[i]))*rho.c[i], # clinical part
                         rep(sigma.c,n.clin[i])%*%t(rep(sigma.g,n.gene[i]))*rho.b[i]), # clinical and molecular part, up right part of the covariance matrix
                   cbind(t(rep(sigma.c,n.clin[i])%*%t(rep(sigma.g,n.gene[i])))*rho.b[i], # clinical and molecular part, bottom left part of the covariance matrix
                         rep(sigma.g,n.gene[i])%*%t(rep(sigma.g,n.gene[i]))*rho.g[i])) # molecular part
    diag(Sigma.h)<-c(rep(sigma.c^2,n.clin[i]),rep(sigma.g^2,n.gene[i]))
    if(!is.positive.definite(Sigma.h)) Sigma.h<-make.positive.definite(Sigma.h)

    # generate the data
    tmp<-rmvnorm(n.obs,c(rep(mu.c,n.clin[i]),rep(mu.g,n.gene[i])),Sigma.h)
    {if(n.clin[i]==0) Gene<-cbind(Gene,exp(tmp))
      else
      {
        Clin<-cbind(Clin,tmp[,1:n.clin[i]])
        Gene<-cbind(Gene,exp(tmp[,-c(1:n.clin[i])]))
      }}
  }
  # the genes exceeding those clustered in the groups are generated uncorrelated to any other predictors
  if(tot.genes>sumBs) Gene<-cbind(Gene,sapply((sumBs+1):tot.genes,function(i,mu.g,sigma.g,n.obs) exp(rnorm(n.obs,mu.g,sigma.g)),mu.g=mu.g,sigma.g=sigma.g,n.obs=n.obs))

  # generation of the noise
  # additive noise
  E<-matrix(rnorm(tot.genes*n.obs,nu,tau),ncol=tot.genes,nrow=n.obs)
  # multiplicative noise
  M<-matrix(rnorm(tot.genes*n.obs,0,phi),ncol=tot.genes,nrow=n.obs)
  # observations including the noise
  Gene<-Gene*exp(M)+E

  # thresholding and normalization
  Gene[Gene<10]<-10
  Gene[Gene>16000]<-16000
  Gene<-log(Gene)

  if (!is.null(Clin)) colnames(Clin)<-paste0('clin',1:dim(Clin)[2])
  colnames(Gene)<-paste0('gene',1:dim(Gene)[2])

  return(list(clin=Clin, gene=Gene))
}
