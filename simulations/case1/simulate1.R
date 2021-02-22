rm(list=ls())

library(MASS)
library(far)

# remotes::install_github("isglobal-brge/mgcca")
library(curatedTCGAData)
library(TCGAutils)

# Libraries needed by mgcca package
library(MultiAssayExperiment)
library(RSpectra)
library(impute)

# Benchmark
library(microbenchmark)
library(mgcca)


fsim <- function(n, p, q, prop.miss, ncomp, ncomp.true, sigma.noise){

  # n <- 100
  # p <- 5
  # q <- 5
  # prop.miss <- 0.1
  # ncomp <- 2
  # sigma.noise <- 0.001
  # ncomp.true <- 2

  # n = number of individuals
  # p = number of variables in X matrix
  # q = number of variables in Y matrix
  # prop.miss= proportion of missing values
  # ncomp = number of correlation functions
  # sigma.noise = sigma of noise added to Xi.
  # ncomp.true = number of true canonical variables.
  # ... = arguments passed to mgcc function (currently ignored)

  if (p>n | q>n) stop("number of individuals must be larger than the number of variables")

  # simulate data (as in "Generalized canonical correlation analysis of matrices with missing rows: A simulation study" paper)
  # note: we add noise to X1 instead of Y, if noise is added to Ytrue as specified in paper var(Xi) is not singular!!!
  Ysim <- matrix(rnorm(n*ncomp.true), n, ncomp.true) # simulate from a standard normal variables
  Ytrue <- Ysim
  b1 <- matrix(runif(ncomp.true*p,-1,1), ncomp.true, p)
  X1 <- Ytrue%*%b1 + matrix(rnorm(n*p, 0, sigma.noise), n, p) # add noise to X1 instead of Y!!!
  b2 <- matrix(runif(ncomp.true*q,-1,1), ncomp.true, q)
  X2 <- Ytrue%*%b2 + matrix(rnorm(n*q, 0, sigma.noise), n, q) # add noise to X2 instead of Y!!!

  rownames(X1) <- rownames(X2) <- 1:n
  X<-Xori <- list(X1,X2)

  ####### introduce individuals missings at random
  n<- nrow(X[[1]])
  Xmiss <- Xori
  k <- length(Xmiss)
  which.miss.prev <- integer()
  for (i in 1:k){
    Xi <- X[[i]]
    n <- nrow(Xi)
    p <- ncol(Xi)
    which.miss.i <- sample(1:n,round(prop.miss*n))
    which.miss.i <- which.miss.i[!which.miss.i%in%which.miss.prev]
    if (length(which.miss.i)>0) Xi[which.miss.i,] <- NA
    Xmiss[[i]] <- Xi
    which.miss.prev <- which.miss.i
  }

  ## average imputation
  Ximp <- lapply(Xmiss, function(Xi){
    apply(Xi,2, function(xx){
      m <- mean(xx,na.rm=TRUE)
      ifelse(is.na(xx), m , xx)
    })
  })

  #### fit mgcc to different methods

  ## gold (all cases)
  resultAll <- mgcca(Xori, nfac = ncomp, scores=TRUE)

  ## average imputation
  resultImp <- mgcca(Ximp, nfac = ncomp, scores=TRUE)

  ## complete cases
  common <- Reduce(intersect,lapply(lapply(Xmiss, na.omit),rownames))
  Xcommon <- lapply(Xmiss, function(xi) xi[common,,drop=FALSE])
  resultCommon <- mgcca(Xcommon, nfac = ncomp, scores=TRUE)
  #cat("n common=", nrow(Xcommon[[1]]),"\n") # number of individuals with complete data

  ## mcca
  Xr <- lapply(Xmiss, na.omit)
  resultMGCCA <- mgcca(Xr, nfac = ncomp, scores=TRUE)

  c(
    mean(rowSums((resultMGCCA$Y - resultAll$Y)^2)),
    mean(rowSums((resultImp$Y - resultAll$Y)^2)),
    mean(rowSums((resultCommon$Y - resultAll$Y[rownames(resultCommon$Y),])^2))
  )

}

## fixed parameters

n <- 500
ncomp <- 2
ncomp.true <- 2
nsim <- 100

## variable parameters
esc <- expand.grid(vars=c(50, 100), prop.miss=c(0.1,0.2,0.3), sigma.noise=c(0.125, 0.250))

res <- list()

set.seed(123456)

for (i in 1:nrow(esc)){
  cat("----Escenari",i,"\n-sigma.noise=",esc$sigma.noise[i],"\n-p,q=",esc$vars[i],"\n-prop.miss=",esc$prop.miss[i],"\n\n")
  iter <- 0
  res[[i]] <- replicate(nsim,{
    iter <<- iter+1
    if (iter%%10==0) cat("iteration",iter,"\n")
    res <- try(fsim(n, p=esc$vars[i], q=esc$vars[i], prop.miss=esc$prop.miss[i], ncomp, ncomp.true, sigma.noise=esc$sigma.noise[i]), silent=TRUE)
    if (inherits(res, "try-error"))
      return(c(NA,NA,NA))
    else
      return(res)
  })
  cat("\n\n\n")
}


save(res, file="sim1res.rda")


