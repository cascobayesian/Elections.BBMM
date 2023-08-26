library(MultiRNG)


init.model.mix<-function(data.set,K=5,alpha=NULL,C=1000)
{
  N <- dim(data.set)[1] ## number of samples (towns)
  Q <- dim(data.set)[2] ## number of questions / assays
  S <- dim(data.set)[3] ## number of species, 2 (or 3) for voter data
  alpha <- array(0,dim=c(K,Q,S)) ## parameters for DM
  eta <- rdiric(1,rep(1,K)) ## prior on r
  
  
  #if (is.null(alpha))
  #{
    for(k in 1:K)
    {
      for(q in 1:Q)
      {
        alpha[k,q,] <- rgamma(S,1,0.1)
      }
    }
  #}
  
  check = 1
   
  r <- rep(1:K,1000)[1:N]
  sample <- sample(1:N,N,replace=FALSE)
  r <- r[sample]
  
  ## prior on the multinomial component
  delta <- rep(1,S)
  
  ## a.0, b.0 : priors on the dispersion component
  a.0 <- 10
  b.0 <- 0.1
  
  lambda <- alpha ##multinomial component of DM
  for(k in 1:K)
  {
    for(q in 1:Q)
    {
      lambda[k,q,] <- alpha[k,q,]/sum(alpha[k,q,])
    }
  }
  
  ## gibbs sampling augmentation parameter
  
  Z <- array(0,dim=c(K,Q,S))
  Z <- gibbsZ(r, data.set,alpha)
  
  model <- list(alpha,delta,lambda,eta,N,K,Q,S,r,Z,a.0,b.0,NULL,NULL,NULL)
  names(model) <- c("alpha","delta","lambda","eta","N","K","Q","S","r","Z","a.0","b.0","llk.mat","llk","lpr")
  
  model$lpr <- calc.lpr(model)
  model$llk.mat <- calc.llk.mat(data.set,model$alpha)
  model$llk <- sum(diag(model$llk.mat[1:model$N,model$r]));
 
  return(model)
}


## dirichlet-multinomial distribution
dirichlet.multinomial <- function(x,alpha)
{
  N <- sum(x); A <- sum(alpha);
  if (is.na(N))
  {
    return(0)
  }else
  {
     v <- lfactorial(N) + lgamma(A)-lgamma(N+A) + sum(lgamma(x+alpha)-lgamma(alpha)-lfactorial(x))
   }
  return(v)
}

calc.alpha.lpr <- function(model)
{
  
  A <- 0*model$lambda[,,1]
  for(s in 1:model$S)
  {
    A <- A + model$alpha[,,s]
  }
  #print(dim(A))
  #print(A)
  #print(dgamma(A,model$a.0,model$b.0,log=TRUE))
  ## assumes symmetric dirichlet prior is 0!! (not great coding)
  return(sum(dgamma(A,model$a.0,model$b.0,log=TRUE)))
}

calc.r.lpr <- function(model)
{
  return(dmultinom(table(factor(model$r,1:model$K)),prob=model$eta,log=TRUE))
}

calc.eta.lpr <- function(model)
{
  return(0)#og(ddirichlet(model$eta,rep(1,model$K))))
}


calc.lpr <- function(model)
{
  #	log.z <- dmultinom(table(factor(model$z,0:model$K)),rep(1/model$N,N),log=TRUE)
  
  a.lpr <- calc.alpha.lpr(model)
  
  ## calculates the r prior

  r.lpr <- 0#calc.r.lpr(model)
  
  ## calculates the eta prior
  
  eta.lpr <- calc.eta.lpr(model)
 
  lpr <- a.lpr + r.lpr + eta.lpr;
  return(lpr)
}

calc.llk.i <- function(i,mat,pos)
{
  return(mat[i,pos[i]])
}

calc.llk <- function(model)
{
  lpr <- sum(as.vector(table(factor(model$r,1:model$K)))*log(model$eta))
  llk <- sum(unlist(lapply(1:length(model$r),calc.llk.i,model$llk.mat,model$r)))
  return(llk+lpr)
}
