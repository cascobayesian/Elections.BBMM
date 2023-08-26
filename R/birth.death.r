
calc.mix.lprob.l <- function(model,l)
{
 # print(c(model$K,l))
  lpr <- (as.vector(table(factor(model$r,1:model$K)))[l]*log(model$eta[l]))
  #print(lpr)
  return(lpr)
}

calc.mix.lprob <- function(model)
{
  lpr <- sum(as.vector(table(factor(model$r,1:model$K)))*log(model$eta))
  return(lpr)
}

kill.k<-function(m,dataSet,k)
{
  model <- m;
  
  if (m$K>1)
  {
    set <- 1:m$K;
    set <- set[-c(k)]
    
    model$K <- m$K-1;
    model$Z <- array(0,dim=c(model$K,model$Q,model$S))
    model$alpha <- array(0,dim=c(model$K,model$Q,model$S))
    model$alpha[1:model$K,,] <- m$alpha[set,,]
    model$Z[1:model$K,,] <- m$Z[set,,]
    
    model$r <- m$r;
    model$r[m$r==k] = (-1)
    model$r[m$r>k] = model$r[m$r>k]-1
    model$r[m$r==k] = sample(1:model$K,sum(m$r==k),replace=TRUE)
    
    model$eta <- m$eta[set]/(1-m$eta[k])
    model$llk.mat <- calc.llk.mat(dataSet,model$alpha)
    model$llk <- calc.llk(model)
    
  }else
  {
    print('Issue: K == 1 and trying to die.')
  }
  return(model)
  
}

d.k = function(model,dataSet)
{
  death.models <- rep(list(),model$K)
  
  for(k in 1:model$K)
  {
    death.models[[k]] <- kill.k(model,dataSet,k)
  }
  return(death.models)
  
}

birth <- function(m,dataSet)
{
  model<-m;
  
  new.alpha <- array(rgamma(model$Q*model$S,model$a.0,model$b.0),dim=c(model$Q,model$S))
  new.Z     <- array(sample(1:100,model$Q*model$S,replace=TRUE),dim=c(model$Q,model$S))
  
  model$K <- m$K+1;
  model$Z <- array(0,dim=c(model$K,model$Q,model$S))
  model$alpha <- array(0,dim=c(model$K,model$Q,model$S))
  model$alpha[1:m$K,,] <- m$alpha[1:m$K,,]
  model$Z[1:m$K,,] <- m$Z[1:m$K,,]
  
  model$alpha[m$K+1,,] <- new.alpha;
  model$Z[m$K+1,,] <- new.Z;
  
  model$r <- m$r;
  new.r <- rbinom(model$N,1,prob=c(1/(model$K)))
  model$r[new.r==1] <- model$K
  
  model$eta <- c(m$eta,runif(1))
  model$eta <- model$eta/sum(model$eta)
  
  model$llk.mat <- calc.llk.mat(dataSet,model$alpha)
  model$llk <- calc.llk(model)
  return(model)
}

calc.death.prob = function(death.models,old.model,beta=1)
{
  log.bd  =  log(beta) - log(old.model$K) 
  log.kpr =  dpois(old.model$K-1,10,log=TRUE) - dpois(old.model$K,10,log=TRUE) 
  
  log.llk <- rep(0,old.model$K)
  log.mix <- rep(0,old.model$K)
  
  for(k in 1:old.model$K)
  {
     log.llk[k] <- death.models[[k]]$llk - old.model$llk;
     log.mix[k] <- calc.mix.lprob(death.models[[k]])  - calc.mix.lprob(old.model)
                      
    #                
    
  }
  
  log.death.rates <-log.llk+log.mix + log.bd + log.kpr#,old.model$K)
  #log.llk+log.mix;
  #print(log.llk)
  #print(log.mix)
  #print(log.death.rates)
  p = exp(log.death.rates);
  p[p == Inf] = 1e300;
  p[p == 0]   = 1e-300;
  return(p)
  
}


birth.death <- function(model,dataSet,beta=1)
{
  if (model$K>1)
  {
    death.models <- d.k(model,dataSet)
    death.prob <- calc.death.prob(death.models,model,beta)
    
    death.rate <- sum(death.prob)
    birth.rate <- beta;
    
    #print(death.prob)
    #print(c(birth.rate,death.rate))
    t <- rexp(1,death.rate+birth.rate)
    u <- rbinom(1,1,birth.rate/(birth.rate+death.rate))
    
    #print(birth.rate)
    if (u==1)
    {
      model <- birth(model,dataSet)
      return(list(model,t))
    }else
    {

      k <- sample(1:model$K,1,replace=TRUE,death.prob)
      return(list(death.models[[k]],t))
    }
    
  }else
  {
    t <- rexp(1,1)
    model <- birth(model,dataSet)
    return(list(model,t))
  }
}

