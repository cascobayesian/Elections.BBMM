mh <- function(model,dataSet,T=1,diff=10,verbose=FALSE)
{
  for(t in 1:T)
  {
    #pri
    if (verbose==TRUE)
    {
      print(t)
    }
    if((t%%(diff)==0))
    {
      model$eta <- gibbsEta(model$K,model$r,1);
      model$r<- gibbsComponent(model$alpha,dataSet,model$eta);
    }
   # print(table(r.true,model$r))
    
    model$alpha <- mhAlpha(model$r,model$Z,model$alpha,dataSet);
    model$alpha <- mhOmega(model$r,model$Z,dataSet,model$alpha);
    model$alpha <- gibbsLambda(model$alpha,model$Z);
    model$Z <- gibbsZ(model$r,dataSet,model$alpha);
  
   # model$llk.mat <- calc.llk.mat(dataSet,model$alpha)
   # model$llk <- calc.llk(model)
  #  model$llk.mat <- calc.llk.mat(dataSet,model$alpha)
  #  model$llk <- calc.llk(model)
    
  #  print(model$llk)
    
  }
  model$llk.mat <- calc.llk.mat(dataSet,model$alpha)
  model$llk <- calc.llk(model)
  
 # print(model$llk)
  #model<-remove.elements(model)
  return(model)
}





target.omega <- function(o,Z,N,S)
{
  n <- length(N);
  
  target = 0;
  target <- target + n*lgamma(o)-sum(lgamma(o+N));
  
  target <- target + Z*log(o);
  target <- target + dgamma(o,shape=S,scale=100,log=TRUE)
  return(target)
  
}

gibbs.grid.omega.pos <- function(Z,omega,N,S,eps =0.1)
{
  log.target <- unlist(lapply(seq(0,200,by=eps),target.omega,Z,N,S))
  log.target <- log.target-max(log.target)
  prob.target <- exp(log.target)
  prob.target[prob.target<1e-300]=1e-300
  #sample(seq(0,200,by=eps),1,prob=prob.target))
}

#gibbs.grid.omega <- function(Z,)
mh.omega.iter <- function(Z,omega,N,S,rate=2,T=100)
{
  accept= 0;
  n <- length(N)
  for(t in 1:T)
  {
    omega.prime <- omega*runif(1,rate/(1+rate),(1+rate)/rate);
    
    log.m.1 <- target.omega(omega.prime,Z,N,S) - target.omega(omega,Z,N,S);
    log.m.2 <- 0;
    log.h <- 0;#dgamma(omega.prime,1,log=TRUE)-dgamma(omega,1,log=TRUE);
  
    #print(omega.prime)
   # print(c(log.m.1,log.h))
    log.mh <- log.m.1+log.m.2+log.h;
    #print(log.mh)
    if (log(runif(1))<log.mh)
    {
      accept=accept+1;
      omega=omega.prime;
    }
  }
  #print(accept)
  return(omega)
}

mh.omega <- function(r,Z,alpha,dataSet)
{
  K <- dim(alpha)[1]; Q <- dim(Z)[2]; S<-dim(Z)[3];
  
  for(k in 1:K)
  {
    for(q in 1:Q)
    {
      #print(c(k,q,dim(alpha)))
      omega <- sum(alpha[k,q,]);
      if (sum(r==k)>1)
      {
        N <- rowSums(dataSet[r==k,q,])
      }else
      {
        if(sum(r==k)==1)
        {
          N <- sum(dataSet[r==k,q,])
        }else
        {
          N <- 0;
        }
      }
     # print(c(k,q,s))
     # print(N)
     # print(omega)
      
      omega <- mh.omega.iter(sum(Z[k,q,]),omega,N,S,rate=8)  
      
    #  print('hi')
     # print(alpha[1,1,])
      alpha[k,q,] <- alpha[k,q,]/sum(alpha[k,q,])*omega;
    #  print('hi')
    }
  }
  return(alpha)
}

gibbs.Z<-function(r, alpha,dataSet)
{
  dims <- dim(dataSet)
  Z = array(0,dim=dim(alpha));
  
  numData = dim(dataSet)[1];
  
  for(i in 1:dims[1])
  {
    for(q in 1:dims[2])
    {
      for( s in 1:dims[3])
      {
        #print(c(i,q,s,r[i]))
        Z[r[i],q,s] = Z[r[i],q,s] + rbinom.sub(dataSet[i,q,s], alpha[r[i],q,s]);
      }
    }
  }
  return(Z);
}


rbinom.sub <- function(datum, alpha)
{
  totalZ=0;
  if (datum==0)
  {
    totalZ = 0;	
  }else
  {
    ps <- unlist(lapply(alpha/(alpha+1:datum-1),min,1.0))
    z  =  rbinom(length(ps),1,prob=ps);
    totalZ = sum(z);
  }
  return(totalZ);
};

mh.alpha.iter <- function(Z,alpha,notAlpha,Ns,T=100)
{
  for(t in 1:T)
  {
    alpha <- mh.alpha.pos(Z,alpha,notAlpha,Ns)
  }
  return(alpha)
}

mh.alpha.pos <- function(Z,alpha,notAlpha,Ns)
{
  n <- length(Ns)  
  
  omega <- alpha + notAlpha;
  
  alpha.prime <- rgamma(1,shape=1,rate=0.05);
  
  omega.prime <- alpha.prime + notAlpha;
  
  
  log.m <- dgamma(alpha.prime,shape=1+Z,scale=100,log=TRUE);
  log.m <- log.m - dgamma(alpha,shape=1+Z,scale=100,log=TRUE);
  log.m <- log.m + n*lgamma(omega.prime) - sum(lgamma(omega.prime+Ns)) 
  log.m <- log.m - n*lgamma(omega) + sum(lgamma(omega+Ns));
  
  log.h <- dgamma(alpha,shape=1,rate=0.05,log=TRUE);
  log.h <- dgamma(alpha.prime,shape=1,rate=0.05,log=TRUE);
  
  log.mh <- log.m+log.h;
  
  if (is.na(log.mh))
  {
    print('issue?')
    print(log.mh)
    print(Z)
    print(alpha.prime)
    print(Ns)
    print(alpha)
    print(omega.prime)
    print(omega)
  }
  #print(c(alpha.prime,log.mh))
  if (log(runif(1))<log.mh)
  {
    return(alpha.prime)
  }else
  {
    return(alpha)
  }
}

mh.alpha <- function(r,Z,alpha,dataSet)
{
  K <- dim(alpha)[1]; Q <- dim(alpha)[2]; S <- dim(alpha)[3];
  
  #print(c(K,Q,S))
  new.alpha = alpha*0;
  
  #print(alpha[1,1,])
  for(k in 1:K)
  {
    for(q in 1:Q)
    {
      Ns <- 0;
      if (sum(r==k)>1)
      {
        Ns <- rowSums(dataSet[r==k,q,]);
      }else if (sum(r==k)==1)
      {
          Ns <- sum(dataSet[r==k,q,])
      }
    #  print(alpha[k,q,])
      for(s in 1:S)
      {
        notAlpha <- sum(alpha[k,q,-c(s)])
    #    print(c(k,q,s,notAlpha))
    #    print(c(Z[k,q,s],alpha[k,q,s]))
        new.alpha[k,q,s] = mh.alpha.iter(Z[k,q,s],alpha[k,q,s],notAlpha,Ns);
      }
    }
  }
  return(new.alpha)
}

gibbs.lambda <- function(Z,alpha)
{
  dims <- dim(alpha)
  K <- dims[1]; Q <- dims[2]; 
  
  for(k in 1:K)
  {
    #if (sum(r==k)>0)
    {
      for(q in 1:Q)
      {
        alpha[k,q,] <- sum(alpha[k,q,])*rdiric(1,Z[k,q,]+1);
      }
    }#else
    {
     # i <- sample(1:length(r),1)
    #  Z.i <- dataSet[i,q,]
      for(q in 1:Q)
      {
     #   alpha[k,q,] <- sum(alpha[k,q,])*rdiric(1,Z.i+1);
      }
    }
  }
  return(alpha)
}

gibbs.eta <- function(K,r,eta.prior)
{
  return(rdiric(1,eta.prior + table(factor(r,1:K))));
}


mh.r <- function(L,dataSet,r=NULL,T=1000,thin=10)
{
  model <- init.model.mix(dataSet,L)
  if (is.null(r))
  {
    r <- model$r
  }
  alpha <- model$alpha;
  Z <- gibbs.Z(r,alpha,dataSet);
  
  for(t in 2:T)
  {
    if((t%%20)==0)
    {cat("\n")}
    else{cat(t,",")}
    alpha <- mh.alpha(r,Z,alpha,dataSet);
    alpha <- mh.omega(r,Z,dataSet,alpha);
    alpha <- gibbs.lambda(Z,alpha);
    Z <- gibbs.Z(r,alpha,dataSet);
  }
  return(list(Z,alpha))
}


gibbs.component <- function(r,K,alpha,dataSet,eta)
{
  llk.mat <- calc.llk.mat(dataSet,alpha)
  
  #print(llk.mat[1,1])
  N <- dim(dataSet)[1]

  print(N)
  
  new.r <- r*0;
  for(i in 1:N)
  {
    probs <- exp(llk.mat[i,] -max(llk.mat[i,]))
    probs[probs<(1e-300)]<-1e-300;
   # print(probs)
  #  print(eta)
    new.r[i]<-sample(1:K,1,prob=(probs*eta))
  }
  return(new.r)
  
}

mh.empty.alpha <- function(k,r,alpha,dataSet)
{
  N <- dim(dataSet)[1]
  this.one <- sample(1:N,1)
  
  Q <- dim(alpha)[2]
  S <- dim(alpha)[3]
  
  logM <- 0; logH <- 0;
  
  new.alpha <- array(0,dim=c(Q,S))
  
  #print('thisone:')
  #print(this.one)
  factor <- 1
  for(q in 1:Q)
  {
    A.kq <- sum(alpha[k,q,])
    lambda.kq.old <- alpha[k,q,]/A.kq
    D <- sum(dataSet[this.one,q,])
    lambda.kq.new <- rdiric(1,dataSet[this.one,q,]/D*factor+1)
    new.alpha[q,] <- lambda.kq.new*A.kq;
    logM <- logM + sum(dgamma(new.alpha[q,],1,100,log=TRUE))
    logM <- logM - sum(dgamma(alpha[k,q,],1,100,log=TRUE))
    #print(lambda.kq.new)
    a <- as.vector(dataSet[this.one,q,]/D*factor+1)
    x <- as.matrix(unlist(lambda.kq.old))
    #print(x)
    #print(a)
    logH <- logH + ddirch(t(x),a,log=TRUE)
    x <- as.matrix(unlist(lambda.kq.new))
    logH <- logH - ddirch(t(x),a,log=TRUE)
  }
  logMH <- logM + logH
  
 # print(c(logM,logH))
 # print(logMH)
  if (log(runif(1))< logMH )
  {
    print('accept!')
    alpha[k,,] <- new.alpha
  }
  
  return(alpha)
}


## calculates the llk for a particular question and town, 
## specifying the culture, k
calc.llk.pos.k <- function(k, q, i,data.set,model)
{
  llk <- dirichlet.multinomial(x=data.set[i,q,],model$alpha[k,q,])
  return(llk)
}

calc.llk.pos <- function(q,k,i,data.set,alpha)
{
  #print(dim(alpha))
  #print(dim(data.set))
  #print(c(i,k,q))
  llk <- dirichlet.multinomial(x=data.set[i,q,],alpha[k,q,])
  return(llk)
}

## calculates the llk for a particular question and town, 
## not specifying the culture, using r[i]
calc.llk.question <- function(k, i,data.set,alpha)
{
  Q <- dim(alpha)[2]
  llk <- sum(as.vector(unlist(lapply(1:Q,calc.llk.pos,k,i,data.set,alpha))))
  return(llk)
}


calc.llk.vec.culture <- function(i,data.set, alpha)
{
  K <- dim(alpha)[1]
  return(as.vector(unlist(lapply(1:K,calc.llk.question,i,data.set,alpha))))
}

calc.llk.mat <- function(data.set,alpha)
{
  N <- dim(data.set)[1]
  K <- dim(alpha)[1]
  return(matrix(unlist(lapply(1:N,calc.llk.vec.culture,data.set,alpha)),ncol=K,byrow=TRUE))
}

