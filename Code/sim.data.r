sim.data.bb.mix <- function(C,Q,N,L,S=2,alpha=0.01)
{
	## C :number of counts
	## Q :number of questions
	## N :number of towns
	## L :number of cultures
  
  print(c(C,N,L,Q,S))
  
	mix <- rdiric(N,rep(alpha,L));
	data.set <- array(0,dim=c(N,Q,S))

	print(c(L,Q,S))
	p <- array(0,dim=c(L,Q,S))
	a <- p;

	for(q in 1:Q)
	{
	  a[,q,] <- sim.alpha.vecs(S,K = L)
		for(l in 1:L)
		{
			p[l,q,] <- a[l,q,]/sum(a[l,q,])
		}
	}
	for(i in 1:N)
	{
		for(q in 1:Q)
		{
			for(l in 1:L)
			{
			  data.set[i,q,] <- data.set[i,q,] + Dirichlet.multinomial(round(mix[i,l]*C),a[l,q,]);
			}
		}
	}
	sims <- list(data.set,p,a,alpha,mix)
	names(sims) <- c("data.set","p","a","alpha","mix")
	return(sims)
}

sim.data.dm <- function(N,alpha,C)
{
    
    data.set <- Dirichlet.multinomial(rep(C,N),alpha)
    return(data.set)
}

sim.alpha.vecs <- function(S = 15, K = 5,eps = 0.5, alpha = 0.3)
{
  check <- 1
  if ((K>1)&(S>1))
  {
    while(check)
    {
      p <- rdiric(K,rep(alpha,S))
      d.m <- dist(p)
     # print(max(d.m))
      #if (max(d.m)>(1/S))
      #{
        check = 0
      #}
    }
  #  print("done")
    a <- p*0
    
    for(k in 1:K)
    {
      a[k,]<-p[k,]*rgamma(1,S,0.1)
    }
  }else
  {
    p <- rdiric(K,rep(alpha,S))
    a <- p*rgamma(1,10,0.1)+1;
  }
  
  return(a)
}


sim.data.dm.mix <- function(N,S,C=1000,K=5,alpha=0.01)
{
  data.set <- array(0,dim=c(N,S))
  alphas <- sim.alpha.vecs(S,K)
  
  mix <- rdiric(N,rep(alpha,K))
  
  for(i in 1:N)
  {
    for(k in 1:K)
    {
      data.set[i,]  <- data.set[i,] + Dirichlet.multinomial(round(C*mix[i,k]),alphas[k,])
    }
  }
  sims <- list(data.set,mix,alphas)
  
  names(sims) <- c("data.set","mix","alphas")
  return(sims)
}
