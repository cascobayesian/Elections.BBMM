#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;

// [[Rcpp::export]]
double calculateSum(int K, int Q, const arma::cube A)
{
  int S=A.n_slices,s=0;
  double total=0;
  
  for (s=0;s<S;s++)
  {
    total += A(K,Q,s);
  }
  return(total);
}

// [[Rcpp::export]]
NumericVector rdirichletCpp(NumericVector alpha_m) 
{
  int distribution_size = alpha_m.length();
  NumericVector distribution(distribution_size);
  distribution.fill(0.0);
  
  double sum_term = 0;
  // loop through the distribution and draw Gamma variables
  for (int j = 0; j < distribution_size; ++j) {
    double cur = R::rgamma(alpha_m[j],1.0);
    distribution(j) = cur;
    sum_term += cur;
  }
  // now normalize
  for (int j = 0; j < distribution_size; ++j) {
    distribution(j) = distribution(j)/sum_term;
  }

  return(distribution);
}


// [[Rcpp::export]]
NumericVector gibbsEta(int K,const NumericVector& r,double etaPrior)
{
  int N=r.length(),i=0;
  NumericVector counts(K);
  counts.fill(etaPrior);
  for(i=0;i<N;i++)
  {
    counts[r[i]-1] = counts[r[i]-1]+etaPrior;
  }
  
  return(rdirichletCpp(counts));
}


// [[Rcpp::export]]
arma::cube gibbsLambda(const arma::cube& alpha,const arma::cube& Z)
{
  arma::cube lambda(Z.n_rows,Z.n_cols,Z.n_slices);
  
  //Z.tube(i,q)
  NumericVector test(Z.n_slices);
  NumericVector test2(Z.n_slices);
  test2.fill(0);
  int i, q,s;
  double omega;
  for( i = 0; i < Z.n_rows; i++)
  {
    for(q = 0; q < Z.n_cols; q++)
    {
      for(s = 0; s< Z.n_slices;s++)
      {
        test2[s] = Z(i,q,s)+1.0;
      }
      
      //Rcpp::NumericVector x = 
      test = rdirichletCpp(test2);//gibbsLambdaVector(d, test2);
      
      omega = calculateSum(i, q, alpha);
      for(s =0; s < Z.n_slices;s++)
      {
        lambda(i,q,s) = test[s]*omega;
      }
    }
  }
  return(lambda);
};

// omega
// [[Rcpp::export]]
NumericVector calculateNs(const NumericVector& r, int k, int q,const arma::cube& dataSet)
{
  int N = r.length(), i =0, total=0, s=0, S=dataSet.n_slices;
  
  for(i =0 ; i<N; i++)
  {
    if((r[i]-1)==k)
    {
      total += 1;
    }
    
  }
  //Rcpp::Rcout<<"Total:"<<total<<arma::endl;
  
  NumericVector Ns(total);
  Ns.fill(0);
  
  total = 0;
  double dataTotal;
  
  for(i =0 ; i<N; i++)
  {
    if((r[i]-1)==k)
    {
      dataTotal=0;
      //Rcpp::Rcout<<" "<<dataTotal<<arma::endl;
      for(s =0; s<S; s++)
      {
        //Rcpp::Rcout<<dataSet(i,q,s)<<" ";
        dataTotal += dataSet(i,q,s);
        //Rcpp::Rcout<<" "<<dataTotal<<arma::endl;
      }
      //Rcpp::Rcout<<" "<<dataTotal<<arma::endl;
      Ns[total] = dataTotal;
      total = total+1;
    }
  }
  return(Ns);  
}



double targetOmega(double omega,const double& Z,const NumericVector& N,int S)
{
  int i = 0, n = N.length();
  double target = 0;
  
  for(i=0; i < n; i++)
  {
    target += lgamma(omega)-lgamma(omega+N[i]);
  }
  //Rcpp::Rcout<<ps<<arma::endl;
  target +=  Z*log(omega);
  target +=  R::dgamma(omega,S,100,true);
  return(target);
}


//gibbs.grid.omega <- function(Z,)
// [[Rcpp::export]]
double mhOmegaIter(const double& Z,double omega,const NumericVector& N,int S,double rate=2,int T=20)
{
  int t=0;
  
  double omegaPrime=0,logM=0,logH=0,logMH=0;
  
  for(t =0; t<T; t++)
  {
    omegaPrime = omega*(Rcpp::runif(1,rate/(1+rate),(1+rate)/rate)[0]);
    logM = targetOmega(omegaPrime,Z,N,S) - targetOmega(omega,Z,N,S);
    logH = 0;
   
    //Rcpp::Rcout<<omegaPrime<<arma::endl<<logM<<" "<<logH<<arma::endl;
    //print(c(logM,logH));
    logMH = logM+logH;
   
     if (log(Rcpp::runif(1,0,1)[0])<logMH)
     {
        omega=omegaPrime;
      }
   }
// Rcpp::Rcout<<accept<<arma::endl;
   return(omega);
}

// [[Rcpp::export]]
arma::cube& mhOmega(NumericVector r,const arma::cube& Z,const arma::cube& dataSet,arma::cube& alpha)
{
  int k=0,q=0,s =0,K = Z.n_rows, Q = Z.n_cols, S = Z.n_slices;

  double omega = 0, A =0,z=0,A2=0;
  NumericVector N;
  
  for(k=0;k<K;k++)
  {
    for(q=0;q<Q;q++)
    {
      
      
      N =calculateNs(r, k, q, dataSet);
      A2 =0;
      for(s =0;s<S;s++)
      {
        A2 += alpha(k,q,s);
      }
      
      A = calculateSum(k,q,alpha);
      z = calculateSum(k,q,Z);
    //  Rcpp::Rcout<<A2<<" "<<A<<arma::endl;
      omega = mhOmegaIter(z,A,N,S);  
      //omega = A;
    //  Rcpp::Rcout<<"Pair:"<<k<<","<<q<<arma::endl;
    //  Rcpp::Rcout<<"Values:"<<omega<<" "<<A<<arma::endl;
     // Rcpp::Rcout<<k<<" "<<" "<<q<<" "<<calculateSum(k,q,Z)<<" "<<A<<" "<<omega<<arma::endl;
      
      for(s=0; s<S;s++)
      {
        alpha(k,q,s) = alpha(k,q,s)*omega/A;
      }
      A = calculateSum(k,q,alpha);
  //    Rcpp::Rcout<<"Sum:"<<A<<arma::endl;
      
    }
  }
  return(alpha);
}

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
// alpha : mh updates
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

// [[Rcpp::export]]
double sumAlpha(const arma::cube& alpha,int k,int q,int s)
{
  double notAlpha =0;
  int S = alpha.n_slices;
  
  for(int r=0; r < S; r++)
  {
    if (r!=s)
    {
      
      notAlpha += alpha(k,q,r);
    }
  }
  return(notAlpha);
}


// [[Rcpp::export]]
double gammaCheck(double kappa, double theta)
{
  return(R::dgamma(10,kappa,theta,true));
}

// [[Rcpp::export]]
double mhAlphaPosition(const double& Z, const double& alpha, const double& notAlpha,const NumericVector& Ns,int T=20)
{
  int n = Ns.length(),t=0,i; 
  
  double omega = alpha + notAlpha, omegaP =0;
  double alphaP =0,a=0;
  double logM =0, logH =0, logMH =0;
  double u=0;
  
  
  a = alpha;
  for(t = 0; t < T; t++)
  {
    omega = a + notAlpha;
    alphaP = rgamma(1,1,20)[0];
    omegaP = alphaP + notAlpha;
    
    logM  = R::dgamma(alphaP,1+Z,100,true);
    logM += (-R::dgamma(a,1+Z,100,true));
    
    logH = 0;
    
    if (n >0)
    {
      for(i =0; i < n; i++)
      {
        logH += (lgamma(omegaP) - lgamma(omegaP+Ns[i]));
        logH += ((-lgamma(omega)) + lgamma(omega+Ns[i]));
      }
    }
    
    logH += R::dgamma(a,1,20,true)-R::dgamma(alphaP,1,20,true);
    
    logMH = logM+logH;
    u = runif(1,0,1)[0];
    
    if (log(u)<logMH)
    {
      a = alphaP;
    }
  }
  return(a);
}


// k,q,s, :: Z, alpha, not.alpha, Ns :: data Set
// [[Rcpp::export]]
arma::cube& mhAlpha(const NumericVector& r,const arma::cube& Z,arma::cube& alpha,const arma::cube& dataSet,int T=20)
{
  int K = Z.n_rows, Q = Z.n_cols, S = Z.n_slices;
  int k=0,q=0,s=0;
  
  double notAlpha;
  
 //arma::cube newAlpha = alpha;

  NumericVector Ns;
  // calculatue the Ns
  
  for(k = 0; k < K; k++)
  {
    for(q = 0; q < Q; q++)
    {
      
      //k = 0; q = 0;
      //Rcpp::Rcout<<"Pos:"<<k<<" "<<q<<":"<<arma::endl;
      Ns = calculateNs(r, k,q,dataSet);
      //Rcpp::Rcout<<Ns[0]<<arma::endl;
      for(s = 0; s < S; s++)
      {
        //s = 1;
        
        notAlpha = sumAlpha(alpha,k,q,s);
        alpha(k,q,s) = mhAlphaPosition(Z(k,q,s), alpha(k,q,s),notAlpha,Ns);
      }
    }
  }
  return(alpha);
}

///////////////////////////////////////////
///////////////////////////////////////////
// Z updates //////////////////////////////
///////////////////////////////////////////
///////////////////////////////////////////


// [[Rcpp::export]]
double rbinomSub(double datum, double alpha)
{
  double totalZ=0;
  
  if (datum==0)
  {
    totalZ = 0;	
  }else
  {
    
    double ps;
    Rcpp::NumericVector z(1);
    
    for(int k=1; k<=datum; k++)
    {
      ps = std::min((alpha/(alpha+k-1)),1.0);
      z  =  rbinom(1,1,ps);
      totalZ += z[0];
    }
      //Rcpp::Rcout<<ps<<arma::endl;
    
  }
  return(totalZ);
}


// [[Rcpp::export]]
arma::cube gibbsZ(NumericVector r, const arma::cube& dataSet,const arma::cube& alpha)
{
  int i,q,s;
  arma::cube Z(alpha.n_rows,alpha.n_cols,alpha.n_slices);
  Z.zeros();
  
  //Rcpp::Rcout<<Z.n_rows<<" "<<Z.n_cols<<" "<<Z.n_slices<<arma::endl;
  for( i = 0; i < dataSet.n_rows; i++)
  {
    for(q = 0; q < alpha.n_cols; q++)
    {
      for( s = 0; s < alpha.n_slices; s++)
      {
       // Rcpp::Rcout<<r(i)-1<<" "<<q<<" "<<s<<std::endl;
        Z(r(i)-1,q,s) += rbinomSub(dataSet(i,q,s), alpha(r(i)-1,q,s));
        
      }
    }
  }
  return(Z);
}

//////////////////////////////////////////////
//////////////////////////////////////////////
//// GIBBS COMPONENT ///////////////////////////////
//////////////////////////////////////////////
//////////////////////////////////////////////

// [[Rcpp::export]]
double dirichletMultinomial(NumericVector dataSet, NumericVector alpha)
{
  
  double N = 0, A = 0,total=0; 
  
  int check = 0;
  for(int k=0; k<dataSet.length();k++)
  {
    if (NumericVector::is_na(dataSet[k]))
    {
      check=1;
    }

  }
  if (check ==0)
  {
    for(int k =0; k < alpha.length(); k++)
    {
      total += lgamma(dataSet[k]+alpha[k]);
      total -= Rcpp::internal::lfactorial(dataSet[k]);
      total -= lgamma(alpha[k]);
      A += alpha[k];
      N += dataSet[k];
    }
    total += lgamma(A);
    total += Rcpp::internal::lfactorial(N);
    total -= lgamma(N+A);
  }else
  {
    total = 0;
    
  }
  return(total);
}

NumericMatrix calcLLKMat(const arma::cube& alpha,const arma::cube& dataSet, const NumericVector eta)
{
  int i = 0, q =0, k = 0, s =0, N = dataSet.n_rows, K=alpha.n_rows, Q = alpha.n_cols, S = alpha.n_slices;
  
  NumericMatrix llkMat(N,K);
  NumericVector d(S), a(S);
  llkMat.fill(0);
  for(i = 0; i < N; i++)
  {
    for(k =0; k <K; k++)
    {
      for(s =0; s <S; s++)
      {
        d(s) = dataSet(i,q,s);
      }
      for(q=0; q<Q; q++)
      {
        for(s =0; s <S; s++)
        {
          a(s) = alpha(k,q,s);
        }
        //a=alpha(k,q,_);
        llkMat(i,k)+= dirichletMultinomial(d,a);
      }
    }
  }
  return(llkMat);
}

// [[Rcpp::export]]
int sampleInt(int K) 
{
  Rcpp::IntegerVector pool = Rcpp::seq(0, K-1);
  std::random_shuffle(pool.begin(), pool.end());
  return pool[0];
} 

NumericVector createProbVector(const NumericVector& logProb,const NumericVector& probs)
{
  int k, K = logProb.length();
  NumericVector returnProb(K);
  returnProb.fill(0);
  double maxValue = R_NegInf, total=0;
  
  for(k =0; k <K; k++)
  {
      if (maxValue < logProb[k])
      {
        maxValue = logProb[k];
      }
  }
  for(k =0; k <K; k++)
  {
    returnProb[k] = exp(logProb[k]-maxValue)*probs[k];
    if (returnProb[k]<1e-300)
    {
      returnProb[k] = 1e-300;
    }
    total += returnProb[k];
  }
  for(k =0; k < K; k++ )
  {
    returnProb[k] = returnProb[k]/total;
  }

  return(returnProb);
}

// [[Rcpp::export]]
NumericVector gibbsComponent(const arma::cube& alpha,const arma::cube& dataSet, const NumericVector eta)
{
  int i = 0, q =0, k = 0, s =0, N = dataSet.n_rows, K=alpha.n_rows, Q = alpha.n_cols, S = alpha.n_slices;
  
  NumericMatrix llkMat(N,K);
  NumericVector d(S), a(S), r(N), vec(K), probs(K),llkVec(K);
  llkMat.fill(0); llkVec.fill(0);
  r.fill(0);
  
  vec = Rcpp::seq(1, K);
  
  for(i = 0; i < N; i++)
  {
    for(k =0; k <K; k++)
    {
      for(q=0; q<Q; q++)
      {
        for(s =0; s <S; s++)
        {
          d(s) = dataSet(i,q,s);
          a(s) = alpha(k,q,s);
        }
       // Rcout <<i <<" "<<k<<" "<<q<<" "<<s<<"\n";
       // Rcout << "The value of a : " << a << "\n";
       // Rcout << "The value of d : " << d << "\n";
       // Rcout << dirichletMultinomial(d,a) <<"\n";

        llkMat(i,k)+= dirichletMultinomial(d,a);
       // Rcout <<llkMat(i,k)<<"\n";

      }
      llkVec(k) = llkMat(i,k);
    }
    probs = createProbVector(llkVec,eta);
  
    r(i) = RcppArmadillo::sample(vec,1,true,probs)[0];
  }
  
  return(r);
}

