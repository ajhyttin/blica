blica_createdata<-function(n=2,u=40,nsources=n,N=Inf,M=NULL,
                           pmvnorm.algorithm=Miwa(steps=4098/2),
                           pairwise=TRUE,verbose=TRUE,kappalim=20) {
  # n - the number of observed variables
  # u - the number of segments
  # nsources - the number of latent sources
  # N - number of samples per segment (default=Inf). Inf means
  #    that we output a distribution that has no small sample problems.
  # M - model to create data from (default=NULL)
  # pmvnorm.algorithm - the algorithm used for multivarite normal CDF calculation  (default=Miwa(steps=4098/2))
  #
  # If N = Inf the output is a list D such that:
  # D$pairwise = TRUE (indicating that only pairwise distributions are given)
  # D$confs  rows have the configurations for a binary variables)
  # D$pairs has the all pairs of n variables for which the distribution is defined
  # 
  # D$X are empty lists delete these.
  # D$Y

  # D$counts has the counts of binary observed variables such that
  # first index is the segment
  # second index is the pair is question, i.e., row of D$pairs
  # third index is the assignments, i.e., row of D$confs
  # 
  # D$type = "infinite" indiciates that we have infinite sample limit data i.e. 
  # distributions instead sample data
  # 
  # D$M defines the model
  #
  # D$M$A the mixing matrix 
  # D$M$mu has means for sources, u (number of segments) x n_sources
  # D$M$sigma log variances for sources, u (number of segments) x n_sources
  # D$M$p the previous parameters in vector form
  #
  # If N is a finite number we get instead:
  # D$pairwise=FALSE as pairwise counts not supported 
  # D$type = "sample", this is sampleR data
  # D$X is a list of data
  # D$X[[i]] is the sample data for segment i (N x n)
  
  
  
  if ( !is.null(M) ) {
    n<-nrow(M$A)
    nsources<-ncol(M$A)
    u<-nrow(M$sigma)
  }
  if ( is.finite(N) ) pairwise=FALSE
  #<-1
  #diag(A)<-c(1,1)
  if ( is.null(M) ) {
    if ( is.na(kappalim) ) {
      kappalim<-est.kappa(n,nsources,times=1000)
    }
    while ( TRUE ) { #used to be without condition number check
      A<-array(runif(n*nsources,min=-3,max=3),c(n,nsources)) #used to be 3
      kappa<-kappa(A)
      if ( kappa <= kappalim ) break;
      cat('rejected kappa=',kappa,'\n')
    }
    cat('kappa:',kappa,'\n')
    #A<-diag(diag(A))
    #A[upper.tri(A)]<-0
    #A[lower.tri(A)]<-0
  #browser()
    mu=array(runif(u*n,min=-.5,max=.5),c(u,nsources)) # used to be 2
    sigma=log(array(runif(u*n,min=0.5,max=3)^2,c(u,nsources))) #following vitoria, used to be without square
  } else {
    n<-nrow(M$A)
    nsources<-ncol(M$A)
    A<-M$A
    mu<-M$mu
    sigma<-M$sigma
    u<-nrow(M$mu)
  }
  
  mus<-array(0,c(0,n))
  
  us<-1:u;
  #R<-array(0,c(length(us),8))
  
  D<-list()
  D$pairwise<-pairwise
  
  if ( pairwise ) {
    D$confs<-array(0,c(2^2,2))
    for ( i in 1:nrow(D$confs)) {
      D$confs[i,]<-dec.to.bin(i-1,2)
    }
    D$pairs<-combinations(n,2)
  } else if ( n <= 10 ) {
    D$confs<-array(0,c(2^n,n))
    for ( i in 1:nrow(D$confs)) {
      D$confs[i,]<-dec.to.bin(i-1,n)
    }
  }
  D$X<-list()
  if ( pairwise ) {
    D$counts<-array(NA,c(u,nrow(D$pairs),4))
  } else {
   # D$counts<-array(NA,c(u,2^n))
  }
  for ( u in us ) {
    #cat(u,'\n')
    if (verbose ) cat(u,'/',max(us),'\n')
    if ( !is.infinite(N) )  {
      Z<-array(0,c(N,nsources))
      for ( i in 1:nsources ) {
        Z[,i]<-rnorm(N,mean=mu[u,i],sd=sqrt(exp(sigma[u,i])))
      }
      mZ<-t(A%*%t(Z))
      prob<-array(0,c(N,0))
      for ( i in 1:n ) {
        prob<-cbind(prob,pnorm(sqrt(pi/8)*mZ[,i]))
        
      }
  
      X<-array(0,c(N,n))
      for (i in 1:N ) {
        for (j in 1:n ) {
          X[i,j]<-sample(c(1,0),1,prob=c(prob[i,j],1-prob[i,j]))
        }
      }
      D$X[[u]]<-X
      D$type="sample"

    } else {
      mean<-as.vector((-1)*sqrt(pi/8)*A%*%mu[u,])
      Su<-diag(exp(sigma[u,]))
      if ( nsources == 1) Su=exp(sigma[u,])
      var<-diag(n)+(pi/8)*A%*%Su%*%t(A)
 
      if ( pairwise ) {
        for ( j in 1:nrow(D$pairs) ) {
          pair<-D$pairs[j,]
          for ( i in 1:nrow(D$confs) ) {
            conf<-D$confs[i,]
            upper<-lower<-rep(0,2)
            upper[conf==0]<-Inf
            lower[conf==1]<-(-Inf)
            #print(dim(D$counts))
            #print(c(u,j,i))
            D$counts[u,j,i]<-pmvnorm(upper=upper,lower=lower,mean=mean[pair],sigma=var[pair,pair],algorithm=pmvnorm.algorithm)
          }
        }
        D$type='infinite'
      } else { #not pairwiseD$
        for ( i in 1:nrow(D$confs) ) {
          conf<-D$confs[i,]
          upper<-lower<-rep(0,n)
          upper[conf==0]<-Inf
          lower[conf==1]<-(-Inf)
          D$counts[u,i]<-pmvnorm(upper=upper,lower=lower,mean=mean,sigma=var,algorithm=pmvnorm.algorithm)
        }
      }
    }
  }

  D$M<-list(A=A,mu=mu,sigma=sigma)
  D$M$p<-blica_M2p(D$M)
  D
}
