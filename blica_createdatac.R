blica_createdatac<-function(n=2,u=40,nsources=n,N=Inf,M=NULL,noisevar=0.1) {
  #Creates a continuous data set.
  #Used for example for initializing the optimization.
  #<-1
  if ( any(is.finite(N)) && any(is.infinite(N))) stop('Mixture of infinite and finite sample data not supported.')
  
  if (length(N) == 1) N<-rep(N,u)
  
  #diag(A)<-c(1,1)
  if ( is.null(M) ) {
    A<-array(0,c(n,n))
    A<-array(3*runif(n*n,min=-1,max=1),c(n,nsources))
    #A<-diag(diag(A))
    #A[upper.tri(A)]<-0
    #A[lower.tri(A)]<-0
    #browser()
    mu=array(runif(u*n,min=-2,max=2),c(u,nsources))
    sigma=log(array(runif(u*n,min=0.5,max=3),c(u,nsources)))
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
  R<-array(0,c(length(us),8))
  
  D<-list()

  D$mux<-array(0,c(length(us),n))
  D$sigmax<-array(0,c(length(us),n,n))
#  D$X<- no need to use samples
  for ( ui in us ) {
    D$mux[ui,]<-A%*%mu[ui,]
    D$sigmax[ui,,]<-A%*%diag(exp(sigma[ui,]))%*%t(A)+noisevar*diag(n)
    if ( !is.infinite(N[ui]) )  {
      #in this case just rewrite the sufficient statistics
      browser()
      Z<-rmvnorm(N[ui],mean=D$mux[ui,],sigma=D$sigmax[ui,,])
      D$mux[ui,]<-colMeans(Z)
      D$sigmax[ui,,]<-cov(Z)*(N[ui]-1)/N[ui] #correcting here to finally get sufficient statistics
      D$X<-rbind(D$X,cbind(Z,ui))
    }
  }
  q=array(0,c(u,n))
  D$infinite<-any(is.infinite(N))
  D$N<-N #this was fixed to 10 earlier when variable sample segments were not supported
  D$continuous=TRUE
  D$M<-list(A=A,mu=mu,sigma=sigma,q=q)
  D$M$p<-blica_M2p(D$M,scale=TRUE)
  D
}