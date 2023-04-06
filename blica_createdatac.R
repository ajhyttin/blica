blica_createdatac<-function(n=2,u=40,nsources=n,N=Inf,M=NULL,noisevar=0.1) {
  #Creates a continuous data set.
  #Used for example for initializing the optimization.
  #<-1
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
  for ( u in us ) {
    D$mux[u,]<-A%*%mu[u,]
    D$sigmax[u,,]<-A%*%diag(exp(sigma[u,]))%*%t(A)+noisevar*diag(n)
    if ( !is.infinite(N) )  {
      #in this case just rewrite the sufficient statistics
      Z<-rmvnorm(N,mean=D$mux[u,],sigma=D$sigmax[u,,])
      D$mux[u,]<-colMeans(Z)
      D$sigmax[u,,]<-cov(Z)*(N-1)/N #correcting here to finally get sufficient statistics
      D$X<-rbind(D$X,cbind(Z,u))
    }
  }
  q=array(0,c(u,n))
  D$infinite<-is.infinite(N)
  D$N<-10
  D$continuous=TRUE
  D$M<-list(A=A,mu=mu,sigma=sigma,q=q)
  D$M$p<-blica_M2p(D$M,scale=TRUE)
  D
}