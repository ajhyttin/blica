scaledgaussian_p2Mu<-function(p,n,u,nsources=n,scale=FALSE,mean=TRUE) {
  # Turns a parameter vector p into model defined by a list.
  # p - parameter vector for the Blica model.
  # n - number of observed variables
  # u - number of segments
  # nsources - number of latent sources
  # mean - whether means are included in p.
  # scale - tells whether a scaling terms are included in p.
  us<-1:u
  total_us<-u
  M<-list(sigma=array(0,c(u,nsources)),
          q=array(0,c(u,n)),
          A=array(0,c(n,nsources)) )
  for ( u in us ) { #loop through the different segments
    
    #determine the p for tphis one
    if ( scale && !mean ) {
      pindex<- c( ((u-1)*n+1):(u*n),  #scaling factors
                  total_us*n+ (((u-1)*nsources+1):(u*nsources)), #variances
                  (length(p)-n*nsources+1):length(p)) #a matrix  
      
    } else if ( scale && mean ) {
      pindex<- c( ((u-1)*nsources+1):(u*nsources),  #means
                  total_us*nsources+ (((u-1)*nsources+1):(u*nsources)), #variances
                  2*total_us*nsources+ (((u-1)*n+1):(u*n)), #scales
                  (length(p)-n*nsources+1):length(p)) #the mixing matrix
      
    } else {
      pindex<- c( ((u-1)*nsources+1):(u*nsources),  #means
                  total_us*nsources+ (((u-1)*nsources+1):(u*nsources)), #variances
                  (length(p)-n*nsources+1):length(p)) #a matrix
    }
    pu<-p[pindex] #parameters for this particular segment
    Mu<-scaledgaussian_p2M(pu,n,nsources=nsources,scale=scale,mean=mean)

    M$sigma[u,]<-diag(Mu$sigma)
    M$q[u,]<-diag(Mu$q)
    M$A<-Mu$A

  }
  M
}


scaledgaussian_M2p<-function(M,scale=FALSE,mean=FALSE) {
  #Turns a model into a parameter vector.
  # M - model as a list.
  # n - number of observed variables
  # u - number of segments
  # nsources - number of latent sources
  # mean - whether means are included in p.
  # scale - tells whether a scaling terms are included in p.
  
  if ( scale && !mean ) {
    c(t(M$q),t(M$sigma),as.vector(M$A))
  } else if (scale && mean ) {
    c(t(M$mu),t(M$sigma),t(M$q),as.vector(M$A))
  } else {
    c(t(M$mu),t(M$sigma),as.vector(M$A))
  }
}

scaledgaussian_p2M<-function(p,n,nsources=n,scale=FALSE,mean=TRUE) {
  M<-list()
  if ( mean ) M$mu<-p[1:nsources]

  M$sigma<-array(0,c(nsources,nsources))
  diag(M$sigma)<-p[(nsources+1):(2*nsources)]
  
  if ( scale && mean ) {
    M$q<-array(0,c(n,n))
    diag(M$q)<-p[(2*nsources+1):(2*nsources+n)]
  } else if ( scale && !mean ) {
    #browser()
    M$q<-array(0,c(n,n))
    diag(M$q)<-p[1:n]   #had to change this from nsources to n when mean=FALSE
    M$sigma<-array(0,c(nsources,nsources))
    diag(M$sigma)<-p[(n+1):(n+nsources)]  
    #browser()
  }
  M$A<-array(p[(length(p)-n*nsources+1):length(p)],c(n,nsources))
  M
}



scaledgaussian_likelihood<-function(p,D,gradient=TRUE,verbose=FALSE,us=1:nrow(D$M$mu),scale=FALSE,noisevar=8/pi,mean) {
  # Scaled gaussian log-likelihood 
  # p - parameter vector of the model
  # D - continuous data set
  # gradient - whether gradient should be returned
  # verbose - print information or not
  # us - indexes of the segments
  # scale - scaling terms included or not
  # noisevar - the used variance (default=8/pi, this is due to coefficient sqrt(pi/8))
  # mean - whether mean is included
  n<-nrow(D$M$A)
  nsources<-ncol(D$M$A)
  total_us<-nrow(D$M$mu)
  grad<-rep(0,length(p))
  l<-0
  for ( u in us ) { #loop through the different segments
    
    #determine the p for this segment
    if ( scale && !mean ) {
      pindex<- c( ((u-1)*n+1):(u*n),  #scaling factors
                  total_us*n+ (((u-1)*nsources+1):(u*nsources)), #variances
                  (length(p)-n*nsources+1):length(p)) #a matrix  
      
    } else if ( scale && mean ) {
      pindex<- c( ((u-1)*nsources+1):(u*nsources),  #means
                  total_us*nsources+ (((u-1)*nsources+1):(u*nsources)), #variances
                  2*total_us*nsources+ (((u-1)*n+1):(u*n)), #scales
                  (length(p)-n*nsources+1):length(p)) #the mixing matrix
      
    } else {
      pindex<- c( ((u-1)*nsources+1):(u*nsources),  #means
                  total_us*nsources+ (((u-1)*nsources+1):(u*nsources)), #variances
                  (length(p)-n*nsources+1):length(p)) #a matrix
    }
    pu<-p[pindex] #parameters for this particular segment
    
    mux<-D$mux[u,]  #data has D$mux[u,] as mean of x and D$sigmax[u,,] as the covariance/correlation matrix
    if (!mean) mux=NA
    
    ll<-scaledgaussian_likelihood_p(pu,mux=mux,sigmax=D$sigmax[u,,],N=D$N,gradient=gradient,verbose=verbose,X=NULL,nsources,scale=scale,noisevar=noisevar,mean=mean)
    #cat('u:',u,'ll:',ll,'\n')
    if ( is.na(ll) ) return(NA)
    
    l<-l+ll  #here could apply the weighting of different data sets
    if (gradient) grad[pindex]<-grad[pindex]+attr(ll,'gradient')
    
  }
  if ( gradient ) attr(l,'gradient')<-grad
  l
}

scaledgaussian_likelihood_p<-function(p,mux=colMeans(X),sigmax=(N-1)/N*cov(X),N=nrow(X),gradient=TRUE,verbose=FALSE,X=NULL,nsources=n,noisevar=0.1,scale=FALSE,mean=TRUE) {
  #Scaled Gaussian log-likelihood for one segment
  #The gradient calculation is based on formula derivation using the Matrix Cookbook,
  #checked against numerical derivative estimates in development phase.
  n<-nrow(sigmax)
  M<-scaledgaussian_p2M(p,n,nsources,scale=scale,mean=mean)
  A<-M$A
  if (mean) mu_u<-M$mu
  sigma_u<-diag(exp(diag(M$sigma)))

  if ( scale ) {
    q_u<-diag(exp(diag(M$q)))
  }

  if ( mean ) mu_y<-A%*%mu_u  
  sigma_y<-A%*%sigma_u%*%t(A)+noisevar*diag(n)
  
  if ( scale ) {
    sigma_y<-q_u%*%sigma_y%*%q_u
    if (mean) mu_y<-q_u%*%mu_y
  }

  if (!mean) mu_y=NA
  l<-scaledgaussian_likelihood_plain(mu_y,sigma_y,mux=mux,sigmax=sigmax,N=N,gradient=gradient,verbose=FALSE,mean=mean)
  if ( is.na(l)) return(NA)

  if (gradient) {
    sigmag<-attr(l,'sigmag')  #sigmag is symmetric
    mug<-attr(l,'mug')
    
    if ( mean ) {
      g1<-as.vector(t(A)%*%mug)
      if ( scale ) g1<-t(A)%*%q_u%*%mug
    }    
    g2<-diag(t(A)%*%(sigmag+t(sigmag))%*%A)/2*diag(exp(M$sigma))
    if ( scale )  g2<-diag(t(A)%*%q_u%*%sigmag%*%q_u%*%A)*diag(exp(M$sigma))

    if ( scale ) {
      gq1<-2*diag( q_u%*%A%*%sigma_u%*%t(A)%*%q_u%*%t(sigmag) )
      gq2<-diag(2*noisevar*q_u*q_u)*diag(sigmag)
      if (mean) {
        gq3<-as.vector(diag(q_u)*mug*A%*%mu_u)
      } else {
        gq3<-0
      } 
      gq<-gq1+gq2+gq3
    } else {
      gq<-c()
    }
    
    if ( mean ) {
      g3<-as.vector(t(t(mug%*%t(mu_u)+sigmag%*%A%*%sigma_u+t(sigmag)%*%A%*%sigma_u)))
      if (scale ) g3<-as.vector(t(t(q_u%*%mug%*%t(mu_u)+q_u%*%sigmag%*%q_u%*%A%*%sigma_u+q_u%*%t(sigmag)%*%q_u%*%A%*%sigma_u)))
    } else {
      g3<-as.vector(t(t(sigmag%*%A%*%sigma_u+t(sigmag)%*%A%*%sigma_u)))
      if (scale ) g3<-as.vector(t(t(q_u%*%sigmag%*%q_u%*%A%*%sigma_u+q_u%*%t(sigmag)%*%q_u%*%A%*%sigma_u)))
    }
    
    if ( !mean && scale ) {
      attr(l,'gradient')<-c(gq,g2,g3)
    } else {
      attr(l,'gradient')<-c(g1,g2,gq,g3)      
    }
    
    attr(l,'mug')<-NULL
    attr(l,'sigmag')<-NULL
    
  }
  l
}


scaledgaussian_likelihood_plain<-function(mu,sigma,mux=colMeans(X),sigmax=(N-1)/N*cov(X),N=nrow(X),
                                            gradient=TRUE,verbose=FALSE,X=NULL,mean=TRUE) {

  n<-length(mu)
  if (mean) mu<-as.vector(mu)

  
  if (nrow(sigma) != ncol(sigma) ) sigma<-diag(as.vector(sigma))

  #N<-nrow(X)
  
  #isigma<-solve(sigma)
  isigma<-mpinv(sigma)
  #l<-sum(dmvnorm(X,mean=mu,sigma=sigma,log=TRUE))
  
  #alternative way of calculating the above quantity, this is useful for obtaining the derivatives
  if ( mean) {
    corr<-t((mux-mu))%*%isigma%*%(mux-mu)
  } else {
    corr<-0
  }
  

  if ( nrow(sigma) < 100 ) {
   detsigma<-det(sigma)
   #if( detsigma < exp(-20)) cat('Reg sigma.\n')
   detsigma<-max(detsigma,exp(-20))
   logdetsigma<-log(detsigma)
  } else {
    logdetsigma<-determinant(sigma,logarithm=TRUE)$modulus
  }

  l<-as.vector( (-N/2)*(n*log(2*pi)+logdetsigma+sum(diag(sigmax%*%isigma))+ corr ) )
  
  if (gradient) {
    if (mean) attr(l,'mug')<-(-N/2)*(-2*isigma%*%(mux-mu))
    #-diag(diag(isigma))
    if ( mean) {
      attr(l,'sigmag')<-(-N/2)*(isigma-isigma%*%sigmax%*%isigma-isigma%*%(mux-mu)%*%t(mux-mu)%*%isigma)
    } else {
      attr(l,'sigmag')<-(-N/2)*(isigma-isigma%*%sigmax%*%isigma)
    }
  }
  l
}
