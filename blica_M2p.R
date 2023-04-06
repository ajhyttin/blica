blica_M2p<-function(M,scale=FALSE) {
  # Turns model parameters FOR A SINGLE SEGMENT - define by a list M with 
  # source means M$mu, source variances M$sigma, and mixing matrix M$A
  # into a vector such that it can be optimized with general purpose methods.
  # if scale=TRUE also scale terms are added.
  if ( scale ) {
    c(t(M$mu),t(M$sigma),t(M$q),as.vector(M$A))
  } else {
    c(t(M$mu),t(M$sigma),as.vector(M$A))
  }
}

blica_p2M<-function(p,n,nsources=n,scale=FALSE) {
  # Turns a parameter vector into a model defined as a list FOR A SINGLE SEGMENT.
  # p - the parameter vector
  # n - the number of observed variables
  # nsources - number of source variables
  # scale - whether scale parameters should be included or not
  M<-list()
  M$mu<-p[1:nsources]
  M$sigma<-array(0,c(nsources,nsources))

  diag(M$sigma)<-p[(nsources+1):(2*nsources)]
  
  if ( scale ) {
    M$q<-array(0,c(n,n))
    diag(M$q)<-p[(2*nsources+1):(2*nsources+n)]
  }

  M$A<-array(p[(length(p)-n*nsources+1):length(p)],c(n,nsources))
  M
}

blica_p2Mu<-function(p,n,u,nsources=n,scale=FALSE) {
  # Turns a parameter vector into a model defined as a list FOR ALL SEGMENTS.
  # p - the parameter vector
  # n - the number of observed variables
  # nsources - number of source variables
  # scale - whether scale parameters should be included or not
  us<-u
  M<-list()
  M$mu<-array(NA,c(us,nsources))
  M$sigma<-array(NA,c(us,nsources))
  if ( scale) M$q<-array(NA,c(us,n))
  M$A<-array(NA,c(n,nsources))
  for ( u in 1:us ) {
    
    #determine the p for tphis one
    if ( scale ) {
      pindex<- c( ((u-1)*nsources+1):(u*nsources),  #means
                  total_us*nsources+ (((u-1)*nsources+1):(u*nsources)), #variances
                  2*total_us*nsources+ (((u-1)*n+1):(u*n)), #scales
                  (length(p)-n*nsources+1):length(p)) #the mixing matrix
      
    } else {
      pindex<- c( ((u-1)*nsources+1):(u*nsources),  #means
                  total_us*nsources+ (((u-1)*nsources+1):(u*nsources)), #variances
                  (length(p)-n*nsources+1):length(p)) #a matrix
    }
    #determine the p for this one
  #  pindex<- c( ((u-1)*nsources+1):(u*nsources),us*nsources+ (((u-1)*nsources+1):(u*nsources)) , (length(p)-n*nsources+1):length(p))
  #  if ( scale ) 
    #print(pindex)
    #browser()
    #the A matrix is last n*n elements
    pu<-p[pindex]
    Mu<-blica_p2M(pu,n,nsources,scale=scale)
    M$mu[u,]<-Mu$mu
    M$sigma[u,]<-diag(Mu$sigma)
    if ( scale ) M$q[u,]<-diag(Mu$q)
    M$A<-Mu$A
  }

  M
}
