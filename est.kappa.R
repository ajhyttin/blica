est.kappa<-function(n,nsources=n,times=100) {
  kappas<-rep(0,times)
  for ( i in 1:times ) {
    A<-array(runif(n*nsources,min=-3,max=3),c(n,nsources))
    kappas[i]<-kappa(A)
  }
  summary(kappas)
  #hist(kappas)
  quantile(kappas,0.75)
} 