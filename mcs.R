mcs<-function(A1,A2,verbose=FALSE) {
  cost<-abs(crossprod(A1,A2))/(sqrt(apply(A1^2,2,sum))%*%t(sqrt(apply(A2^2,2,sum))))
  R<-lp.assign(cost,direction="max")
  R$objval/ncol(A1)
}


