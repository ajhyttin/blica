blica_writecsv<-function(D,file='blica_data.csv') {
  #writes data created by blica_createdata to csv file format
  if (any(is.infinite(D$N))) stop('Writing infinite data to csv file not possible/supported.')
  for ( ui in 1:length(D$X)) {
    write.table(cbind(D$X[[ui]],ui),file=file,append=(ui != 1),row.names=FALSE,col.names=FALSE,sep=',')
  }  
}


blica_readcsv<-function(file='blica_data.csv') {
  #reads a csv file to the list format used by blica
  D<-list()
  D$X<-list()
  #D$X is a list where D$X[[ui]] are the samples for segment ui
  D$pairwise<-FALSE
  D$infinite<-FALSE
  D$type<-"sample"
  X<-read.table(file=file,header=FALSE,sep=',')
  u<-max(X[ncol(X)])
  D$N<-rep(NA,u)
  for (ui in 1:u ) {
    D$X[[ui]]<-X[(X[,ncol(X)]==ui),1:(ncol(X)-1)]
    D$N[ui]<-nrow(D$X[[ui]])
  }
  D$M<-NULL
  D
}