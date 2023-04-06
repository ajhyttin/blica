#this file runs basic tests and provides reference on how to use blica.

test1<-function() {
  #Basic run that creates appropriate data and runs the Blica algorithm.
  cat('Creating data.\n')
  D<-blica_createdata(n=10,u=20,nsources=5,N=Inf)
  
  cat('Running Blica with LBFGS...\n')
  stic<-proc.time()[3];
  R<-blica(D)
  t<-proc.time()[3]-stic
  cat('Blica done, took:',t,'\n')

}

test2<-function() {
  #Basic test that runs blica from data in a csv file.
  
  
  
  
}

