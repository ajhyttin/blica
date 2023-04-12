#this file runs basic tests and provides reference on how to use blica.

test1<-function() {
  cat('Test1: Infinite data (true distrubutions as input) n=nsources.\n')
  #Basic run that creates appropriate data and runs the Blica algorithm.
  cat('Creating data.\n')
  D<-blica_createdata(n=10,u=10,nsources=10,N=Inf)
  
  cat('Running Blica with LBFGS...\n')
  stic<-proc.time()[3];
  R<-blica(D)
  t<-proc.time()[3]-stic
  cat('Blica done, took:',t,'\n')
}

test2<-function() {
  cat('Test2: Infinite data (true distrubutions as input) n >> nsources.\n')
  #Basic run that creates appropriate data and runs the Blica algorithm.
  cat('Creating data.\n')
  D<-blica_createdata(n=10,u=10,nsources=5,N=Inf)
  
  cat('Running Blica with LBFGS...\n')
  stic<-proc.time()[3];
  R<-blica(D)
  t<-proc.time()[3]-stic
  cat('Blica done, took:',t,'\n')
}

test3<-function() {
  cat('Test3: Sample data. Same number of samples per segments.\n')
  #Basic test that runs blica from data in a csv file.
  #Basic run that creates appropriate data and runs the Blica algorithm.
  cat('Creating data.\n')
  D<-blica_createdata(n=10,u=20,nsources=5,N=1000)
  
  cat('Running Blica with LBFGS...\n')
  stic<-proc.time()[3];
  R<-blica(D)
  t<-proc.time()[3]-stic
  cat('Blica done, took:',t,'\n')  
}

test4<-function() {
  cat('Test3: Sample data. Different number of samples per segments.\n')
  #Basic test that runs blica from data in a csv file.
  #Basic run that creates appropriate data and runs the Blica algorithm.
  cat('Creating data.\n')
  D<-blica_createdata(n=10,u=10,nsources=5,N=(1:10)*1000)
  
  cat('Running Blica with LBFGS...\n')
  stic<-proc.time()[3];
  R<-blica(D)
  t<-proc.time()[3]-stic
  cat('Blica done, took:',t,'\n')  
}

test5<-function() { #writing and loading
  cat('Test5: Saving and loading data from csv file and running inference.\n')
  #Basic test that runs blica from data in a csv file.
  #Basic run that creates appropriate data and runs the Blica algorithm.
  cat('Creating data.\n')
  D<-blica_createdata(n=10,u=10,nsources=5,N=(1:10)*1000)
  
  cat('Writing data to a file.\n') 
  blica_writecsv(D,file="blica_data.csv")
  cat('Reading data from a file.\n')
  DR<-blica_readcsv(file="blica_data.csv")
  
  cat('Running Blica with LBFGS...\n')
  stic<-proc.time()[3];
  R<-blica(DR, nsources=5 ) #recall to give the number of sources here since it is not 
                            #written in the data files and can thus be different
  t<-proc.time()[3]-stic
  cat('Blica done, took:',t,'\n')  
  mcs<-mcs(D$M$A,R$M$A)
  cat('Final MCS against true mixing matrix:',mcs,'\n')
}



