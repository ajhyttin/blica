blica<-function(D,type='pairwise',lseed=NA,cor_method='o',verbose=TRUE,reg=0,
                skip=FALSE, rfile=NULL, delta=0, corfile=NA) {
  #D is the dataset. See for how it is created.
  #type 
  #lseed
  #cor_method
  #verbose - 
  #reg - regularization parameter, use 0 for infinite data or 0.1 for sample data
  #skip -
  #rfile -
  #delta - 
  #corfile -
  
  #first step
  #delta<-0 #0.1 #0.1
  cat('reg=',reg,'\n')
  n<-nrow(D$M$A)
  nsources<-ncol(D$M$A)
  if (D$type=="sample") {
    us<-length(D$X)
  } else {
    us<-nrow(D$M$mu)
  } 
  
  if ( is.na(corfile) ) {
  cat('PHASE 1: Find covariances.\n')
  truecovs_adj<-truecovs<-array(0,c(nrow(D$M$A),nrow(D$M$A),us))
  
  omus<-array(0,c(nrow(D$M$A),nrow(D$counts)))
  qs<-trueqs<-array(1,c(nrow(D$counts),nrow(D$M$A)))
  cat('Estimating correlation matrices:\n')
  tic()
  index<-0
  #continuous data to be used later
  DC<-blica_createdatac(n=n,u=us,nsources=nsources)
  DC$M$A<-D$M$A
  for ( u in 1:us ) {
    
    if ( D$type != "sample" ) {
      covtru<-diag(nrow(D$M$A))+pi/8*D$M$A%*%diag(exp(D$M$sigma[u,]))%*%t(D$M$A)   
      truecovs[,,u]<-covtru
      truecovs_adj[,,u]<-(8/pi* (truecovs[,,u]-diag(nrow(D$M$A))))
    }

    covest<-diag(rep(1,n))
    
    for ( i in 1:n ) {
      for ( j in 1:n ) {
        if ( i >= j ) next;
        if (D$type == "sample" ) {
          counts<-rep(NA,4)
          for ( k in 1:4 ) {
            conf<-dec.to.bin(k-1,2)
            counts[k]<-sum((D$X[[u]][,i]==conf[1]) & (D$X[[u]][,j]==conf[2]) )
          } 
         # print(counts)
        } else if (D$pairwise) {
          ii<-which(D$pairs[,1] == i & D$pairs[,2] == j)
          counts<-D$counts[u,ii,]
        } else {
          counts<-rep(NA,4)
            for ( k in 1:4 ) {
              conf<-dec.to.bin(k-1,2)
              I<-(D$confs[,i]==conf[1]) & (D$confs[,j]==conf[2])
              counts[k]<-sum(D$counts[u,I])
            } 
        }
        #if ( reg == 1 ) counts<-counts+1
        #if (any(counts== 0) && skip ) next
        #if ( any(counts == 0) ) counts<-counts+1
        R<-blica_correlation(counts,cor_method)
        #if ( reg == 2 ) covest<-covest+diag(n)
        #if ( reg == 3 ) covest<-mpinv(mpinv(covest+diag(n)))
        covest[i,j]<-covest[j,i]<-R
      } #for j
    } #for i
    
    cat('segment:',u,'/',us,'\n')
    if (D$type != "sample" ) {
      error<-max(abs(cov2cor(covest)-cov2cor(covtru)))
      cat(' max cor error:',error,'\n')
      #browser()
    }
    #while (TRUE ) {
      #covest<-1/(1+delta)*(covest+delta*diag(n))
      #while ( TRUE ) {
      evs<-eigen(covest)$values
      cat('before eigenvalues:',evs,'\n')      
      deltai<-max(c( 0 , (max(evs)-reg*min(evs))/(reg-1) ) )
      #deltai<-0
      covest<-1/(1+deltai)*(covest+deltai*diag(n))
      evs<-eigen(covest)$values
      cat('after eigenvalues:',evs,'\n')      

      if (D$type != "sample" ) {
        error<-max(abs(cov2cor(covest)-cov2cor(covtru)))
        cat(' max cor error:',error,'\n')
      }
      if (any(evs < 0) ) { #some regularization needed
        stop('negative eigenvalues after regularization.')
    #    browser()
    #    deltai<-max(c(0,(max(evs)-100*min(evs))/99))
    #    covest<-1/(1+deltai)*(covest+deltai*diag(n))
        #covest<-covest +(abs(min(evs))+0.1)*diag(n)
    #    evs<-eigen(covest)$values
    #    cat(' eigenvalues:',evs,'\n')
      }
    #}
    #print(covest)
    DC$sigmax[u,,]<-covest
    print(round(covest,2));
    DC$mux[u,]<-0 #setting means to zero!
    cat('t:',toc(),'\n')
  } # if ( is.na(corfile) )
  tphase1<-toc()
  cat('Took:',tphase1,'\n')
  cat('Time after PHASE 1:',tphase1,'\n')

#  cat('saving covs.\n')
#  save(DC,file=paste(rfile,'.cor',sep=''))

  } #else {
    
  cat('PHASE 2: Running LBFGS on scaled Gaussian likelihood:\n')
  if (!is.na(lseed) ) {
    set.seed(lseed);
  }
  noisevar=8/pi;
  cat('Noise level:',noisevar,'\n')

  Mstart<-list(mu=array(0,dim(DC$M$mu)),sigma=log(array(runif(prod(dim(DC$M$sigma)),0.5,1.5),dim(DC$M$sigma))),
                 q=array(0,c(nrow(DC$M$mu),nrow(DC$M$A))),
                 A=array(3*runif(n*nsources,min=-1,max=1),c(n,nsources))
               )
  mean=FALSE
  #browser()
  start<-scaledgaussian_M2p(Mstart,scale=TRUE,mean=mean)
  f_upper<-0 #blica_upper(DC)
  lasttoc1<<-lasttoc2<<-(-20)
  interval<-3
  append<-FALSE
  rfile<<-rfile;
  #function
  f<-function(x) {
    fval<-f_upper-scaledgaussian_likelihood(x,DC,gradient=FALSE,noisevar=noisevar,scale=TRUE,mean=mean)

    Aest<-array(x[(length(x)-n*nsources+1):length(x)],c(n,nsources))
    if (verbose && (abs(toc()-lasttoc1) > interval) ) { 
      lasttoc1<<-toc()
      mcs<-mcs(D$M$A,Aest)
      kappa<-kappa(Aest)
      cat('t:',lasttoc1,'l:',format(fval,digits=15),'kappa:',kappa,'mcs:',mcs,'\n');
      print(round(Aest,2));
      cat('Variance summary:\n');
      print(summary(exp(x[(u*n+1):(u*n+u*nsources)])));
      #browser()
      
      
      if ( !is.null(rfile) ) {
        if ( !is.null(rfile) ) write.table(Aest,file=paste(rfile,'.mix',sep=''),row.names=FALSE,col.names=FALSE)
        cat(lasttoc1,format(fval,digits=15),mcs,'\n',file=rfile,append=append);
        append<<-TRUE; 
      }
    }
    fval
  } 
  
  
  #gradient
  g<-function( x) {
    v<-(-1)*attr(scaledgaussian_likelihood(x,DC,gradient=TRUE,noisevar=noisevar,scale=TRUE,mean=mean),'gradient')
    if ( verbose && (abs(toc()-lasttoc2) > interval)  ) { 
      lasttoc2<<-toc()
      cat('t:',lasttoc2,'max abs grad:',format(max(abs(v)),digits=15),'\n')
    }
    v
  }

  #browser()
  
  R<-optim(start,f,g,method="L-BFGS-B",control=list(maxit=100000,factr=0,pgtol=1e-10))#list(trace=6,REPORT=1)
  cat('Done!\n')
  print(R)
  tphase2<-toc()
  cat('Took:',tphase2-tphase1,'\n')
  cat('Time after PHASE 2:',tphase2,'\n')
  #browser()
  lasttoc1<<-lasttoc2<<-(-20)
  f(R$par);g(R$par);
  
  x<-R$par;
  Aest<-array(x[(length(x)-n*nsources+1):length(x)],c(n,nsources))
  if ( !is.null(rfile) ) write.table(Aest,file=paste(rfile,'.mix',sep=''),row.names=FALSE,col.names=FALSE)

  R$l<-scaledgaussian_likelihood(R$par,DC,gradient=FALSE,noisevar=noisevar,scale=TRUE,mean=mean)
  if ( verbose ) cat('Final likelihood obtained:',R$l,'\n')
  #browser()
  R$M<-scaledgaussian_p2Mu(R$par,nrow(DC$M$A),nrow(DC$M$mu),ncol(DC$M$A),scale=TRUE,mean=mean)
  invisible(R)  
}
