blica<-function(D,type='pairwise',lseed=NA,cor_method='o',verbose=TRUE,reg=NA,
                rfile=NULL, n=NA, nsources=NA) {
  #D is the dataset. See for how it is created.
  #type 
  #lseed
  #cor_method
  #verbose - 
  #reg - regularization parameter, use 0 for infinite data or 0.1 for sample data
  #rfile - file stub into which results are saved
  #n - number of observed variables (default: determined from the input)
  #nsources - number of latent sources (default: determined from the input if possible,
  #           otherwise need to specify)  
  
  if ( any(is.finite(D$N)) && any(is.infinite(D$N))) stop('Mixture of infinite and finite sample data not supported.')
  
  
  true_model_given = !is.null(D$M)
  
  #first step
  if ( is.na(reg) && D$type == "sample" ) reg<-10 #regularization needed only for sample data
  if ( is.na(reg) && D$type == "infinite" ) reg<-NA
  if ( verbose ) cat('reg =',reg,' (NA means no regularization of estimated correlation matrices.)\n')

  if ( is.na(n) && D$type == "sample" ) n<-ncol(D$X[[1]])
  if ( is.na(n) && D$type == "infinite" ) n<-max(D$pairs)
  if ( verbose ) cat('n =',n,'\n')
  
  if (is.na(nsources) && true_model_given ) nsources<-ncol(D$M$A) #number of sources can be determined only 
                                                              #if the true model is given
  if ( is.na(nsources)) stop('Please give the number of sources as input blica(D,nsources=nsources).')
  if ( verbose ) cat('nsources =',nsources,'\n')
  
  
  if (D$type=="sample") u<-length(D$X)
  if (D$type == "infinite" )  u<-dim(D$counts)[1]
  if ( verbose ) cat('segments=',u,'\n')
  
  
  if ( verbose ) cat('PHASE 1: Estimate pairwise correlations from binary data.\n')
  
  if ( true_model_given ) { #these are used for outputs only when true model is present
    truecovs_adj<-truecovs<-array(0,c(nrow(D$M$A),nrow(D$M$A),u))
    omus<-array(0,c(nrow(D$M$A),nrow(D$counts)))
    qs<-trueqs<-array(1,c(nrow(D$counts),nrow(D$M$A)))
  }
  tic()
  index<-0
  #create the continuous data structure into whcih the estimates will be written
  DC<-blica_createdatac(n=n,u=u,nsources=nsources,N=D$N)
  DC$N[is.infinite(DC$N)]<-10 #put here sample size to 10 if infinite sample data

  us = 1:u
  for ( ui in us ) {
    
    if ( D$type != "sample" && true_model_given) {
      covtru<-diag(nrow(D$M$A))+pi/8*D$M$A%*%diag(exp(D$M$sigma[ui,]))%*%t(D$M$A)   
      truecovs[,,ui]<-covtru
      truecovs_adj[,,ui]<-(8/pi* (truecovs[,,ui]-diag(nrow(D$M$A))))
    }

    covest<-diag(rep(1,n))
    
    for ( i in 1:n ) {
      for ( j in 1:n ) {
        if ( i >= j ) next; #estimate each correlation only once
        
        #Calculate the pairwise distribution if not given as input or fetch it.
        if (D$type == "sample" ) {
          counts<-rep(NA,4)
          for ( k in 1:4 ) {
            conf<-dec.to.bin(k-1,2)
            counts[k]<-sum((D$X[[ui]][,i]==conf[1]) & (D$X[[ui]][,j]==conf[2]) )
          } 
        } else if (D$pairwise) {
          ii<-which(D$pairs[,1] == i & D$pairs[,2] == j)
          counts<-D$counts[ui,ii,]
        } else {
          counts<-rep(NA,4)
            for ( k in 1:4 ) {
              conf<-dec.to.bin(k-1,2)
              I<-(D$confs[,i]==conf[1]) & (D$confs[,j]==conf[2])
              counts[k]<-sum(D$counts[ui,I])
            } 
        }

        R<-blica_correlation(counts,cor_method)
        
        covest[i,j]<-covest[j,i]<-R
      } #for j
    } #for i
    
    if ( verbose ) cat('segment:',ui,'/',u,'\n')
    if (D$type != "sample"  && true_model_given ) {
      error<-max(abs(cov2cor(covest)-cov2cor(covtru)))
      if ( verbose ) cat(' max cor error before regularization:',error,'\n')
    }

    evs<-eigen(covest)$values
    if ( verbose ) cat('min eigenvalue before regularization:',min(evs),'\n')      
    
    if ( !is.na(reg) ) { 
      deltai<-max(c( 0 , (max(evs)-reg*min(evs))/(reg-1) ) )
  
          
      covest<-1/(1+deltai)*(covest+deltai*diag(n))
      evs<-eigen(covest)$values
      if ( verbose ) cat('min eigenvalue after regularization:',min(evs),'\n')      
  
      if (D$type != "sample" && true_model_given) {
        error<-max(abs(cov2cor(covest)-cov2cor(covtru)))
        if ( verbose ) cat(' max cor error after regularization:',error,'\n')
      }
    } else {
      cat('No regularization done.\n')
    }
    if (any(evs < 0) ) { #some regularization needed
      stop('Negative eigenvalues (after regularization). Use a higer regularization parameter, e.g., current+10.')
    }

    DC$sigmax[ui,,]<-covest
    DC$mux[ui,]<-0 #setting means to zero!
    if ( verbose ) cat('t:',toc(),'\n')
  }
  tphase1<-toc()
  if ( verbose ) cat('Took:',tphase1,'\n')
  if ( verbose ) cat('Time after PHASE 1:',tphase1,'\n')

  if ( verbose ) cat('PHASE 2: Running LBFGS on scaled Gaussian likelihood and continuous data:\n')
  if (!is.na(lseed) ) {
    set.seed(lseed);
  }
  noisevar=8/pi;
  if ( verbose ) cat('Noise level:',noisevar,'\n')
  Mstart<-list(mu=array(0,dim(DC$M$mu)),sigma=log(array(runif(prod(dim(DC$M$sigma)),0.5,1.5),dim(DC$M$sigma))),
                 q=array(0,c(nrow(DC$M$mu),nrow(DC$M$A))), #segments times sources
                 A=array(3*runif(n*nsources,min=-1,max=1),c(n,nsources))
               )
  mean=FALSE
  #"browser()
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
      mcs<-NA
      if ( true_model_given) mcs<-mcs(D$M$A,Aest)
      kappa<-kappa(Aest)
      cat('t:',lasttoc1,'l:',format(fval,digits=15),'kappa:',kappa,'mcs:',mcs,'\n');
      #print(round(Aest,2));
      cat('Variance summary:\n');
      print(summary(exp(x[(u*n+1):(u*n+u*nsources)])));

      if ( !is.null(rfile) ) {
        write.table(Aest,file=paste(rfile,'.mix',sep=''),row.names=FALSE,col.names=FALSE)
        cat(lasttoc1,format(fval,digits=15),mcs,'\n',file=rfile,append=append);
        append<<-TRUE; 
      }
    }
    fval
  } 
  
  #gradient
  g<-function(x) {
    v<-(-1)*attr(scaledgaussian_likelihood(x,DC,gradient=TRUE,noisevar=noisevar,scale=TRUE,mean=mean),'gradient')
    if ( verbose && (abs(toc()-lasttoc2) > interval)  ) { 
      lasttoc2<<-toc()
      cat('t:',lasttoc2,'max abs grad:',format(max(abs(v)),digits=15),'\n')
    }
    v
  }
  
  R<-optim(start,f,g,method="L-BFGS-B",control=list(maxit=100000,factr=0,pgtol=1e-10))
  if ( verbose ) cat('Done!\n')
  if ( verbose ) print(R)
  tphase2<-toc()
  if ( verbose ) cat('Took:',tphase2-tphase1,'\n')
  if ( verbose ) cat('Time after PHASE 2:',tphase2,'\n')
  #browser()
  lasttoc1<<-lasttoc2<<-(-20)
  f(R$par);g(R$par);
  
  x<-R$par;
  Aest<-array(x[(length(x)-n*nsources+1):length(x)],c(n,nsources))
  if ( !is.null(rfile) ) write.table(Aest,file=paste(rfile,'.mix',sep=''),row.names=FALSE,col.names=FALSE)

  R$l<-scaledgaussian_likelihood(R$par,DC,gradient=FALSE,noisevar=noisevar,scale=TRUE,mean=mean)
  if ( verbose ) cat('Final likelihood value obtained:',R$l,'\n')
  #browser()
  R$M<-scaledgaussian_p2Mu(R$par,nrow(DC$M$A),nrow(DC$M$mu),ncol(DC$M$A),scale=TRUE,mean=mean)
  invisible(R)  
}
