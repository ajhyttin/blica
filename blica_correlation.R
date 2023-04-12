blica_correlation<-function( counts, method='o',transform=FALSE,verbose=FALSE,cortru=NA) {
  #calculates the pairwise correlation for two binary variables as defined in the paper, i.e.
  #which correlation the variables have before binarization.
  #browser()

  #if ( method!= "weighted") counts<-counts+1
  if ( method=="weighted" && sum(counts) > 2) counts<-counts+1e-2#it seems the weights do not work if counts are 0!
  #cat('Counts:',counts,'\n')
  if ( any(counts==0)) {
    #cat('Counts:',counts,'\n')
    #counts<-counts+1
  }
  if ( method == 'first' ) {
    Dp<-blica_createdata(n=2,u=1,pairwise=FALSE) #now need to update counts
    Dp$counts[1,]<-counts
    R<-blica_lbfgs(Dp,verbose=TRUE)
    covesti<-diag(nrow(Dp$M$A))+pi/8*R$M$A%*%diag(exp(R$M$sigma[1,]))%*%t(R$M$A)
    coresti<-cov2cor(covesti)
    print(coresti[1,2])
    coresti[1,2]
  } else {
    eps<-1e-10
    p1<-(sum(counts[3:4]))/(sum(counts))
    p2<-(sum(counts[c(2,4)]))/(sum(counts))
    if (p1 > 1-eps) p1<-1-eps
    if (p1 < eps) p1<-eps
    if (p2 > 1-eps) p2<-1-eps
    if (p2 < eps) p2<-eps
    mu<-rep(NA,2)  
    mu<-c(qnorm(p1),qnorm(p2))*(-1)
    #print(mu)
    #browser()
    f_upper<-sum(counts*log(counts/sum(counts)))
    #eps<-1e-2
    
    #p0<-0
    alpha0<-p0<-counts2cor(counts)
    if ( !transform ) {
      #cat('NOT transformed.\n')
      #R<-optim(par=p0,fn=blica_correlation_f,gr=blica_correlation_g,method="L-BFGS-B",
      #       control=list(maxit=100000),
      #       upper=0.999,lower=-0.999,
      #       counts=counts,mu=mu,transform=transform,f_upper=0)
      if ( max(counts) > 1 ) { #for sample data
        tol<-1e-2
        #cat('lower')
      } else { #for infinite data
        tol<-1e-5
        #cat('higher')
      }
      interval=c(-1+1e-5,1-1e-5)
    #if (method == "weighted") interval=c(-1,1)
      R<-optimize(f=blica_correlation_f,interval=interval,tol=tol,
                  counts=counts,mu=mu,transform=transform,f_upper=0)
      #  10.198  secs
      R$par<-R$minimum   
      if (method == "weighted") {
        #browser()
        #fpar<-blica_correlation_f(R$par,counts=counts,mu=mu,f_upper=0)
        #lower<-blica_correlation_f(max(c(-1,R$par-0.1)),counts=counts,mu=mu,f_upper=0)
        #upper<-blica_correlation_f(min(c(1,R$par+0.1)),counts=counts,mu=mu,f_upper=0)
        #w<-max(abs(lower-fpar),abs(fpar-upper))
       # print(c(lower,upper))
        #blica_correlation_f(R$par-1e-3,counts=counts,mu=mu,f_upper=0)
        #a<-seq(from=-1+1e-5,to=1-1e-5,length.out=100)
        #a<-c(R$par-0.01,R$par,R$par+0.01)
        a<-sort(R$par+rnorm(100,sd=0.1))
        a<-a[a>-1]
        a<-a[a<1]
        fval<-rep(0,length(a))
        for ( iii in 1:length(a) ) {
          fval[iii]<-blica_correlation_f(a[iii],counts=counts,mu=mu,f_upper=0,transform=transform)-R$objective
        }
        asq<-(a-R$par)^2
        M<-lm(fval~asq-1)
        #browser()
        M$coefficients<-c(R$par^2,-2*R$par,1)*M$coefficients    #(a-R$par)^2 = a^2 - 2*R$par*a+R$par^2
        #instead we could fit the coefficients such that R$par is the optimum
        #c[2]+2*c[3]*a = 0
        # w*(a-R$par)^2 = w(a^2-2*a*R$par+R$par^2)
        #asq2<-(a-R$par)^2  #= a^2 -2aR$par + R$par^2
        #M2<-lm(fval~asq2)
        #M2$coefficients[1]+ TOOK OUT FROM FIRST!
        #M$coefficients<-c(M2$coefficients[2]*R$par^2,
         #                 M2$coefficients[2]*(-2)*R$par,
        #                  M2$coefficients[2])
        #fval<-fval-M$coefficients[1]
        #intercept<-M$coefficients[1]
        #M$coefficients[1]<-0
        #test with unit
        #M$coefficients<-c(R$par^2,(-2)*R$par,1)
        R$par2<-(-1)*M$coefficients[2]/(2*M$coefficients[3])
        if ( !is.na(cortru ) ) {

          loss<-sum(M$coefficients*c(1,cortru,cortru^2))
          error<- abs(R$par-cortru)
          cat('counts:',paste(counts,collapse='-'),
              'cor:', R$par,'cor2:',R$par2,'truloss:',loss,'error:',error,'coef:',M$coefficients[3],'\n')
          #browser()
        }
        # c[1]+c[2]*a+c[3]*a^2
        #     c[2]+2*c[3]*a = 0
        # a = -c[2]/(2*c[3])
        #print(counts)
        R$par2<-(-1)*M$coefficients[2]/(2*M$coefficients[3])
        #print(c(R$par,R$par2, abs(R$par-R$par2) ) )
      # browser()
        if (any(is.na(M$coefficients))) browser()
        return(c(M$coefficients,R$par,R$par2))
      }
      
      
      if ( FALSE ) {
        #browser()
        a<-seq(from=-1+1e-5,to=1-1e-5,by=0.01)
        fval<-rep(0,length(a))
        for ( iii in 1:length(a) ) {
          fval[iii]<-blica_correlation_f(a[iii],counts=counts,mu=mu,f_upper=0)
        }
        plot(a,-fval)
        asq<-a^2
        M<-lm(fval~a+asq)
        fval2<-cbind(1,a,asq)%*%M$coefficients
        points(a,-fval2,col='red')
        browser()
      }
    } else {
      #cat('Transformed.\n')
      p0<-log((p0+1)/2)-log((1-(p0+1)/2))
      R<-optimize(f=blica_correlation_f,interval=c(-15.293,15.293),
               counts=counts,mu=mu,transform=transform,f_upper=0)
      # 11.363 secs
      R$par<-R$minimum
      browser()
    }

    alpha<-R$par
    if ( method == "weighted") attr(alpha,'weight')<-w
    #alpha<-R$minimum
    if (transform ) { #this is wrong we need to transform the other direction

      alpha<-2*exp(R$par)/(1+exp(R$par))-1
    }

    if (verbose ) cat('spearman:',alpha0,'qcor:',alpha,'grad:',df,'\n')
    #browser()
    alpha
  }
}

blica_correlation_f<-function(x,counts,mu,transform=TRUE,f_upper) {
  #x is now the only parameter
  #transform x
  xorig<-x
  if (transform) {
    x<-2*exp(x)/(1+exp(x))-1
  }
  #dalpha/dx = (2*exp(x)*(1+exp(x))-exp(x)*2*exp(x) )/(1+exp(x))^2
  #  = (2*exp(x) )/(1+exp(x))^2
  # (1+exp(x))(alpha+1)<-2*exp(x)
  
  
  sigma<-diag(2)
  sigma[1,2]<-sigma[2,1]<-x
  fval<-correlation_likelihood(mu,sigma,counts,gradient=FALSE)
  #gradient for x is directly the 1/2 element of sigmag
  #cat('xorig',xorig,'x',x,'fval:',fval,'\n')
  f_upper-fval
}

blica_correlation_g<-function(x,counts,mu,verbose=FALSE,transform=TRUE,f_upper=0) {
  xorig<-x
  if (transform) x<-2*exp(x)/(1+exp(x))-1
  sigma<-diag(2)
  sigma[1,2]<-sigma[2,1]<-x
  fval<-correlation_likelihood(mu,sigma,counts,gradient=TRUE)
  if (verbose) print(fval)
  gval<-attr(fval,'sigmag')[1,2]
  if (transform) gval<-gval*(2*exp(x))/(1+exp(x))^2
  #if (verbose) 
  #cat('xorig',xorig,'x',x,'fval:',fval,'gval:',gval,'\n')
  gval*(-1)
}

#Miwa(4098/2)mvnorm.algorithm=GenzBretz()
correlation_likelihood<-function(mu,sigma,counts,gradient=TRUE,verbose=FALSE,pmvnorm.algorithm=Miwa(4098/2)) {
  #print(mu)
  #print(sigma)
  #for one component
  n<-length(mu)
  mu<-as.vector(mu)
  #print(mu)
  #print(sigma)
  if (nrow(sigma) != ncol(sigma) ) sigma<-diag(as.vector(sigma))
  #browser()
  if ( any( eigen(sigma)$values < 0) ) {
    cat('sigma not posdef!!!\n')
    print(sigma)
    print(NA)
    return(NA)
  }
  
  isigma<-solve(sigma)
  
  p<-rep(0,2^n)
  tmeans<-array(0,c(length(counts),n))
  tvars<-array(0,c(length(counts),n,n))
  
  if ( !is.null(pmvnorm.algorithm)) {
    #confs<-
    for ( i in 1:length(counts) ) {
      conf<-dec.to.bin(i-1,2) #no need to carry confs
      #conf<-confs[i,]
      count<-counts[i]      #sum(apply(t(X) == conf,2,all))
      if ( count == 0 ) next #not working if we normalize
      upper<-lower<-rep(0,n)
      upper[conf==0]<-Inf
      lower[conf==1]<-(-Inf)
      
      
      #gather the probabilities, sum up later
      #set.seed(0);
      #browser()
      pp<-pmvnorm(upper=upper,lower=lower,mean=mu,sigma=sigma,algorithm=pmvnorm.algorithm)
      pp<-as.vector(pp)
      p[i]<-pp
      
      if ( gradient ) {
        #set.seed(0);
        moments<-mtmvnorm(mu, sigma, lower=lower, upper=upper,pmvnorm.algorithm =pmvnorm.algorithm)
        
        tmeans[i,]<-moments$tmean
        tvars[i,,]<-moments$tvar
      }
    }
    
    #  browser()
  } else {
    #set.seed(0);
    #browser()
    N<-1000
    S<-rmvnorm(N,mean=mu,sigma=sigma)
    for ( i in 1:length(counts) ) {
      conf<-confs[i,]
      
      upper<-lower<-rep(0,n)
      upper[conf==0]<-Inf
      lower[conf==1]<-(-Inf)
      #
      I<-apply( ( t(S)< upper ) & (t(S) > lower) ,2,all)  #for configuration of all ones
      
      if ( any(I) ) {
        
        pp<-sum(I)/N
        
        p[i]<-pp
        
        if ( gradient ) {
          # browser()
          tmeans[i,]<-colMeans(S[I,,drop=FALSE])
          tvars[i,,]<-cov(S[I,,drop=FALSE])
        }
      } else {
        p[i]<-1e-10
        tmeans[i,]<-rep(0,n)
        tvars[i,,]<-diag(rep(1,n))
      }
    }
  }
  
  
  #browser()
  
  #calculating here: since there is some noise in the probabilities, good to make sure they 
  #actually sum up to 1
  #same should be done for the covariance matrix and mean:
  #on average they should match mean and covariance of the true model.
  if ( any(p <= 0) ) {
    #cat('Watch out negative probabilities?\n')
    #browser()
    p[p<=0]<-1e-20 #for low precision estimates negatives may be producd
    
  }
  p<-p/sum(p)
  
  l<-0
  mugs<-rep(0,n)
  sigmags<-array(0,c(n,n))
  if (gradient ) {
    for ( i in 1:length(counts) ) {
      count<-counts[i]      #sum(apply(t(X) == conf,2,all))
      tmean<-tmeans[i,]
      tvar<-tvars[i,,]
      
      apu<-tmean-mu
      correction<-apu%*%t(apu)
      #browser()
      mug<-as.vector(isigma%*%(tmean-mu))
      #browser()
      sigmag<-(t(isigma)%*%(tvar+correction )%*%t(isigma)-t(isigma))
      
      #ah, this is because sigma is symmetric
      diag(sigmag)<-1/2*diag(sigmag)
      mugs<-mugs+count*mug
      sigmags<-sigmags+count*sigmag
    }
  }#if gradient  
  
  #cat('mugs at plain:\n')
  #print(mugs)
  #print(p)
  l<-l+sum(counts*log(p))
  #browser()
  if (gradient) {
    attr(l,'mug')<-mugs
    attr(l,'sigmag')<-sigmags
  }
  l
}


counts2cor<-function(counts) {
  if (any(counts==0)) counts<-counts+1
  if ( sum(counts) < 2 ) counts<-round(counts*10000)
  D<-array(0,c(sum(counts),2))
  D[(counts[1]+1):(counts[1]+counts[2]),2]<-1
  D[(counts[1]+counts[2]+1):sum(counts),1]<-1
  D[(counts[1]+counts[2]+counts[3]+1):sum(counts),2]<-1
  cor(D,method='spearman')[1,2]
}

