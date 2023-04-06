# Blica
Code for the Blica method published at UAI2022.

Note that the iVAE based method described in the paper can be found in a 
separate repo: 
https://github.com/vitoriapacela/iVAE

Binary Independent Component Analysis: A Non-stationarity-based Approach. Conference on Uncertainty in Artificial Intelligence 2022.
https://arxiv.org/abs/2111.15431


> source('load.R')
> loud()

#create a data set, infinite number of samples
>D<-blica_createdata(n=2,u=5)

#Upper bound on the likelihood is available by
> blica_upper(D)

#run lbfgs:
>blica(D)   

#or you can run gradient descent, uses Miwa() by default
>sense_gd(D)

#With faster but more stochastic alternative for cdf and mvnorm and 
tmvnorm
>sense_gd(D,pmvnorm.algorithm=GenzBretz())
