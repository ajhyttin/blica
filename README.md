# Blica
Code for the Blica method published at UAI2022:

Binary Independent Component Analysis: A Non-stationarity-based Approach. Conference on Uncertainty in Artificial Intelligence 2022.
https://arxiv.org/abs/2111.15431

Note that the iVAE-based method described in the paper can be found in a 
separate repo: 
https://github.com/vitoriapacela/iVAE

Details:
--------

This the R code for the method by the original authors. This code is research code, and thus aim is to provide an implementation to 
method in the paper. While much of the method is likely easy to reproduce from the paper, i.e. in Python, the code provides:

-How to calculate the derivatives of the scaled Gaussian likelihood. This is needed for LBFGS to work. They are present in scaledgaussian.R

-An algorithm to calculate the correlations from binary data. This is in blica_correlation.R.

-Data creating from the Blica model at infinite sample limit in blica_createdata.R. Only distributions for each pair are created as the Blica
method uses only those.

Further improvements planned: 
-----------------------------

-Weighing different segments by their sample size (instead on weighing each segment equally 
regardless of sample count).

-Loading sample data in a format that can be used by Blica.

-Reproducing plots similar in the paper.

How to use:
-----------

To start go to the directory of the code and start R. Type the following to load the code:

> source('load.R')
> loud()

Run tests in test.R:

> test1()

Install any packages that are reported missing.

How to use Blica in your use case should be appararent from the code of the tests in test.R.

Contact:
--------

If you cannot solve any problems yourself 

Antti Hyttinen
ajhyttin@gmail.com



