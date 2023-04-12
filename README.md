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

-How to calculate the derivatives of the scaled Gaussian likelihood. This is needed for LBFGS to work. They are present in scaledgaussian.R.

-An algorithm to calculate the correlations from binary data. This is in blica_correlation.R.

-Data creating from the Blica model at infinite sample limit in blica_createdata.R. Only distributions for each pair are created as the Blica
method uses only those.

Further improvements planned: 
-----------------------------

-Reproducing plots (similar to the ones) in the paper.

-Report the packages needed.

-Improving clarity, style, usability and adding comments.

-Outputting sources in a sample data form.

-Running several different seeds.

-Adding the log-likelihood function calculations to run a likelihood based method (which is inefficient in practice).

How to use:
-----------

To start go to the directory of the code and start R. Type the following to load the code:

> source('load.R')

> loud()

Install any packages that are reported missing here or during the test runs. The particular packages the code extensively uses are 'mvtnorm' for truncated normal distribution calculations and 'lpSolve' calculating the MCS taking account the order indeterminacy of the sources which is common in ICA literature.

Run tests in test.R:

> test1()

> test2()

> test3()

> test4()

> test5()

Test should finish in minutes and give MCS value close to 1. These tests run only a single seed: you might want to call blica several (e.g. 3) times and select the best run according to likelihood value obtained, to disregard runs that converge to a local optima.

How to use Blica in a use case should be apparent from the code of the tests in test.R.

Contact:
--------

If you cannot solve any problems yourself contact:

Antti Hyttinen
ajhyttin@gmail.com



