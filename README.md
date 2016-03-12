# PQSQ-DataApproximators

This is a set of Matlab procedures for performing PQSQ-based data approximation.

In particular,

PQSQmean - computes the mean value with PQSQ approximation error 

pcaPQSQ - computes PCA with PQSQ approximation error 

pcaL1 - computes approximative l1-based PCA

For large data sets, we include wrappers of the PQSQ-based procedures implemented in Java, providing several fold gain in the computational time.

Simplest examples of use:

>> x = load('test.txt'); % rows of x are objects (data points), and columns are variables (data space coordinates)

>>[V,U,C] = pcaL1(x,2); plot(U(:,1),U(:,2),'ko'); axis equal;
% computes 2 first approximative l1-based principal components and plots the distribution of points

>>[V,U,C] = pcaPQSQ(x,2,@L2); 
% computes 2 first l2-based error principal components

>>[V,U,C] = pcaPQSQ(x,3,@L1,'intervals',defineIntervals(x,10)); 
% computes 3 first l1-based error principal components, with increased accuracy of l1 approximation (10 intervals instead of default 5)
