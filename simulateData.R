##Created on 1/17/2013 by plawrence
##This script is intended to create simulation data for assessing VRA adaptive management techniques

##Load required packages
library(geoR)
library(spatstat)

##Set the seed
set.seed(359)

PrintModelList() ## the complete list of implemented models - choose one!

##Create an empty dataset
model <- "stable"
mean <- 0
variance <- 10
nugget <- 0 #noise around the structure
scale <- 10  #structure/patchiness
alpha <- 2 ## see help("CovarianceFct") for additional
## parameters of the covariance functions
step <- 1 ## nicer, but also time consuming if step <- 0.1
x <- seq(0, 100, step)
y <- seq(0, 100, step)

f <- GaussRF(x=x, y=y, model=model,grid=T, param=c(mean, variance, nugget,
                                                   scale, alpha))   #

par(mfcol=c(1,2))
image(x, y, f)
contour(f) 
