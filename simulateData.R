##Created on 1/17/2013 by plawrence
##This script is intended to create simulation data for assessing VRA adaptive management techniques

##Load required packages
library(geoR)
library(RandomFields)

##Set the seed
set.seed(359)

PrintModelList() ## the complete list of implemented models - choose one!

##Create an empty dataset
model <- "exponenti"
mean <- 0
variance <- 30
scale <- .5  #structure/patchiness
range <- 40

#Set the number of field columns and rows
fieldcol=40
fieldrow=40

## parameters of the covariance functions
step <- 1 ## nicer, but also time consuming if step <- 0.1
x <- seq(0, fieldcol-1, step)
y <- seq(0, fieldrow-1, step)

#Additional parameters for the "stable" model
alpha <- 2 ## see help("CovarianceFct") for additional
nugget <- 0 #noise around the structure

#Setup for the "stable model" - ignore
f <- GaussRF(x=x, y=y, model=model,grid=T, param=c(mean, variance, nugget,
                                                   scale, alpha))   #

#simulate the exponential model
soilN <- GaussRF(x=x, y=y, model=model,grid=T, param=c(mean,variance,scale,range))

#Make realistic values for soil N, simulate uniform application and yield response.  Sd quickly generates noise depending on
#the underlying variation
soilN = soilN*2+40
appN = matrix(100,nrow=fieldrow,ncol=fieldcol)
yield = .5*(soilN+appN)+matrix(rnorm(fieldcol*fieldrow,mean=0,sd=2),ncol=fieldcol,nrow=fieldrow)

#Plot it
par(mfcol=c(1,2))
image(x, y, soilN)
image(x, y, yield)
contour(f) 

#Divide the soil N into 3 bins that will be used for experimentation
bin1 <- soilN[soilN<quantile(soilN,.33)]
bin2 <- soilN[soilN<quantile(soilN,.66) & soilN>quantile(soilN,.33)]
bin3 <- soilN[soilN>quantile(soilN,.66)]

bins = matrix(0,ncol=fieldcol,nrow=fieldrow)
quant.33 <- quantile(soilN,.33)
quant.66 <- quantile(soilN,.66)
for (i in 1:fieldrow){
  for (j in 1:fieldcol){
    if (soilN[i,j] < quant.33){
      bins[i,j]=1
    } else if (soilN[i,j]<quant.66 & soilN[i,j]>quant.33){
      bins[i,j]=2
    } else {
      bins[i,j]=3
    }
  }
}

image(x,y,bins)

#Find areas in each bin that are 3 cells in height, then randomly apply 4 treatments

#function to search individual columns for areas that match a particular bin and are of sufficient size
findTreats = function(colvec,bin,binsize){
  binvec = match(colvec,bin)
  matchvec = vector(length=length(binvec))
  count = 0
  
  #For each column, count the numbers of consecutive within-bin cells
  for (i in 1:length(binvec)){
    if (!is.na(binvec[i])){
      count =count+1
      matchvec[i] = count
    } else {
      count = 0
      matchvec[i] = 0
    }
    #Reset when the desired binsize has been met
    if (count==binsize){
      count = 0
    }
  }
  locationvec = which(matchvec==binsize)
  #Set the corresponding locations to the number of the bin
  treatmentvec = rep(0,length(binvec))
  treatmentvec[locationvec]=bin
  #Make sure that the binsize-vector cells are all coded
  binsize = binsize - 1
  treatmentvec[locationvec-binsize] = bin
  while (binsize >0){
    treatmentvec[locationvec-binsize] = bin
    binsize = binsize - 1
  }
  return (treatmentvec)    
}

#Function to apply the random experimental treatments to each of the bins - specify binsize, # of reps
applyExperiment = function(dataset,bin,binsize,reps,numlevels){
  exptvec = matrix(0,ncol=fieldcol,nrow=fieldrow)
  for (i in 1:fieldcol){
    colvec = as.vector(dataset[,i])
    treatvec = findTreats(colvec,bin,binsize)
    exptvec[,i] = treatvec
  }
  exptlist = which(exptvec==bin)
  expt.out = vector()
  for (i in 1:(reps*numlevels)){
    surrounded = FALSE
    while (surrounded == FALSE){
      draw = sample(exptlist,1)
      if (((draw - 1) %in% exptlist & !((draw - 1) %in% expt.out)) & ((draw+1) %in% exptlist & !((draw+1) %in% expt.out))){
        expt.out <- c(expt.out,draw)
        surrounded = TRUE
      }
    }
  }
  #Background level of 200 lbs applied
  outmat <- matrix(200,nrow=fieldrow,ncol=fieldcol)
  for (i in 1:length(expt.out)){
    outmat[expt.out[i]] = bin
  }
  levelvec <- c(0,40,80,120)
  outmat = applyTreatLevels(outmat,levelvec,reps)
  return (outmat)
}

par(mfrow=c(2,2))
image(x,y,bins)
image(x,y,applyExperiment(bins,3,3,4,4))
image(x,y,applyExperiment(bins,2,3,4,4))
image(x,y,applyExperiment(bins,1,3,4,4))

#Function to apply treatment levels - make sure the number of treat areas is divisible by the number of levels
#Assuming a bin size of 3 for now
applyTreatLevels <- function(inmat,levelvec,reps){
  possiblelocs = which(inmat>0)
  levrep <- rep(levelvec,reps)
  treatsamplevec <- sample(possiblelocs)
  count = 1
  for (i in levrep){
    treatloc = treatsamplevec[count]
    inmat[treatloc] = i
    inmat[treatloc-1] = i
    inmat[treatloc+1] = i
    count = count+ 1
    #possiblelocs = possiblelocs[-which(possiblelocs==treatloc)]
  }
  return(inmat)
}

image(x,y,applyTreatLevels(,levelvec,4))
