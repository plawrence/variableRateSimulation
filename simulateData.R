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
#f <- GaussRF(x=x, y=y, model=model,grid=T, param=c(mean, variance, nugget, scale, alpha))   

#simulate the exponential model for both precip and soil texture 
soilN <- GaussRF(x=x, y=y, model=model,grid=T, param=c(mean,variance,scale,range))
set.seed(333)
variance <- 20
scale <- .4
range <- 50
soilTex <- GaussRF(x=x, y=y, model=model, grid=T, param=c(mean,variance,scale,range))

##Set up curves for the different bin areas
##Grabbed these coefficients off of the extension page http://landresources.montana.edu/fertilizerfacts/17_Predicting_Spring_Wheat_Yield.htm
bin1beta0 = 21
bin2beta0 = 24.4
bin3beta0 = 30.2
bin1beta1 = .18
bin2beta1 = .36
bin3beta1 = .31
bin1beta2 = -.00051
bin2beta2 = -.00113
bin3beta3 = -.0003
bincoefmat = matrix(0,nrow=3,ncol=3)
bincoefmat[,1] = c(bin1beta0,bin1beta1,bin1beta2)
bincoefmat[,2] = c(bin2beta0,bin2beta1,bin2beta2)
bincoefmat[,3] = c(bin3beta0,bin3beta1,bin3beta2)

#Make realistic values for soil N, simulate uniform application and yield response.  Sd quickly generates noise depending on
#the underlying variation - NOTE: no differential bin responses to start
soilN = soilN*2+40

##Plot the initial soil N
par(mfrow=c(1,4))
image(x,y,soilN)

#Initial precip
precipYr = rnorm(1,7,2)
precipYr
soilTexWat = precipYr*soilTex

appN = matrix(200,nrow=fieldrow,ncol=fieldcol)
yield = .2*(soilN+appN)+.1*soilTexWat + matrix(rnorm(fieldcol*fieldrow,mean=0,sd=.7),ncol=fieldcol,nrow=fieldrow)

##Calculate remaining soil N
nremcoef=.2
calcSoilN = function(soilN,yield,nremcoef){
  soilN = soilN - nremcoef*yield
  return(soilN)
}
#Don't calc - do it at the end of the page in the loop
#soilN = calcSoilN(soilN,yield,nremcoef)

#Plot it
image(x, y, yield)
image(x, y, soilN)
image(x, y, soilTex)
#contour(f) 

#Divide the yield into 3 bins that will be used for experimentation
calcBins <- function(fieldcol,fieldrow,yield){
  bins = matrix(0,ncol=fieldcol,nrow=fieldrow)
  quant.33 <- quantile(yield,.33)
  quant.66 <- quantile(yield,.66)
  for (i in 1:fieldrow){
    for (j in 1:fieldcol){
      if (yield[i,j] < quant.33){
        bins[i,j]=1
      } else if (yield[i,j]<quant.66 & yield[i,j]>quant.33){
        bins[i,j]=2
      } else {
        bins[i,j]=3
      }
    }
  }
  return(bins)
}
bins = calcBins(fieldcol,fieldrow,yield)

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
      if (((draw - 1) %in% exptlist & !((draw - 1) %in% expt.out)) & ((draw+1) %in% exptlist & !((draw+1) %in% expt.out)) & !((draw-2) %in% expt.out) & !((draw+2) %in% expt.out)){
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

#Function to apply treatment levels - make sure the number of treat areas is divisible by the number of levels
#Assuming a bin size of 3 for now
applyTreatLevels <- function(inmat,levelvec,reps){
  possiblelocs = which(inmat<200)
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

#par(mfrow=c(2,2))
#image(x,y,bins)
#image(x,y,applyExperiment(bins,3,3,4,4))
#image(x,y,applyExperiment(bins,2,3,4,4))
#image(x,y,applyExperiment(bins,1,3,4,4))

##Store the experimental N treatments
bin1treat = applyExperiment(bins,1,3,4,4)
bin2treat = applyExperiment(bins,2,3,4,4)
bin3treat = applyExperiment(bins,3,3,4,4)

applymin <- function(bin1treat,bin2treat,bin3treat){
  outmat = matrix(0,nrow=fieldrow, ncol=fieldcol)
  for (i in 1:fieldrow){
    for (j in 1:fieldcol){
      outmat[i,j] = min(bin1treat[i,j],bin2treat[i,j],bin3treat[i,j])
    }
  }
  return(outmat)
}

fertcollection = applymin(bin1treat,bin2treat,bin3treat)

cellYieldCalc <- function(row,col,bins,Nmat,soilTexWat){
  cell = c(row,col)
  cellbin = bins[cell[1],cell[2]]
  cellN = Nmat[cell[1],cell[2]]
  bincoef = bincoefmat[,cellbin]
  soilTW = soilTexWat[cell[1],cell[2]]
  outyld = bincoef[1]+bincoef[2]*cellN+bincoef[3]*(cellN^2)+.2*soilTW
  return(outyld)
}

newYieldCalc <- function(bins,fertmat,soilN,precipYr){
  outmat = matrix(0,nrow=fieldrow,ncol=fieldcol)
  Nmat = fertmat+soilN
  soilTexWat = precipYr*soilTex
  #Didn't want to loop but have to since we're dealing with a function applied to multiple matrices
  for (row in 1:nrow(outmat)){
    for (col in 1:ncol(outmat)){
      outmat[row,col] = cellYieldCalc(row,col,bins,Nmat,soilTexWat)+rnorm(1,mean=0,sd=2)
    }
  }
  #outmat=outer(1:nrow(outmat),1:ncol(outmat),FUN=function(r,c) cellYieldCalc(r,c,bins,Nmat))
  #outmat=outer(1:nrow(outmat),1:ncol(outmat),FUN=function(r,c) tmpfun(r,c,bins[r,c]))
  #outmat = apply(outmat,c(1,2),function(x) tmpfun(x)#cellYieldCalc(c(r,c),bins,Nmat))
  return(outmat)
}

##Continue for successive rounds - soil N is calculated from the last year's yield N removal
roundFun = function(yield,soilN,fieldcol,fieldrow){
  bins = calcBins(fieldcol,fieldrow,yield)
  image(x,y,soilN,main=as.character(precipYr))
  image(x,y,bins)
  
  #Calc new soil N for new round
  soilN = calcSoilN(soilN,yield,nremcoef)
  
  bin1treat = applyExperiment(bins,1,3,4,4)
  bin2treat = applyExperiment(bins,2,3,4,4)
  bin3treat = applyExperiment(bins,3,3,4,4)
  fertcollection = applymin(bin1treat,bin2treat,bin3treat)
  image(x,y,fertcollection)
  precipYr = rnorm(1,7,2)
  precipYr
  outyld = newYieldCalc(bins,fertcollection,soilN,precipYr)
  image(x,y,outyld,main=as.character(precipYr))
  #hist(outyld)
  plot(outyld~fertcollection,col=bins)
  outlist = list(outyld,soilN)
  return(outlist)
}


numRounds = 3  ##Number of experimental treatment-years to run
par(mfrow=c(3,5))
for (i in 1:numRounds){
  outlist = roundFun(yield,soilN,fieldcol,fieldrow)
  yield = outlist[[1]]
  soilN = outlist[[2]]
}

##Next up, run regressions on each of the bin areas from the experiment-year, then apply N accordingly