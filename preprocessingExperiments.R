require(signal)
require(robfilter)
library(xts)

set.seed(0)

silencingPeriod <- 15

# define PB sandbox params
pbParams <- list()
pbParams$filterCoeff_b <- c(0.1652,0.2617,0.1652) #moving average coeff
pbParams$filterCoeff_a <- c(1.0000,-0.6671,0.5034) #autoregressive coeff
pbParams$plots <- 'sortAndFit_img' # {sortAndFit_img, sortAndFit, rgb_hsv, filter, ts} comparative plots to show
pbParams$deriv <- 1 # {0,1,2} = {x, d(x), d^2(x)}
pbParams$analyte <- 'K.pneumoniae'
pbParams$trialIdx <- 15

# define spots to be excluded
excludedSpots <-  c(1, 10, 12, 19, 62, 71, 80)
# define RGB indices
excludedIndicatorsB 	<- 	excludedSpots*3
excludedIndicatorsG 	<-  (excludedIndicatorsB -1)
excludedIndicatorsR 	<-  (excludedIndicatorsG -1)

# assemble parameter list for manual road testing - these will be passed on automatically
parameters       			        <- 	list()
parameters$centralRepoRoot 		<-	'/Users/michaelwalton/workspace/growthDetectionExperiments'
parameters$selectedStudy 		  <-	'2015-03--Sepsis-TF'
parameters$spotLightVersion   <- 	'03-30-2015_15-14-06'

# go to root
setwd(parameters$centralRepoRoot)

# go to study
setwd(parameters$selectedStudy)

if(parameters$selectedStudy=='2014-04--Sepsis-TF'){
  parameters$inputFile           <- '2014-12-23_11-27-15_2014-04--Sepsis-TF_YC.csv'
  parameters$spotLightVersion  	 <- 	'12-17-2014_12-27-56'
} else if (parameters$selectedStudy=='2015-03--Sepsis-TF'){
  parameters$inputFile           <- '2015-05-21_09-35-46_2015-03--Sepsis-TF_MW.csv'
  parameters$spotLightVersion    <- 	'03-30-2015_15-14-06'
}

metadata 	<- read.table(paste(parameters$inputFile, sep='/'), sep=',', header=TRUE, quote='"', comment.char="")
metadata 	<- metadata[which(metadata$ArrayType == "SLII-156"),]
metadata 	<- metadata[-which(metadata$Ignore == TRUE),]

folderStrings 	<- metadata$Folder
folderStrings 	<- gsub("\\\\","/",folderStrings)
folderStrings 	<- gsub("^/", "", folderStrings)
metadata$Folder <- folderStrings
rm(folderStrings)

doHsv <- function(csaData, maxColorValue = 4096) {
  hsvMatrix <- matrix(0, nrow=(dim(csaData)[1]), ncol=dim(csaData)[2])
  for(t in seq(dim(csaData)[1])) {
    # reshape to 3 x 80, where rows of rgb are r,g,b components
    rgb <- matrix(csaDataCropped[t,], nrow=3)
    
    # map RGB ~> HSV color space
    hsv <- rgb2hsv(rgb, maxColorValue = maxColorValue)
    hsvMatrix[t,] <- matrix(hsv, ncol = dim(hsvMatrix)[2])
  }
  return(hsvMatrix)
}

doFilter <- function(csaData, a, b){
  filtered <- matrix(0, nrow=(dim(csaData)[1]), ncol=dim(csaData)[2])
  for(c in seq(dim(csaData)[2])){
    filtered[,c] <- filter(filt=b, a=a, x=csaData[,c])
  }
  return(filtered)
}

doDeltaI_I0 <- function(csaData) {
  refImage	<- csaData[1,]
  dIoIMatrix  <- t(csaData)-refImage
  dIoIMatrix  <- t(dIoIMatrix/refImage)
  return(dIoIMatrix)
}

doDeriv <- function(csaData, order = 1) {
  if(order == 0) {
    return(csaData)
  }
  
  csaData <- doDeriv(diff(csaData), order - 1)
}

selectedTrial <- metadata[which(metadata$Analyte == pbParams$analyte),][pbParams$trialIdx,]
raw <- read.csv(file=paste(selectedTrial$Folder, parameters$spotLightVersion, 'spots.txt', sep='/'),sep="\t",header=1)
rawMat <- as.matrix(sapply(raw, as.numeric))
# make sure additional NA columns are stripped
csaData <- rawMat[which(is.na(apply(rawMat,1, sum))==FALSE),]
csaDataCropped <- rawMat[,-1]
csaDataCropped <- csaDataCropped[,-c(excludedIndicatorsR,excludedIndicatorsG,excludedIndicatorsB)]

rgbData <- csaDataCropped
hsvData <- doHsv(csaDataCropped)

# dI/I0 >> slope >> filter
dsfRgb <- doFilter(doDeriv(doDeltaI_I0(rgbData), pbParams$deriv), pbParams$filterCoeff_a,pbParams$filterCoeff_b)
dsfHsv <- doFilter(doDeriv(doDeltaI_I0(hsvData), pbParams$deriv), pbParams$filterCoeff_a,pbParams$filterCoeff_b)
# raw
rawRgb <- doDeriv(doDeltaI_I0(rgbData), pbParams$deriv)
rawHsv <- doDeriv(doDeltaI_I0(hsvData), pbParams$deriv)

if (pbParams$plots == 'filter') {
  end <- dim(dsfRgb)[1]
  
  for(c in seq(dim(dsfRgb)[2])){
    plot(raw[silencingPeriod:end,c], type='l',col='red')
    lines(dsfRgb[silencingPeriod:end,c], type='l',col='blue')
    legend("topright", legend=c("Raw", "Filtered"), cex=0.75, lwd=2, col=c("red", "blue"))
  }
} else if (pbParams$plots == 'rgb_hsv') {
  end <- dim(dsfHsv)[1]
  
  for(i in seq(dim(dsfHsv)[2] / 3)) {
    idx <- (i - 1) * 3
    
    hsvlims <- c(min(dsfHsv[,idx + c(1,2,3)]), max(dsfHsv[,idx + c(1,2,3)]))
    rgblims <- c(min(dsfRgb[,idx + c(1,2,3)]), max(dsfRgb[,idx + c(1,2,3)]))
    
    # Multi-plot; add extra space to right of plot area; change clipping to figure
    attach(mtcars)
    par(mfrow=c(2,1))
    
    # PLOT HSV
    plot(dsfHsv[,idx + 1], type='l',col='darkorange', ylim = hsvlims, main = paste('Indicator #', i))
    lines(dsfHsv[,idx + 2], type='l',col='purple')
    lines(dsfHsv[,idx + 3], type='l',col='black')
    legend("bottomright", legend=c("Hue", "Saturation", "Value"), cex=0.5, lwd=2, col=c("darkorange", "purple","black"))
    
    # PLOT RGB
    plot(dsfRgb[,idx + 1], type='l',col='red', ylim = rgblims)
    lines(dsfRgb[,idx + 2], type='l',col='darkgreen')
    lines(dsfRgb[,idx + 3], type='l',col='blue')
    legend("bottomright", legend=c("Red", "Green", "Blue"), cex=0.5, lwd=2, col=c("red", "darkgreen","blue"))
  }
} else if (pbParams$plots == 'sortAndFit' || pbParams$plots == 'sortAndFit_img') {
  greyScale <- grey(seq(0, 1, length = 256))
  indicatorIndex <- seq(dim(dsfRgb)[2])
  regressionSlope <- matrix(nrow = dim(dsfRgb)[1], ncol = 2)
  
  for(t in seq(dim(dsfRgb)[1])) {
    rgbSorted <- dsfRgb[,order(dsfRgb[t,])]
    hsvSorted <- dsfHsv[,order(dsfHsv[t,])]
    
    rgbFit <- lm(rgbSorted[t,] ~ indicatorIndex)
    hsvFit <- lm(hsvSorted[t,] ~ indicatorIndex)
    
    regressionSlope[t,1] <- rgbFit$coefficients[1]
    regressionSlope[t,2] <- hsvFit$coefficients[1]
    
    if (pbParams$plots == 'sortAndFit_img') {
      img <- as.matrix(rgbSorted[silencingPeriod:t,])
      image(img, main = paste('Indicators sorted by t = ', t), col = greyScale)
    }
  }
  
  plot(regressionSlope[silencingPeriod:dim(regressionSlope)[1],1], type='l', col = 'red')
  #plot(regressionSlope[,2], type='l', col = 'blue')
  
} else if (pbParams$plots == 'ts') {
  ts.plot(dsfHsv, col=rainbow(dim(dsfHsv)[2]))
  ts.plot(dsfRgb, col=rainbow(dim(dsfRgb)[2]))
} else if (pbParams$plots == 'histo') {
  
  inputXLim <- c(min(dsfRgb), max(dsfRgb))
  inputYLim <- c(0,dim(dsfRgb)[2])
  
  for(t in seq(dim(dsfRgb)[1])) {
    hist(dsfRgb[t,], breaks = 10, main = "RGB Indicator Distribution", 
         xlab = paste(pbParams$mode, '[', t, ']'), ylab = 'n indicators', xlim = inputXLim, ylim = inputYLim)
  }
}
