#####
# Installing packages
install.packages(c('raster','readr','tidyr','dplyr','RColorBrewer','GISTools','sp',
                   'sf','ggplot2','bioimagetools','plyr','orcutt','lmtest','readxl'), dependencies = T)

library(raster)
library(readr)
library(tidyr)
library(plyr)
library(dplyr)
library(RColorBrewer)
library(GISTools)
library(sp)
library(sf)
library(ggplot2)
library(bioimagetools)
library(orcutt)
library(lmtest)
library(readxl)

#####
# Define stacks in my folder
rainanom_ <- paste0('~/R/PWR/tls_rainanom/rainanom_',0:443,".tif")
rainanom <- stack(rainanom_)

sstanom_ <- paste0('~/R/PWR/tls_sstanom/sstanom_',0:443,".tif")
sstanom <- stack(sstanom_)

s <- stack(rainanom, sstanom)

a <- length(rainanom_)
b <- length(sstanom_)

#####
# Slope Function
## Linear Model for Time-Series Data
fun=function(x) { if (is.na(x[1])){ NA } else { cochrane.orcutt(lm(x[1:a] ~ x[(1+a):(a+b)]))$coefficients[2] }}
slope <- calc(s, fun)

# Saving Slope into geotiff
writeRaster(slope, filename = 'tls_slope_rainanom.tif')

# Open the geotiff file slope
slope <- raster('./tls_slope_rainanom.tif')

# Slope Visualization
par(mar=c(1,1,1,1))
plot(slope, col = c(rgb(112/255,23/255,29/255),
                    rgb(213/255,72/255,47/255),
                    rgb(237/255,145/255,79/255),
                    rgb(248/255,203/255,111/255),
                    rgb(255/255,253/255,187/255),
                    rgb(204/255,204/255,204/255)),
     asp = 1,
     box = F,
     ylab = NA,
     xlab = NA,
     axes = F,
     legend = F)
legend(x = 'topleft',
       legend = c('-10 to 0','-20 to -10','-30 to -20','-40 to -30','-50 to -40','< -50'),
       fill = c(rgb(204/255,204/255,204/255),
                rgb(255/255,253/255,187/255),
                rgb(248/255,203/255,111/255),
                rgb(237/255,145/255,79/255),
                rgb(213/255,72/255,47/255),
                rgb(112/255,23/255,29/255)), 
       bty = 'n')

#####
# P-Value of the Slope

fun1=function(x) { if (is.na(x[1])){ NA } else { cochrane.orcutt(lm(x[1:a] ~ x[(1+a):(a+b)]))$DW[4] }}
pval <- calc(s, fun1)

plot(pval)

#####
# Gridcorts Function

gridcorts <- function(rasterstack, method, type=c("corel","pval","both")){
  # Values for (layers, ncell, ncol, nrow, method, crs, extent) come straight from the input raster stack
  # e.g. nlayers(rasterstack), ncell(rasterstack)... etc.
  print(paste("Start Gridcorts:",Sys.time()))
  print("Loading parameters")
  layers=nlayers(rasterstack);ncell=ncell(rasterstack);
  ncol=ncol(rasterstack);nrow=nrow(rasterstack);crs=crs(rasterstack);
  extent=extent(rasterstack);pb = txtProgressBar(min = 0, max = ncell, initial = 0)
  print("Done loading parameters")
  mtrx <- as.matrix(rasterstack,ncol=layers)
  empt <- matrix(nrow=ncell, ncol=2)
  print("Initiating loop operation")
  if (type == "corel"){
    for (i in 1:ncell){
      setTxtProgressBar(pb,i)
      if (all(is.na(mtrx[i,1:(layers/2)])) | all(is.na(mtrx[i,((layers/2)+1):layers]))){ 
        empt[i,1] <- NA 
      } else 
        if (sum(!is.na(mtrx[i,1:(layers/2)]/mtrx[i,((layers/2)+1):layers])) < 4 ){
          empt[i,1] <- NA 
        } else 
          empt[i,1] <- as.numeric(cor.test(mtrx[i,1:(layers/2)], mtrx[i,((layers/2)+1):layers],method=method)$estimate)
    }
    print("Creating empty raster")
    corel <- raster(nrows=nrow,ncols=ncol,crs=crs)
    extent(corel) <- extent
    print("Populating correlation raster")
    values(corel) <- empt[,1]
    print(paste("Ending Gridcorts on",Sys.time()))
    corel
  } 
  else
    if (type == "pval"){
      for (i in 1:ncell){
        setTxtProgressBar(pb,i)
        if (all(is.na(mtrx[i,1:(layers/2)])) | all(is.na(mtrx[i,((layers/2)+1):layers]))){ 
          empt[i,2] <- NA 
        } else 
          if (sum(!is.na(mtrx[i,1:(layers/2)]/mtrx[i,((layers/2)+1):layers])) < 4 ){
            empt[i,2] <- NA 
          } else 
            empt[i,2] <- as.numeric(cor.test(mtrx[i,1:(layers/2)], mtrx[i,((layers/2)+1):layers],method=method)$p.value)
      }
      pval <- raster(nrows=nrow,ncols=ncol,crs=crs)
      extent(pval) <- extent
      print("Populating significance raster")
      values(pval) <- empt[,2]
      print(paste("Ending Gridcorts on",Sys.time()))
      pval
    }
  else
    if (type == "both"){
      for (i in 1:ncell){
        setTxtProgressBar(pb,i)
        if (all(is.na(mtrx[i,1:(layers/2)])) | all(is.na(mtrx[i,((layers/2)+1):layers]))){ 
          empt[i,] <- NA 
        } else 
          if (sum(!is.na(mtrx[i,1:(layers/2)]/mtrx[i,((layers/2)+1):layers])) < 4 ){
            empt[i,] <- NA 
          } else {
            empt[i,1] <- as.numeric(cor.test(mtrx[i,1:(layers/2)], mtrx[i,((layers/2)+1):layers],method=method)$estimate) 
            empt[i,2] <- as.numeric(cor.test(mtrx[i,1:(layers/2)], mtrx[i,((layers/2)+1):layers],method=method)$p.value)
          }
      }
      c <- raster(nrows=nrow,ncols=ncol,crs=crs)
      p <- raster(nrows=nrow,ncols=ncol,crs=crs)
      print("Populating raster brick")
      values(c) <- empt[,1]
      values(p) <- empt[,2]
      brk <- brick(c,p)
      extent(brk) <- extent
      names(brk) <- c("Correlation","Pvalue")
      print(paste("Ending Gridcorts on",Sys.time()))
      brk
    }
}

#####
# Correlation

correlation <- gridcorts(rasterstack = s, method = "pearson", type = "corel")

# Correlation Plot
plot(correlation)

# Setting directory
getwd()
setwd('~/R/PWR')

# Saving Correlation into geotiff
writeRaster(correlation, filename = 'tls_correlation.tif')

# Open the correlation geotiff file
k1 <- raster('~/R/PWR/tls_correlation.tif')
par(mar=c(1,1,1,1))
plot(k1)

# Correlation Plot Visualization
par(mar=c(1,1,1,1))
plot(-correlation,breaks = c(0.18,0.23,0.26,0.3,0.35), 
     col = brewer.pal(5,'Reds'),
     asp = 1,
     box = F,
     ylab = NA,
     xlab = NA,
     axes = F,
     legend = F)
legend(x = 'topleft',
       legend = c('< 0.2','0.2 to 0.23','0.23 to 0.26','0.26 to 0.3','> 0.3'),
       fill = brewer.pal(5,'Reds'), 
       bty = 'n')

