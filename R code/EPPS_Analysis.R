#######################################
#
# Using EPPS to Analyze the LHS Data
#
######################################
#
# May 23, 2020
#
rm(list = ls())

# Load R Libraries
library(MASS)
library(foreach)
library(doParallel)

# Set up the number of cores used in parallel computing 
numCores <- detectCores()
registerDoParallel(numCores)

# Set up the path and source the function for one replication
path = "/ubclxing/EPPS_April2021/"
setwd(path)
source("OneReplicate.R")

# Load Data and Winsorize the outcome variable due to its long tails 
mydata <- get(load("Data/LHS_5000SNP_DataAnalysis.rda"))
dim(mydata)
#[1] 1774 5004

cut_point_top <- quantile(mydata$fev1.change, 0.98, na.rm=T)
cut_point_bottom <- quantile(mydata$fev1.change, 0.02, na.rm=T)
mydata$fev1.change[mydata$fev1.change > cut_point_top] = cut_point_top
mydata$fev1.change[mydata$fev1.change < cut_point_bottom] = cut_point_bottom

# Set up parameter values, which were obtained based on optimizing the estimated power.  
myprop = 0.5   # split proportion, \pi
myalpha = 0.00001  # screening threshold, \alpha_{1}

# Run EPPS 
set.seed(2021)
myresult <- OneReplicate(iprop = myprop, ialpha = myalpha, irepli = 100, idata = mydata)

# Save Results
save(myresult, file = paste("Data/myanalysis_result.rda", sep = ""))
