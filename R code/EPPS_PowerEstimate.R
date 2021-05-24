######################################
#
# Analysis of Lung Health Study Data
#
######################################
#
# Feb 25th, 2021
#
rm(list = ls())
T1 = proc.time()

# Getting parameter values   
args <- commandArgs(trailingOnly = TRUE)
myprop = as.numeric(args[1])    # the split proportion for screening and test data, \pi,    
myalpha = as.numeric(args[2])   # the threshold for screening, \alpha_{1} 

# Load R Libraries
library(MASS)
library(foreach)
library(doParallel)

# Detect Cores in the Current Node 
numCores <- detectCores()
registerDoParallel(numCores) 
# Or specify a number for your PC
#registerDoParallel(8)

path = "./ubclxing/EPPS_April2021/"
setwd(path)
source("PowerFunction.R")
mydata <- get(load("Data/LHS_5000SNP_DataSS.rda"))
dim(mydata)
#[1] 1774 5004

# Winsorize for the change of FEV1 due to its long tails
cut_point_top <- quantile(mydata$fev1.change, 0.98, na.rm=T)
cut_point_bottom <- quantile(mydata$fev1.change, 0.02, na.rm=T)

mydata$fev1.change[mydata$fev1.change > cut_point_top] = cut_point_top
mydata$fev1.change[mydata$fev1.change < cut_point_bottom] = cut_point_bottom


# Run Power Estimate function
set.seed(2021)
myresult <- mypower(iprop = myprop, ialpha = myalpha, irepli = 100, idata = mydata)

# Save Result 
save(myresult, file = paste("Data/Myprop_", myprop,
                            "_Alpha_", myalpha, ".rda", sep = ""))
T2 = proc.time()
print(T2-T1)





