#####################################
#
# SNP-wise tests for LHS data
#
#####################################
#
# April, 2021
# Li Xing
#
#
rm(list = ls())

path = "/ubclxing/EPPS_April2021/"
setwd(path)

# LHS_SNP_Data is the original data, which is not attached here.  
mydata <- get(load("Data/LHS_SNP_Data.rda"))

# Winsorize the outcome variable
cut_point_top <- quantile(mydata$fev1.change, 0.98, na.rm=T)
cut_point_bottom <- quantile(mydata$fev1.change, 0.02, na.rm=T)
mydata$fev1.change[mydata$fev1.change > cut_point_top] = cut_point_top
mydata$fev1.change[mydata$fev1.change < cut_point_bottom] = cut_point_bottom

n.len = ncol(mydata)
snp.effect <- data.frame(matrix(NA, nrow = n.len-13, ncol = 4))

for (ii in 14:n.len){
#  ii: column number of a SNP
   
   mysub <- mydata[, c(1:13, ii)]
   mylm = lm(formula = mysub$fev1.change ~  mysub[, 14] + mysub$age + mysub$sex + mysub$group + 
                mysub$smoke + mysub$bmi + mysub$fev.post.s2 + mysub$PC1 + mysub$PC2 + 
                mysub$PC3 + mysub$PC4 + mysub$PC5 + mysub$group*mysub[, 14], na.action = na.exclude)
   snp.effect[(ii-13), 1] <- colnames(mysub)[14]
   temp <- coef(summary(mylm))
   snp.effect[(ii-13), c(2:4)] <- temp[2, c(1, 2, 4)]
}

save(snp.effect, file = paste(path, "Data/snp.rda", sep = ""))


