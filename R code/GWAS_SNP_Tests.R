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
path = "./ubclxing/EPPS_April2021/"
setwd(path)

mydata <- get(load("Data/LHS_SNP_Data.rda"))

# Winsorize the outcome variable due to its long tails
cut_point_top <- quantile(mydata$fev1.change, 0.98, na.rm=T)
cut_point_bottom <- quantile(mydata$fev1.change, 0.02, na.rm=T)
mydata$fev1.change[mydata$fev1.change > cut_point_top] = cut_point_top
mydata$fev1.change[mydata$fev1.change < cut_point_bottom] = cut_point_bottom


n.len = ncol(mydata)

interact.effect = snp.effect = age.effect = sex.effect = intercept.effect = 
group.effect = smoke.effect = baselinefev1.effect <- data.frame(matrix(NA, nrow = n.len-13, ncol = 4))
error.sd <- rep(NA, n.len-13)
Rsquared <- rep(NA, n.len-13)

for (ii in 14:n.len){

   mysub <- mydata[, c(1:13, ii)]
   mylm = lm(formula = mysub$fev1.change ~  mysub[, 14] + mysub$age + mysub$sex + mysub$group + 
                mysub$smoke + mysub$bmi + mysub$fev.post.s2 + mysub$PC1 + mysub$PC2 + 
                mysub$PC3 + mysub$PC4 + mysub$PC5 + mysub$group*mysub[, 14], na.action = na.exclude)
   snp.effect[(ii-13), 1] = age.effect[(ii-13), 1] = sex.effect[(ii-13), 1] = 
      intercept.effect[(ii-13), 1] = group.effect[(ii-13), 1] = 
      interact.effect[(ii-13), 1] = smoke.effect[(ii-13), 1] = baselinefev1.effect[(ii-13), 1] <- colnames(mysub)[14]
   temp <- coef(summary(mylm))
   intercept.effect[(ii-13), c(2:4)] <- temp[1, c(1, 2, 4)]
   snp.effect[(ii-13), c(2:4)] <- temp[2, c(1, 2, 4)]
   age.effect[(ii-13), c(2:4)] <- temp[3, c(1, 2, 4)]
   sex.effect[(ii-13), c(2:4)] <- temp[4, c(1, 2, 4)]
   group.effect[(ii-13), c(2:4)] <- temp[5, c(1, 2, 4)]
   smoke.effect[(ii-13), c(2:4)] <- temp[6, c(1, 2, 4)]
   baselinefev1.effect[(ii-13), c(2:4)] <- temp[8, c(1, 2, 4)]
   interact.effect[(ii-13), c(2:4)]<- temp[14, c(1, 2, 4)]
   
   k <- length(mylm$coefficients)-1 #Subtract one to ignore intercept
   SSE <- sum(mylm$residuals**2)
   n <- length(mylm$residuals)
   error.sd[ii-13] <- sqrt(SSE/(n-(1+k))) #Residual Standard Error
   Rsquared[ii-13] <- summary(mylm)$r.squared
}

save(Rsquared, file = paste(path, "Data/Rsquared.rda", sep = ""))
save(interact.effect, file = paste(path, "Data/intact.rda", sep = ""))
save(snp.effect, file = paste(path, "Data/snp.rda", sep = ""))
save(error.sd, file = paste(path, "Data/errorsd.rda", sep = ""))
save(intercept.effect, file = paste(path, "Data/int.rda", sep = ""))
save(age.effect, file = paste(path, "Data/age.rda", sep = ""))
save(sex.effect, file = paste(path, "Data/sex.rda", sep = ""))
save(group.effect, file = paste(path, "Data/group.rda", sep = ""))
save(smoke.effect, file = paste(path, "Data/sex.rda", sep = ""))
save(baselinefev1.effect, file = paste(path, "Data/group.rda", sep = ""))
