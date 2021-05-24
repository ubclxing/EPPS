########################################################
#
# Power Function to Estimate Power in EPPS 
#
########################################################
#
# April, 2021
# 
library(dplyr)

mypower <- function(iprop = myprop, ialpha = myalpha, irepli = 5, idata = mydata){
  
  # SNPs are from 14th column to the end 
  snpindex <- 14:dim(idata)[2]
  snpnames <- colnames(idata)[snpindex]
  
  # Create a data frame to contain results
  myoutput <- data.frame(snpid = integer(), p.val = double())

  # Screening Process   
  for (jj in 1:irepli){

    myresult.list <- foreach (ii = snpindex) %dopar% {
      
      idata.sub <- idata[, c(1:13, ii)]
      set.seed(jj)
      
      # Stratified sampling based on quartiles of change in FEV1 and SNP levels 
      mytrain1 <- idata.sub %>% 
       mutate(quartile = ntile(fev1.change, 4)) %>%
       group_by(!!as.symbol(snpnames[ii-13]), quartile)%>%
       slice_sample(prop = iprop)
      
      # Cells with less than 4 observations are removed from screening to avoid outliers 
      mytrain <- mytrain1 %>% 
        group_by(!!as.symbol(snpnames[ii-13]), quartile)%>%
        mutate(number = n())%>% 
        filter(number > 4)
      
      # Screening model 
      mytemp <- lm(mytrain$fev1.change ~ mytrain[, 14][[1]] +  mytrain$group + mytrain$sex + mytrain$age +
                     mytrain$smoke + mytrain$bmi + mytrain$fev.post.s2 + mytrain$PC1 + mytrain$PC2 + 
                     mytrain$PC3 + mytrain$PC4 + mytrain$PC5 + mytrain[, 14][[1]]*mytrain$group)
      
      # P-value for the SNP main effect 
      pp <- summary(mytemp)$coef[2, 4]
      c(ii, pp)}
    
    
    myresult <- data.frame(matrix(unlist(myresult.list), byrow = T, ncol = 2))
    ind <- (myresult[, 2] < ialpha)
    if (sum(ind) > 0) { mytt <- myresult[ind, ]
    myoutput <- rbind(mytt, myoutput)}
  }
  
  # SNPs passed the screen once 
  SNPids <- unique(unlist(myoutput[, 1]))
  M = length(SNPids)
  
  if(M == 0 ) {return(0)}    # no SNPs passed the screening, 0 is return
  
  # M SNPs passed screening, calculate the estimated Power based the math formula in the paper    
  if(M != 0 ) {
    icoef <- matrix(NA, nrow = M, ncol = irepli)
    ise <- matrix(NA, nrow = M, ncol = irepli)
    
    # Obtain estimate and its SE for each screening data for the selected SNPs 
    for (jj in 1:irepli){
      
      myresult.list <- foreach (ii = SNPids) %dopar% {
        
        idata.sub <- idata[, c(1:13, ii)]
        set.seed(jj)
        mytrain1 <- idata.sub %>% 
          mutate(quartile = ntile(fev1.change, 4)) %>%
          group_by(!!as.symbol(snpnames[ii-13]), quartile)%>%
          slice_sample(prop = iprop)
        
        mytrain <- mytrain1 %>% 
          group_by(!!as.symbol(snpnames[ii-13]), quartile)%>%
          mutate(number = n())%>% 
          filter(number > 4)
        
        mytemp2 <- lm(mytrain$fev1.change ~ mytrain[, 14][[1]] +  mytrain$group + mytrain$sex + mytrain$age +
                        mytrain$smoke + mytrain$bmi + mytrain$fev.post.s2 + mytrain$PC1 + mytrain$PC2 + 
                        mytrain$PC3 + mytrain$PC4 + mytrain$PC5 + mytrain[, 14][[1]]*mytrain$group)
        c(coef(summary(mytemp2))[2, c(1, 2)])
      }
      
      myresult <- data.frame(matrix(unlist(myresult.list), byrow = T, ncol = 2))
      icoef[, jj] <-  myresult[, 1] 
      ise[, jj ] <-  myresult[, 2] 
    }
    
    ivec <- matrix(1, nrow = irepli, ncol = 1)
    iBeta <- icoef %*% ivec
    iSE <- rep(NA, M)
    ilen <- dim(idata)[1]
    iodd = (iprop*ilen - 4)/((1-iprop)*ilen-4)
    
    for (kk in 1:M){
      
      temp.vec <- c(ise[kk, ])
      iexd <- expand.grid(temp.vec, temp.vec)
      ivar.multi <- matrix(iexd[, 1]*iexd[, 2], nrow = irepli)
      iSE[kk] <- sqrt(iodd*(2*sum((1-iprop)*ivar.multi[upper.tri(ivar.multi)]) + sum(diag(ivar.multi))))}
    
    # Power estimation
    itest <- iBeta/iSE
    ialpha2 = 0.05
    Z1 <- qnorm(1-ialpha2/(2*M))
    Z2 <- qnorm(ialpha2/(2*M))
    
    # the estimated power
    myreturn <- c(1 - pnorm(Z1-itest) + pnorm(Z2-itest))
    return(myreturn)
  }
}

