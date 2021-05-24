###################################
#
# One Replicate For EPPS Method
#
##################################
#
# May, 2021
#
library(dplyr)

OneReplicate <- function(iprop = myprop, ialpha = myalpha, irepli = 5, idata = mydata){
# iprop: \pi the split proportion; 
# ialpha: \alpha_{1} the screening threshold;
# irepli: the number of replication in EPPS; 
# idata: part of real data subsetted with the top 5000 SNPs. 
 

  # SNP index is from 14th column until the end   
  snpindex <- 14:dim(idata)[2]
  snpnames <- colnames(idata)[snpindex]
  
  # Create a data frame to contain results 
  myoutput <- data.frame(snpid = integer(), p.val = double())
  
  
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
  
  # only keep SNPs passed screening $nn$ times (5% of the number of replication)
  nn <- irepli*0.05
  mylist <- table(unlist(myoutput[, 1]))
  SNPids <- as.numeric(names(mylist)[mylist >= nn ])
  #SNPids <- unique(unlist(myoutput[, 1]))
  M = length(SNPids)
  SNPnames <- colnames(idata)[SNPids]
  
  # no SNPs passed screening, return emperty data frame. 
  if(M == 0 ) {myoutput <- data.frame(SNPnames=character(),
                                      myraw.p=double(),
                                      myadjust.p.fdr=double(),
                                      myadjust.p.bf=double(),
                                      myadjust.p.by=double())
  return(myoutput)
  }
  
  # M SNPs passed the screening, run the ensemble test  
  if(M != 0 ) {
    
    icoef <- matrix(NA, nrow = M, ncol = irepli)
    ise <- matrix(NA, nrow = M, ncol = irepli)
    
    # Obtain estimate and its SE for each test  
    for (jj in 1:irepli){
      
      myresult.list <- foreach (ii = SNPids) %dopar% {

        idata.sub <- idata[, c(1:13, ii)]
        set.seed(jj)
        
        mytrain <- idata.sub %>% 
          mutate(quartile = ntile(fev1.change, 4)) %>%
          group_by(!!as.symbol(snpnames[ii-13]), quartile)%>%
          slice_sample(prop = iprop)
        
        mytest <- anti_join(idata.sub, mytrain, by = 'patient')
        mytemp2 <- lm(mytest$fev1.change ~ mytest[, 14] + mytest$sex + mytest$age + mytest$group +
                        mytest$smoke + mytest$bmi + mytest$fev.post.s2 + mytest$PC1 + mytest$PC2 + 
                        mytest$PC3 + mytest$PC4 + mytest$PC5 + mytest[, 14]*mytest$group)
        c(coef(summary(mytemp2))[2, c(1, 2)])
      }
      
      myresult <- data.frame(matrix(unlist(myresult.list), byrow = T, ncol = 2))
      icoef[, jj] <-  myresult[, 1] 
      ise[, jj ] <-  myresult[, 2] 
    }
    
    
    # Calculate the ensembel test based on mathematical formula provided in the paper
    ivec <- matrix(1, nrow = irepli, ncol = 1)
    iBeta <- icoef %*% ivec
    iSE <- rep(NA, M)
    
    for (kk in 1:M){
      temp.vec <- c(ise[kk, ])
      iexd <- expand.grid(temp.vec, temp.vec)
      ivar.multi <- matrix(iexd[, 1]*iexd[, 2], nrow = irepli)
      iSE[kk] <- sqrt(2*sum((1-iprop)*ivar.multi[upper.tri(ivar.multi)]) + sum(diag(ivar.multi)))}
    
    # P value of the ensemble test 
    itest <- abs(iBeta/iSE)
    myraw.p <- 2*(1-pnorm(itest))
    myadjust.p.fdr <- p.adjust(myraw.p, method = "fdr")
    myadjust.p.bf <- p.adjust(myraw.p, method = "bonferroni")
    myadjust.p.by <- p.adjust(myraw.p, method = "BY")
    
    # Output with identified important SNPs and their p-values
    myoutput <- data.frame(SNPnames, myraw.p, myadjust.p.fdr, myadjust.p.bf, myadjust.p.by)
    colnames(myoutput) <- c("snp", "rawp", "fdr", "bf", "by")
    return(myoutput)
  }
}


# OneReplicate(iprop = 0.5, ialpha = 0.01, irepli = 10, idata = mydata)
