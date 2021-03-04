library(mclust)
library(MixSim)

NN <- 100 # n_1 (n_2) or n/2 sample size
DD <- 2 # Dimension of mixture model
BO <- 0.05 # Difference in mixtures for MixSim function
TrueG <- 4 # Generative number of components

# Generate random mixture model parameters based on settings
Sim_para <- MixSim(BarOmega = BO,K=TrueG,p=DD,PiLow = 0.1)

# Generate two data sets of sizes n_1 = n_2
Data1 <- simdataset(n=NN,Pi=Sim_para$Pi,Mu=Sim_para$Mu,S=Sim_para$S)
Data2 <- simdataset(n=NN,Pi=Sim_para$Pi,Mu=Sim_para$Mu,S=Sim_para$S)

# Initiate null and alternative hypotheses
GG <- 4 # number of null components
LL <- 1 # additional components in alternative

# Fit mixtures under null hypothesis
MC1_null <- densityMclust(Data1$X,G=GG,modelNames = 'VVV')
MC2_null <- densityMclust(Data2$X,G=GG,modelNames = 'VVV')

# Fit mixtures under the alternative hypothesis
MC1_alt <- densityMclust(Data1$X,G=GG+LL,modelNames = 'VVV')
MC2_alt <- densityMclust(Data2$X,G=GG+LL,modelNames = 'VVV')

# Calculate test statistics
V1 <- exp(sum(log(predict.densityMclust(MC2_alt,Data1$X)))-MC1_null$loglik)
V2 <- exp(sum(log(predict.densityMclust(MC1_alt,Data2$X)))-MC2_null$loglik)
V_bar <- (V1+V2)/2

# Calculate P-values
P1 <- min(1/V1,1)
P2 <- min(1/V2,1)
P_bar <- min(1/V_bar,1)

# Print P-values
print(c(P1,P2,P_bar))

# LOOP --------------------------------------------------------------------

library(mclust)
library(MixSim)

# Seed
set.seed(69)

NN <- 1500 # n_1 (n_2) or n/2 sample size
DD <- 2 # Dimension of mixture model
BO <- 0.01 # Difference in mixtures for MixSim function
TrueG <- 5 # Generative number of components
LL <- 1 # additional components in alternative

# Simulation Settings
NumSim <- 100
GGmax <- 5
ResMat <- matrix(0,NumSim,GGmax) 

Fail_para_list <- list(); jj <- 1
Succ_para_list <- list(); kk <- 1

for (ii in 1:NumSim) {
  
  Try <- 1
  
  while (is.null(Try)==FALSE) {
    
    Try <- try({
      
      # Generate random mixture model parameters based on settings
      Sim_para <- MixSim(K=TrueG, p=DD, 
                         PiLow = 0.1, eps = 1e-03, BarOmega = BO)
      
      # Generate two data sets of sizes n_1 = n_2
      Data1 <- simdataset(n=NN,Pi=Sim_para$Pi,Mu=Sim_para$Mu,S=Sim_para$S)
      Data2 <- simdataset(n=NN,Pi=Sim_para$Pi,Mu=Sim_para$Mu,S=Sim_para$S)
      
      
      for (GG in 1:5) {
        
        # Fit mixtures under null hypothesis
        MC1_null <- densityMclust(Data1$X,G=GG,modelNames = 'VVV',verbose=F)
        MC2_null <- densityMclust(Data2$X,G=GG,modelNames = 'VVV',verbose=F)
        
        # Fit mixtures under the alternative hypothesis``
        MC1_alt <- densityMclust(Data1$X,G=GG+LL,modelNames = 'VVV',verbose=F)
        MC2_alt <- densityMclust(Data2$X,G=GG+LL,modelNames = 'VVV',verbose=F)
        
        # Calculate test statistics
        V1 <- exp(sum(log(predict.densityMclust(MC2_alt,Data1$X)))-MC1_null$loglik)
        V2 <- exp(sum(log(predict.densityMclust(MC1_alt,Data2$X)))-MC2_null$loglik)
        V_bar <- (V1+V2)/2
        
        # Calculate P-values
        P1 <- min(1/V1,1)
        P2 <- min(1/V2,1)
        P_bar <- min(1/V_bar,1)
        
        ResMat[ii,GG] <- P_bar
      }
      if (ResMat[ii,5] < 0.05) {
        Fail_para_list[[jj]] <- Sim_para
        jj <- jj + 1
      } else {
        Succ_para_list[[kk]] <- Sim_para
        kk <- kk + 1
      }
      cat("ii = ",ii,": ")
      cat(round(ResMat[ii,],5))
      cat("\n")
    })
  }
}

1-mean(ResMat[,5])

saveRDS(Fail_para_list, file= "Fail_para_list.rds") 
saveRDS(Succ_para_list, file= "Succ_para_list.rds") 

round(ResMat, 5)
