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

