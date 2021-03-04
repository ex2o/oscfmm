# Libaries
library(mclust)
library(MixSim)
library(Rcpp)
library(RcppArmadillo)
library(inline)
library(profvis)

# Load ARMA code
gmm_full_src <- '
using namespace arma;
// Convert necessary matrix objects to arma
mat data_a = as<mat>(data_r);
rowvec pi_a = as<rowvec>(pi_r);
mat mean_a = as<mat>(mean_r);
cube cov_a = as<cube>(cov_r);
int maxit_a = as<int>(maxit_r);
int groups_a = as<int>(groups_r);
// Initialize a gmm_full object
gmm_full model;
// Set the parameters
model.set_params(mean_a, cov_a, pi_a);
model.learn(data_a, groups_a, maha_dist, keep_existing, 0, maxit_a, 2.2e-16, false);
// 
return Rcpp::List::create(
  Rcpp::Named("log-likelihood")=model.sum_log_p(data_a),
  Rcpp::Named("proportions")=model.hefts,
  Rcpp::Named("means")=model.means,
  Rcpp::Named("covariances")=model.fcovs);
'
GMM_arma <- cxxfunction(signature(data_r='numeric',
                                  pi_r='numeric',
                                  mean_r='numeric',
                                  cov_r='numeric',
                                  maxit_r='integer',
                                  groups_r='integer'),
                        gmm_full_src, plugin = 'RcppArmadillo')

# Seed
set.seed(9999)

NN <- 10000 # n_1 (n_2) or n/2 sample size
DD <- 2 # Dimension of mixture model
BO <- 0.01 # Difference in mixtures for MixSim function
TrueG <- 5 # Generative number of components
LL <- 1 # additional components in alternative

# Simulation Settings
NumSim <- 100
GGmax <- 5
ResMat <- matrix(0,NumSim,GGmax) 
# Generate random mixture model parameters based on settings
# Sim_para <- MixSim(BarOmega = BO,K=TrueG,p=DD,PiLow = 1/(2*TrueG),eps = 1e-03,int=c(-10,10))

for (ii in 1:NumSim) {
  
  Try <- 1
  
  while (is.null(Try)==FALSE) {
    
    # # Generate random mixture model parameters based on settings
    Sim_para <- MixSim(BarOmega = BO,K=TrueG,p=DD,PiLow = 1/(2*TrueG),eps = 1e-03,int=c(-100,100))
    
    Try <- try({
      
      # Generate two data sets of sizes n_1 = n_2
      Data1 <- simdataset(n=NN,Pi=Sim_para$Pi,Mu=Sim_para$Mu,S=Sim_para$S)
      Data2 <- simdataset(n=NN,Pi=Sim_para$Pi,Mu=Sim_para$Mu,S=Sim_para$S)
      
      
      for (GG in 1:GGmax) {
        
        # Fit mixtures under null hypothesis
        K_means_init <- kmeans(Data1$X,centers = GG,iter.max = 100, nstar=100) # Initialize a clustering using k-means
        param_init <- me('VVV',Data1$X,z=unmap(K_means_init$cluster),control=emControl(itmax=1)) # Use mclust to convert kmeans to parameters
        MC1_null <- GMM_arma(t(Data1$X),param_init$parameters$pro, param_init$parameters$mean,
                             param_init$parameters$variance$sigma,2000,GG) # Run GMM using arma
        
        K_means_init <- kmeans(Data2$X,centers = GG,iter.max = 100, nstar=100)
        param_init <- me('VVV',Data2$X,z=unmap(K_means_init$cluster),control=emControl(itmax=1))
        MC2_null <- GMM_arma(t(Data2$X),param_init$parameters$pro, param_init$parameters$mean,
                             param_init$parameters$variance$sigma,2000,GG)
        
        # Fit mixtures under the alternative hypothesis
        K_means_init <- kmeans(Data1$X,centers = GG+LL,iter.max = 100, nstar=100)
        param_init <- me('VVV',Data1$X,z=unmap(K_means_init$cluster),control=emControl(itmax=1))
        MC1_alt <- GMM_arma(t(Data1$X),param_init$parameters$pro, param_init$parameters$mean,
                            param_init$parameters$variance$sigma,2000,GG+LL)
        
        K_means_init <- kmeans(Data2$X,centers = GG+LL,iter.max = 100, nstar=100)
        param_init <- me('VVV',Data2$X,z=unmap(K_means_init$cluster),control=emControl(itmax=1))
        MC2_alt <- GMM_arma(t(Data2$X),param_init$parameters$pro, param_init$parameters$mean,
                            param_init$parameters$variance$sigma,2000,GG+LL)
        
        # Calculate test statistics
        # Use arma code to calculate likelihood at specified parameters (i.e., run with zero iterations)
        Log_lik_MC2_alt_on_D1 <- GMM_arma(t(Data1$X),MC2_alt$proportions, MC2_alt$means,
                                          MC2_alt$covariances,0,GG+LL)$`log-likelihood` 
        # Compute likelihood
        V1 <- exp(Log_lik_MC2_alt_on_D1 - MC1_null$`log-likelihood`) # Compute
        Log_lik_MC1_alt_on_D2 <- GMM_arma(t(Data2$X),MC1_alt$proportions, MC1_alt$means,
                                          MC1_alt$covariances,0,GG+LL)$`log-likelihood`
        V2 <- exp(Log_lik_MC1_alt_on_D2 - MC2_null$`log-likelihood`)
        V_bar <- (V1+V2)/2
        
        # Calculate P-values
        P1 <- min(1/V1,1)
        P2 <- min(1/V2,1)
        P_bar <- min(1/V_bar,1)
        
        ResMat[ii,GG] <- P_bar
        
        #print(c(ii,GG,P_bar))
        
        if (GG==5 && P_bar < 0.05) {
          plot(rbind(Data1$X,Data2$X),col=c(Data1$id,Data2$id),main=ii)
        }
      }
      cat("ii = ",ii,": ")
      cat(round(ResMat[ii,],5))
      cat("\n")
    })
  }
  
}

