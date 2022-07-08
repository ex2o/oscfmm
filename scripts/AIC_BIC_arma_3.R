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

loglik_arma3 <- function(NN, DD, BO, TrueG, LL, NumSim, GGmax) { 
  
  p_mat <- matrix(0,GGmax,3)
  loglik <- rep(0,length=GGmax)
  result_list <- list()
  
  for (ii in 1:NumSim) {
    
    Try <- 1
    
    while (is.null(Try)==FALSE) {
      
      # # Generate random mixture model parameters based on settings
      Sim_para <- MixSim(
        BarOmega=BO, K=TrueG, p=DD, PiLow=1/(2*TrueG), eps=1e-03, sph=F)
      
      
      Try <- try({
        
        # Generate two data sets of sizes n_1 = n_2
        Data1 <- simdataset(n=NN,Pi=Sim_para$Pi,Mu=Sim_para$Mu,S=Sim_para$S)
        Data2 <- simdataset(n=NN,Pi=Sim_para$Pi,Mu=Sim_para$Mu,S=Sim_para$S)
        
        # Full configuration snapshot for storing with the results
        config <- list(
          NN=NN
          ,DD=DD
          ,BO=BO
          ,TrueG=TrueG
          ,LL=LL
          ,NumSim=NumSim
          ,GGmax=GGmax
          ,Sim_para=Sim_para
          ,Data1=Data1
          ,Data2=Data2
        )
        
        for (GG in 1:GGmax) {
          
          ts <- loglik_compute_ts(Data1, Data2, GG, LL, TrueG)
          
          loglik_GG <- ts$Log_lik_full
          
          loglik[GG] <- loglik_GG
          
        }
      })
    }
    
    result_list[[ii]] <- structure(0, config=config,
                                   Log_lik_full=loglik)
  }
  
  return(result_list)
}

loglik_compute_ts <- function(Data1, Data2, GG, LL, TrueG) {
  
  ##### Fit mixtures under the null hypothesis #####
  
  if (GG > TrueG) {
    
    K_means_init <- kmeans(Data2$X,
                           centers = GG,
                           iter.max = 100, 
                           nstar = 100)
    param_init <- me(modelName='VVV', data=Data2$X,
                     z=unmap(K_means_init$cluster),
                     control=emControl(itmax=1))
  } else {
    
    param_init <- me(modelName='VVV', data=Data2$X[Data2$id<=GG,],
                     z=unmap(Data2$id[Data2$id<=GG]),
                     control=emControl(itmax=1))
    
  }
  
  ##### Calculate log likelihood on combined data, for AIC and BIC #####
  Data_all_X <- rbind(Data1$X, Data2$X)
  
  Log_lik_full <- GMM_arma(t(Data_all_X),
                           param_init$parameters$pro, 
                           param_init$parameters$mean,
                           param_init$parameters$variance$sigma,
                           500,GG)$`log-likelihood` # fewer iterations, due to double sample size
  
  
  return(list(Log_lik_full = Log_lik_full))
  
  # Log_lik_MC2_alt_on_D1 = Log_lik_MC2_alt_on_D1
  # Log_lik_MC1_alt_on_D2 = Log_lik_MC1_alt_on_D2
  # Log_lik_MC1_null = MC1_null$`log-likelihood`
  # Log_lik_MC2_null = MC2_null$`log-likelihood`
  # Log_lik_MC1_alt = MC1_alt$`log-likelihood`
  # Log_lik_MC2_alt = MC2_alt$`log-likelihood`
}