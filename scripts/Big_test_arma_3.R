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

arma3 <- function(NN, DD, BO, TrueG, LL, NumSim, GGmax) { 
  
  p_mat <- matrix(0,GGmax,3)
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
    
          ts <- compute_ts(Data1, Data2, GG, LL, TrueG)
          
          # Calculate P-values
          P1 <- min(1/ts$V1,1)
          P2 <- min(1/ts$V2,1)
          P_bar <- min(1/ts$V_bar,1)
          
          p_mat[GG,] <- c(P1, P2, P_bar)
          print(as.integer(c(ii,GG,P_bar)))
          
        }
      })
    }
    
    result_list[[ii]] <- structure(p_mat, config=config,
                                   Log_lik_MC2_alt_on_D1 = ts$Log_lik_MC2_alt_on_D1,
                                   Log_lik_MC1_alt_on_D2 = ts$Log_lik_MC1_alt_on_D2,
                                   Log_lik_MC1_null = ts$Log_lik_MC1_null,
                                   Log_lik_MC2_null = ts$Log_lik_MC2_null,
                                   Log_lik_MC1_alt  = ts$Log_lik_MC1_alt,
                                   Log_lik_MC2_alt  = ts$Log_lik_MC2_alt)
  }
  
  return(result_list)
}



compute_ts <- function(Data1, Data2, GG, LL, TrueG) {
  
  ##### Fit mixtures under the null hypothesis #####
  
  if (GG > TrueG) {
    
    K_means_init <- kmeans(Data1$X,
                           centers = GG,
                           iter.max = 100, 
                           nstar = 100)
    param_init <- me(modelName='VVV', data=Data1$X,
                     z=unmap(K_means_init$cluster),
                     control=emControl(itmax=1))
  } else {
    
    # use mclust to convert kmeans to parameters 
    param_init <- me(modelName='VVV', data=Data1$X[Data1$id<=GG,],
                     z=unmap(Data1$id[Data1$id<=GG]),
                     control=emControl(itmax=1)) 
    
  }
  
  MC1_null <- GMM_arma(t(Data1$X),
                       param_init$parameters$pro, 
                       param_init$parameters$mean,
                       param_init$parameters$variance$sigma,
                       2000,GG) # Run GMM using arma
  
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
  
  MC2_null <- GMM_arma(t(Data2$X),
                       param_init$parameters$pro, 
                       param_init$parameters$mean,
                       param_init$parameters$variance$sigma,
                       2000,GG)
  

  ##### Fit mixtures under the alternative hypothesis #####
  
  if (GG + LL > TrueG) {
    K_means_init <- kmeans(Data1$X,
                           centers = GG+LL,
                           iter.max = 100, 
                           nstar = 100)
    param_init <- me(modelName='VVV', data=Data1$X,
                     z=unmap(K_means_init$cluster),
                     control=emControl(itmax=1))
  } else {
     
    
    param_init <- me(modelName='VVV', data=Data1$X[Data1$id<=(GG+LL),],
                     z=unmap(Data1$id[Data1$id<=(GG+LL)]),
                     control=emControl(itmax=1)) 
    
  }
  
  MC1_alt <- GMM_arma(t(Data1$X),
                      param_init$parameters$pro, 
                      param_init$parameters$mean,
                      param_init$parameters$variance$sigma,
                      2000,GG+LL)
  

  
  if (GG + LL > TrueG) {
    
    K_means_init <- kmeans(Data2$X,
                           centers = GG+LL,
                           iter.max = 100, 
                           nstar = 100)
    param_init <- me(modelName='VVV', data=Data2$X,
                     z=unmap(K_means_init$cluster),
                     control=emControl(itmax=1))
  } else {
    
    param_init <- me(modelName='VVV', data=Data2$X[Data2$id<=(GG+LL),],
                     z=unmap(Data2$id[Data2$id<=(GG+LL)]),
                     control=emControl(itmax=1))  
  }
  
  MC2_alt <- GMM_arma(t(Data2$X),
                      param_init$parameters$pro, 
                      param_init$parameters$mean,
                      param_init$parameters$variance$sigma,
                      2000,GG+LL)
  
  ##### Calculate test statistics #####
  
  # Use arma code to calculate likelihood at specified parameters 
  # (i.e., run with zero iterations)
  Log_lik_MC2_alt_on_D1 <- GMM_arma(t(Data1$X),
                                    MC2_alt$proportions, 
                                    MC2_alt$means,
                                    MC2_alt$covariances,
                                    0,GG+LL)$`log-likelihood` 
  # Compute likelihood
  V1 <- exp(Log_lik_MC2_alt_on_D1 - MC1_null$`log-likelihood`)
  
  Log_lik_MC1_alt_on_D2 <- GMM_arma(t(Data2$X),
                                    MC1_alt$proportions, 
                                    MC1_alt$means,
                                    MC1_alt$covariances,
                                    0,GG+LL)$`log-likelihood`
  V2 <- exp(Log_lik_MC1_alt_on_D2 - MC2_null$`log-likelihood`)
  V_bar <- (V1+V2)/2
  
  return(list(V1 = V1, V2 = V2, V_bar = V_bar,
              Log_lik_MC2_alt_on_D1 = Log_lik_MC2_alt_on_D1,
              Log_lik_MC1_alt_on_D2 = Log_lik_MC1_alt_on_D2,
              Log_lik_MC1_null = MC1_null$`log-likelihood`,
              Log_lik_MC2_null = MC2_null$`log-likelihood`,
              Log_lik_MC1_alt = MC1_alt$`log-likelihood`,
              Log_lik_MC2_alt = MC2_alt$`log-likelihood`))
}
