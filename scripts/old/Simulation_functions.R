# Simulation functions -----------------------------------------------------------

mixsim <- function(config) {
  MixSim::MixSim(
    BarOmega=config$ms_grid_point$BO, 
    MaxOmega=config$ms_grid_point$MO, 
    eps=config$ms_grid_point$eps, 
    sph=config$sph,
    K=config$ms_grid_point$TrueG, 
    p=config$ms_grid_point$DD, 
    PiLow=config$ms_grid_point$PiLow)
}

simdata <- function(config) {
  MixSim::simdataset(
    n=config$ds_grid_point$NN, 
    Pi=config$sim_para$Pi, 
    Mu=config$sim_para$Mu, 
    S=config$sim_para$S)
}

# fit mixture models, and compute p-values 

one_g_step <- function(Data1, Data2, GG, config) {
  
  LL <- config$ds_grid_point$LL
  
  # Fit mixtures under null hypothesis
  MC1_null <- mclust::densityMclust(
    Data1$X, G=GG, modelNames=config$modelNames, verbose=F)
  MC2_null <- mclust::densityMclust(
    Data2$X, G=GG, modelNames=config$modelNames, verbose=F)
  
  # Fit mixtures under the alternative hypothesis
  MC1_alt <- mclust::densityMclust(
    Data1$X, G=GG+LL, modelNames=config$modelNames, verbose=F)
  MC2_alt <- mclust::densityMclust(
    Data2$X, G=GG+LL, modelNames=config$modelNames, verbose=F)
  
  # Calculate test statistics
  V1 <- exp(sum(log(
    mclust::predict.densityMclust(MC2_alt, Data1$X))) - MC1_null$loglik)
  V2 <- exp(sum(log(
    mclust::predict.densityMclust(MC1_alt, Data2$X))) - MC2_null$loglik)
  V_bar <- (V1+V2)/2
  
  # Calculate P-values
  P1 <- min(1/V1, 1)
  P2 <- min(1/V2, 1)
  P_bar <- min(1/V_bar, 1)
  
  return(c(P1, P2, P_bar))
}

# Take a draw from simdataset and perform one_procedure

one_procedure <- function(config, prev_results = list()) {
  
  if (!config$fixed_paras) {
    # Generate random mixture model parameters
    config$sim_para <- mixsim(config)
  }
  
  # Generate two data sets of sizes n_1 = n_2
  Data1 <- simdata(config)
  Data2 <- simdata(config)
  
  grid_pos <- unlist(c(config$ds_grid_point, 
                       config$ms_grid_point))
  if (config$verbose) {
    cat("*** Begin procedure ***","on",worker(),"\n")
  }
  
  config$Gmax <- config$ms_grid_point$TrueG + config$Gextra
  
  # Check that enough time remains to execute another procedure.
  config <- is_there_time(config)
  
  GG <- 1
  results <- matrix(NA, nrow=config$Gmax, ncol=3)
  while(config$there_is_time){
    
    if (config$verbose) {
      cat("Testing g =", GG, "against g =", 
          GG+config$ds_grid_point$LL, "at", date(),"on",worker(),"\n")
    }
    
    results[GG, ] <- retry_on_error({
      config$g_stepper(Data1, Data2, GG, config)
    }, max_errors=2, no_error=F)
    
    check <- all(results[GG,] >= 0.05) | (GG == config$Gmax)
    if (check) { break }
    
    config <- is_there_time(config)
    GG <- GG + 1
  }
  
  return(c(prev_results, list(
    structure(
      results, 
      grid_pos=grid_pos,
      sim_para=config$sim_para,
      there_is_time=config$there_is_time
    ))
  ))  
}

# For each ms_grid point, take a draw from MixSim and 
# repeat one_procedure for ds_draws times

one_mixture_grid <- function(config) {
  
  results <- list()
  
  create_rcpp_functions(config)
  
  for (i in 1:nrow(config$ms_grid)) {
    
    config$ms_grid_point <- config$ms_grid[i,]
    
    if (config$verbose) {
      cat("ms_grid_point:",
          paste0(names(config$ms_grid_point),"=",config$ms_grid_point),
          "on",worker(),"\n")
    }
    
    if (config$fixed_paras) {
      config$sim_para <- mixsim(config)
    }
    
    for (j in 1:nrow(config$ds_grid)) {
      
      config$ds_grid_point <- config$ds_grid[j,]
      if (config$verbose) {
        cat("ds_grid_point:",
            paste0(names(config$ds_grid_point),"=",config$ds_grid_point),
            "on",worker(),"\n")
      }
      
      for (k in 1:config$ds_draws) {
        
        config <- is_there_time(config)
        if (!config$there_is_time) {
          if(config$verbose){
            cat("<*><*><*><*> WE'RE OUT OF TIME! <*><*><*><*>\n")}
          return(results)
        }
        
        if (config$verbose) {
          cat("ds_draw",k,"of",config$ds_draws,
              "on",worker(),"\n")
        }
        results <- c(one_procedure(config), results)
        
      }
    }
  }
  
  results
}

# repeat one_mixture_grid for ms_draws times

simulation <- function(config) {
  
  if (config$parallel) {
    
    results <- do.call(foreach::foreach, config$foreach_args) %dopar% {
      retry_on_error({one_mixture_grid(config)}, max_errors=2, no_error=T,
                     config=config)
    }
    
    stopCluster(config$cl)
    
  } else {
    
    results <- list()
    for (i in 1:config$ms_draws) {
      
      if (config$verbose) {
        cat("-----------------------------------------------\n")
        cat("ms_draw", i, "of", config$ms_draws, "at", date(), "\n")
        cat("-----------------------------------------------\n")
      }
      this_result <- retry_on_error({
        one_mixture_grid(config)
      }, max_errors=2, no_error=T)
      results <- c(this_result, results)
      save_results(structure(results, config=config), config$slurm_task_id)
    }
  }
  
  structure(results, config=config,
            there_is_time=is_there_time(config)$there_is_time)
}


# Rcpp functions ----------------------------------------------------------

make_GMM_arma <- function() {
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
  cxxfunction(signature(data_r='numeric',
                        pi_r='numeric',
                        mean_r='numeric',
                        cov_r='numeric',
                        maxit_r='integer',
                        groups_r='integer'),
              gmm_full_src, plugin = 'RcppArmadillo')
  
}

# An rcpp alternative to one_g_step. Creating these inside this function 
# is a work around for Rcpp parallel processing problems related to exporting 
# these functions that depend on GMM_arma to all cores from foreach.
# These functions were not able to be found unless they were
# created on the cores.

create_rcpp_functions <- function(config) {
  
  if (config$use_rcpp) {
    
    if (config$verbose) {cat("Making GMM_arma on",worker(),"\n")}
    
    assign("GMM_arma", make_GMM_arma(), envir=.GlobalEnv)
    
    fit_mixture_rcpp <- function(Data, CC, config) {
      
      # Initialize a clustering using k-means
      K_means_init <- kmeans(Data$X, centers = CC,
                             iter.max=100, nstar=100)
      
      # Use mclust to convert kmeans to parameters
      param_init <- mclust::meVVV(Data$X,
                                  z=mclust::unmap(K_means_init$cluster),
                                  control=mclust::emControl(itmax=1))
      
      # Run GMM using arma
      MC <- GMM_arma(t(Data$X),
                     param_init$parameters$pro,
                     param_init$parameters$mean,
                     param_init$parameters$variance$sigma,
                     2000, CC)
      MC
    }
    
    assign("fit_mixture_rcpp",fit_mixture_rcpp,envir=.GlobalEnv)
    
    compute_likelihood_rcpp <- function(Data, MC_null, MC_alt, CC) {
      
      # Use arma code to calculate likelihood at specified parameters (i.e., run with zero iterations)
      Log_lik <- GMM_arma(t(Data$X),
                          MC_alt$proportions, 
                          MC_alt$means,
                          MC_alt$covariances,
                          0, CC)$`log-likelihood` 
      
      # Compute likelihood
      exp(Log_lik - MC_null$`log-likelihood`)
    }
    
    assign("compute_likelihood_rcpp",compute_likelihood_rcpp,envir=.GlobalEnv)
    
  }
}

one_g_step_rcpp <- function(Data1, Data2, GG, config) {
  
  LL <- config$ds_grid_point$LL
  
  # Fit mixtures under null hypothesis
  MC1_null <- fit_mixture_rcpp(Data1, GG, config)
  MC2_null <- fit_mixture_rcpp(Data2, GG, config)
  
  # Fit mixtures under the alternative hypothesis
  MC1_alt  <- fit_mixture_rcpp(Data1, GG+LL, config)
  MC2_alt  <- fit_mixture_rcpp(Data2, GG+LL, config)
  
  # Calculate test statistics
  V1 <- compute_likelihood_rcpp(Data1, MC1_null, MC2_alt, GG+LL)
  V2 <- compute_likelihood_rcpp(Data2, MC2_null, MC1_alt, GG+LL)
  
  V_bar <- (V1+V2)/2
  
  # Calculate P-values
  P1 <- min(1/V1,1)
  P2 <- min(1/V2,1)
  P_bar <- min(1/V_bar,1)
  
  return(c(P1, P2, P_bar))
}
