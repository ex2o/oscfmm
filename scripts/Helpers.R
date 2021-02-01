# Utilities ---------------------------------------------------------------
load_packages <- Vectorize(function(package) {
  if (!require(package, character.only = T)) {
    install.packages(package)
    require(package, character.only = T)
    return(F)
  }
  return(T)
})

combiner <- function(a, ...) {
  #message("combining: ",length(a))
  #saveRDS(a,file=paste0("result_backup",length(a),".rds"))
  #c(a, list(...))
  c(a, ...)
}

create_config <- function(...) {
  config <- list(...)
  config$grids <- list(names(which(lapply(config, length) > 1)))
  visualise_config(config)
  config
}

# Assuming no more than 100 simultaneous jobs on SLURM
save_results <- function(results) {
  for (i in 1:100) {
    file <- paste0("results",i,".rds")
    if (!file.exists(file)) {
      saveRDS(results, file = file)
      return("Success")
    }
  }
  error("Unable to continue saving results")
}

# only works if config has no grids
visualise_config <- function(config) {
  ngrids <- length(config$grids[[1]])
  if(ngrids != 0) {
    return(F)
  }
  
  # Generate random mixture model parameters based on settings
  Sim_para <- MixSim(BarOmega = config$BO, K=config$TrueG, p=config$DD, 
                     PiLow = config$PiLow)
  
  # Generate example data
  Data <- simdataset(n=config$NN,Pi=Sim_para$Pi,Mu=Sim_para$Mu,S=Sim_para$S)
  
  plot(Data$X, col = Data$id)
  
  return(T)
}


# Sim functions -----------------------------------------------------------

one_g_step <- function(Data1, Data2, GG, LL) {
  # Fit mixtures under null hypothesis
  MC1_null <- mclust::densityMclust(Data1$X,G=GG,modelNames = 'VVV',verbose=F)
  MC2_null <- mclust::densityMclust(Data2$X,G=GG,modelNames = 'VVV',verbose=F)
  
  # Fit mixtures under the alternative hypothesis
  MC1_alt <- mclust::densityMclust(Data1$X,G=GG+LL,modelNames = 'VVV',verbose=F)
  MC2_alt <- mclust::densityMclust(Data2$X,G=GG+LL,modelNames = 'VVV',verbose=F)
  
  # Calculate test statistics
  V1 <- exp(sum(log(
    mclust::predict.densityMclust(MC2_alt,Data1$X)))-MC1_null$loglik)
  V2 <- exp(sum(log(
    mclust::predict.densityMclust(MC1_alt,Data2$X)))-MC2_null$loglik)
  V_bar <- (V1+V2)/2
  
  # Calculate P-values
  P1 <- min(1/V1,1)
  P2 <- min(1/V2,1)
  P_bar <- min(1/V_bar,1)
  
  return(c(P1,P2,P_bar))
}

one_procedure <- function(config, prev_results = list()) {
  
  Gmax <- config$TrueG + config$Gextra
  
  # Generate random mixture model parameters based on settings
  sim_para <- MixSim::MixSim(
    BarOmega = config$BO, K=config$TrueG, p=config$DD, PiLow = config$PiLow)
  
  # Generate two data sets of sizes n_1 = n_2
  Data1 <- MixSim::simdataset(
    n=config$NN, Pi=sim_para$Pi, Mu=sim_para$Mu, S=sim_para$S)
  Data2 <- MixSim::simdataset(
    n=config$NN, Pi=sim_para$Pi, Mu=sim_para$Mu, S=sim_para$S)
  
  results <- matrix(NA, nrow = Gmax, ncol = 3)
  if (config$stop_on_accept) {
    
    GG <- 1
    repeat{
      results[GG, ] <- one_g_step(Data1, Data2, GG, config$LL)
      if ( all(results[GG,] >= 0.05) | GG == Gmax) { break }
      GG <- GG + 1
    }
  } else {
    
    for (GG in 1:Gmax) {
      cat("Testing g =",GG,"against g=",GG+config$LL,"\n")
      results[GG, ] <- one_g_step(Data1, Data2, GG, config$LL)
    }
  }
  
  return(c(prev_results, list(structure(
    results, config = unlist(config[config$grids[[1]]])))))  
}

sim_recursive <- function(config, prev_results = list()) {
  
  if (config$verbose) {
    cat("Went one deeper at",date(),"with length(prev_results) = ",
        length(prev_results),"\n")}
  
  grids <- which(lapply(config, length) > 1)
  if (length(grids) == 1) {
    next_sim <- one_procedure
  } else {
    next_sim <- sim_recursive
  }
  
  next_grid <- config[[grids[1]]]
  results <- prev_results
  for( x in next_grid ) {
    
    configx <- config
    configx[[grids[1]]] <- x
    results <- next_sim(configx, results)
  }
  
  if (config$verbose) {
    cat("Unwound a layer at",date(),"with length(prev_results) = ",
        length(prev_results),"\n")}
  
  return(results)
}


# Sim repetition functions ------------------------------------------------

repeat_sim <- function(config, sim) {
  
  results <- list()
  
  if (config$parallel) {
  
    # Prepare foreach loop arguments
    foreach_config <- list(
      i = 1:config$reps 
      ,.inorder = FALSE 
      ,.combine = combiner, .init = list()
      ,.export = c("sim_recursive","one_procedure", "one_g_step")
    )
    
    # Activate cluster
    no_cores <- parallel::detectCores()
    if (config$free_core) {no_cores <- no_cores - 1}
    doParallel::registerDoParallel(cores=no_cores)
    if (config$verbose) {
      foreach_config$.verbose = T
      cl <- parallel::makeCluster(no_cores, outfile = "foreach_log.txt")
    } else {
      cl <- parallel::makeCluster(no_cores)
    }
      
    message("Num cores = ", no_cores,"\n")
    print(cl)

    # (parallel) foreach loop
    results <- do.call(foreach::foreach, foreach_config) %dopar% {
        sim(config)
    }
    
    stopCluster(cl)
    
  } else {
    
    for (i in 1:config$reps) {
      results <- c(sim(config), results)
    }
  }
  
  structure(results, config = config)
}


# Summary functions -------------------------------------------------------

first_accept  <- function(res) {
  lapply(res, function(x){
    apply(x, MARGIN = 2, FUN = function(x){which(x>=0.05)[1]})
  })
}