# Utilities ---------------------------------------------------------------
load_packages <- Vectorize(function(package) {
  if (!require(package, character.only = T)) {
    install.packages(package)
    require(package, character.only = T)
    return(F)
  }
  return(T)
})

combine_print_save <- function(a, ...) {
  message("combining: ",length(a))
  saveRDS(a,file=paste0("result_backup",length(a),".rds"))
  c(a, list(...))
}

create_config <- function(...) {
  args <- list(...)
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


# Sim functions -----------------------------------------------------------

one_g_step <- function(Data1, Data2, GG, LL) {
  # Fit mixtures under null hypothesis
  MC1_null <- mclust::densityMclust(Data1$X,G=GG,modelNames = 'VVV')
  MC2_null <- mclust::densityMclust(Data2$X,G=GG,modelNames = 'VVV')
  
  # Fit mixtures under the alternative hypothesis
  MC1_alt <- mclust::densityMclust(Data1$X,G=GG+LL,modelNames = 'VVV')
  MC2_alt <- mclust::densityMclust(Data2$X,G=GG+LL,modelNames = 'VVV')
  
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

one_procedure <- function(config) {
  
  Gmax <- config$TrueG + config$Gextra
  
  # Generate random mixture model parameters based on settings
  sim_para <- MixSim::MixSim(
    BarOmega = config$BO, K=config$TrueG, p=config$DD, PiLow = config$PiLow)
  
  # Generate two data sets of sizes n_1 = n_2
  Data1 <- MixSim::simdataset(
    n=config$NN, Pi=Sim_para$Pi, Mu=Sim_para$Mu, S=Sim_para$S)
  Data2 <- MixSim::simdataset(
    n=config$NN, Pi=Sim_para$Pi, Mu=Sim_para$Mu, S=Sim_para$S)
  
  results <- matrix(NA, nrow = Gmax, ncol = 3)
  
  if (config$stop) {
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
  
  return(results)  
}

sim_1D <- function(config, grid_name) {
  
  results <- list()
  grid <- config[[grid_name]]
  for (x in grid) {
    configx <- config
    configx[[grid_name]] <- x
    results <- c(one_procedure(config), results)
  }
}

sim_NN <- function(config) {
  
  
  results <- list()
  for(i in 1:length(config$N)) {
    results[[i]] <- one_procedure(config)
  }
  
  return(structure(results, sim_para = sim_para))
}

sim_BO <- function(config) {
  
  results <- list()
  for(i in 1:length(config$BO)) {
    results[[i]] <- one_procedure(config)
  }
  
  return(structure(results, sim_para = sim_para))
}

sim_TrueG <- function(config) {
  
  results <- list()
  for(i in 1:length(config$TrueG)) {
    results[[i]] <- one_procedure(config)
  }
  
  return(structure(results, sim_para = sim_para))
}

choose_1D_grid_sim <- function(config) {
  
  grids <- which(lapply(config, length) > 1)
  if (length(grids) > 1) {
    stop("Check showed more than one element in config has length ",
         "greater than 1.",
         "There should only be one grid variable in config ",
         "within choose_1D_grid_sim. ",
         "Try using repeat_sim instead.")
  }
  
  grid <- names(grids)[1]
  switch(grid,
         "NN" = sim_NN,
         "BO" = sim_BO,
         "TrueG" = sim_TrueG,
         stop(grid, " grid sim not yet implemented"))
}

repeat_sim <- function(config, result = list()) {
  
  grids <- which(lapply(config, length) > 1)
  if (length(grids) == 1) {
    
    result <- repeat_1D_sim(config)
    
  } else {
    
    first_grid <- config[[grids[1]]]
    for( x in first_grid ) {
      
      configx <- config
      configx[[grids[1]]] <- x
      result <- c(repeat_sim(configx, result), result)
    }
  }
  
  result
}

repeat_1D_sim <- function(config) {
  
  res <- list()
  if (config$parallel) {
    
    # Activate cluster
    no_cores <- parallel::detectCores()
    if (config$free_core) {no_cores <- no_cores - 1}
    doParallel::registerDoParallel(cores=no_cores)  
    cl <- parallel::makeCluster(no_cores)  
    message("Num cores = ", no_cores,"\n")
    print(cl)
    
    res <- foreach::foreach(i = 1:config$reps, 
      .init = list(), .inorder = FALSE, .combine = combine_print_save,
      .export = c("one_procedure", "one_g_step")) %dopar% {
        sim_1D(config, grid_name)
      }
    
    stopCluster(cl)
    
  } else {
    
    for (i in 1:config$reps) {
      message("non-parallel: repeat_sim iteration ",i," of ",config$reps)
      res[[i]] <- sim_1D(config, grid_name)
      if (config$save) {
        saveRDS(res, file = paste0(config$name,"_","res",i,".rds"))
      }
    }
  }
    
  results <- structure(res, config = config)
  if (config$save) {
    saveRDS(results, file = paste0(config$name,"_","results.rds"))
  }
  return(results)
}

# Summary functions -------------------------------------------------------

first_reject <- function(res) {
  lapply(res, function(x){
    apply(x, MARGIN = 2, FUN = function(x){which(x>=0.05)[1]})
  })
}