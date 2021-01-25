
# Utilities ---------------------------------------------------------------
load_packages <- Vectorize(function(package) {
  if (!require(package, character.only = T)) {
    install.packages(package)
    require(package, character.only = T)
    return(F)
  }
  return(T)
})

combine_print_save <- function (a, ...) {
  message("combining: ",length(a))
  saveRDS(a,file=paste0("result_backup",length(a),".rds"))
  c(a, list(...))
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

one_procedure <- function(Sim_para, LL, NN, Gmax, stop) {
  
  TrueG <- length(Sim_para$Pi)
  
  # Generate two data sets of sizes n_1 = n_2
  Data1 <- MixSim::simdataset(n=NN,Pi=Sim_para$Pi,Mu=Sim_para$Mu,S=Sim_para$S)
  Data2 <- MixSim::simdataset(n=NN,Pi=Sim_para$Pi,Mu=Sim_para$Mu,S=Sim_para$S)
  
  results <- matrix(NA, nrow = Gmax, ncol = 3)
  if (stop) {
    GG <- 1
    repeat{
      results[GG, ] <- one_g_step(Data1, Data2, GG, LL)
      if (all(results[GG,] >= 0.05)) { break }
      GG <- GG + 1
    }
  } else {
    for (GG in 1:Gmax) {
      cat("Testing g =",GG,"against g=",GG+LL,"\n")
      results[GG, ] <- one_g_step(Data1, Data2, GG, LL)
    }
  }
  
  return(results)  
}

sim_NN <- function(NN_range, DD, BO, TrueG, LL, Gmax, 
                   parallel = T, stop, free_core) {
  
  # Generate random mixture model parameters based on settings
  sim_para <- MixSim::MixSim(BarOmega = BO, K=TrueG, p=DD, PiLow = 0.1)
  
  results <- list()
  for(i in 1:length(NN_range)) {
    results[[i]] <- one_procedure(sim_para, LL, NN_range[i], Gmax, stop)
  }
  
  return(structure(results, sim_para = sim_para))
}

repeat_sim <- function(sim, params, reps, save, name, parallel) {
  res <- list()
  
  if (parallel) {
    # Activate cluster
    no_cores <- parallel::detectCores()
    if (free_core) {no_cores <- no_cores - 1}
    doParallel::registerDoParallel(cores=no_cores)  
    cl <- parallel::makeCluster(no_cores)  
    message(paste0("Num cores = ", no_cores,"\n"))
    print(cl)
    
    res <- foreach::foreach(
      i = 1:reps, .inorder = FALSE, .combine = combine_print_save,
      .export = c("one_procedure", "one_g_step", 
                  deparse(substitute(sim)))) %dopar% {
        do.call(sim, params)
      }
    
  } else {
    
    for (i in 1:reps) {
      message(paste0("non-parallel: repeat_sim iteration ",i," of ",reps))
      res[[i]] <- do.call(sim, params)
      if (save) {
        saveRDS(res, file = paste0(name,"_","res",i,".rds"))
      }
    }
  }
    
  results <- structure(
    res,
    params = params,
    reps = reps)
  if (save) {
    saveRDS(results, file = paste0(name,"_","results.rds"))
  }
  return(results)
}

# Summary functions -------------------------------------------------------

first_reject <- function(res) {
  lapply(res, function(x){
    apply(x, MARGIN = 2, FUN = function(x){which(x>=0.05)[1]})
  })
}