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
  b <- c(a, ...)
  save_results(b)
  b
}

create_config <- function(...) {
  
  config <- list(...)
  config$modelNames <- if (config$sph) {
    "VII"
    message("spherical model only")
  } else {
    "VVV"
  }
  

  cat("---- Mixture draw Grid ----\n")
  print(config$ms_grid)
  cat("\n")
  cat("---- Dataset draw Grid ----\n")
  print(config$ds_grid)
  cat("\n")

  
  if (config$parallel) {

    # Prepare foreach loop arguments
    config$foreach_args <- list(
      i = 1:config$ms_draws
      ,.inorder = FALSE
      ,.combine = combiner, .init = list()
      ,.export = c("one_mixture_grid","one_procedure", "one_g_step",
                   "retry_on_error")
    )

    # Detect, modify and check number of cores
    config$no_cores <- parallel::detectCores()
    if (config$free_core) {
      config$no_cores <- config$no_cores - 1
    }
    if (config$ms_draws %% config$no_cores != 0) {
      warning(
        "ms_draws=", config$ms_draws,
        " is not a multiple of no_cores=", config$no_cores,
        ", consider increasing ms_draws for free")
    }

    # Activate the cluster
    doParallel::registerDoParallel(cores=config$no_cores)
    if (config$verbose) {
      config$foreach_args$.verbose <- T
      config$cl <- parallel::makeCluster(
        config$no_cores,
        outfile = "foreach_log.txt")
    } else {
      config$cl <- parallel::makeCluster(
        config$no_cores)
    }
    cat("Parallelising at ms_draws level\n")
    cat(capture.output(print(config$cl)),"\n")
  }
  config
}

# Assuming no more than 100 simultaneous jobs on SLURM
save_results <- function(results) {
  for (i in 1:1000) {
    file <- paste0("results",i,".rds")
    if (!file.exists(file)) {
      saveRDS(results, file = file)
      return("Success")
    }
  }
  error("Unable to continue saving results")
}

retry_on_error <- function(
  expr, 
  isError=function(x) "try-error" %in% class(x), 
  maxErrors=10) {
  
  attempts <- 0
  retval <- try(eval(expr))
  while (isError(retval)) {
    attempts <- attempts + 1
    if (attempts >= maxErrors) {
      msg <- sprintf("retry: too many retries [[%s]]\n", 
                     retval[1])
      stop(msg)
    } else {
      msg <- sprintf("retry: error in attempt %i/%i [[%s]]\n", 
                     attempts, 
                     maxErrors, 
                     retval[1])
      warning(msg)
    }
    retval <- try(eval(expr))
  }
  return(retval)
}

# Sim functions -----------------------------------------------------------

# fit mixture models, and compute p-values 

one_g_step <- function(Data1, Data2, GG, LL) {
  
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
    config$sim_para <- MixSim::MixSim(
      BarOmega=config$ms_grid_point$BO, 
      MaxOmega=config$ms_grid_point$MO, 
      eps=config$ms_grid_point$eps, 
      sph=config$ms_grid_point$sph,
      K=config$ms_grid_point$TrueG, 
      p=config$ms_grid_point$DD, 
      PiLow=config$ms_grid_point$PiLow)
  }
  
  # Generate two data sets of sizes n_1 = n_2
  Data1 <- MixSim::simdataset(
    n=config$ds_grid_point$NN, 
    Pi=config$sim_para$Pi, 
    Mu=config$sim_para$Mu, 
    S=config$sim_para$S)
  Data2 <- MixSim::simdataset(
    n=config$ds_grid_point$NN, 
    Pi=config$sim_para$Pi, 
    Mu=config$sim_para$Mu, 
    S=config$sim_para$S)
  
  grid_pos <- unlist(c(config$ds_grid_point, 
                       config$ms_grid_point))
  if (config$verbose) {
    cat("*** Begin procedure ***\n")
  }
  
  Gmax <- config$ms_grid_point$TrueG + config$Gextra
  results <- matrix(NA, nrow=Gmax, ncol=3)
  GG <- 1
  repeat{
    if (config$verbose) {
      cat("Testing g =", GG, "against g =", 
          GG+config$ds_grid_point$LL, "at", date(), "\n")
    }
    
    results[GG, ] <- retry_on_error({
      one_g_step(Data1, Data2, GG, config$ds_grid_point$LL)
    })
    
    check <- all(results[GG,] >= 0.05) | (GG == Gmax)
    if (check) { break }
    GG <- GG + 1
  }
  
  return(c(prev_results, list(
    structure(
      results, 
      grid_pos=grid_pos,
      sim_para=config$sim_para
    ))
  ))  
}

# For each ms_grid point, take a draw from MixSim and 
# repeat one_procedure for ds_draws times

one_mixture_grid <- function(config) {
  
  results <- list()
  
  for (i in 1:nrow(config$ms_grid)) {
    
    config$ms_grid_point <- config$ms_grid[i,]
    
    if (config$verbose) {
      cat("ms_grid_point:",
          paste0(names(config$ms_grid_point),"=",config$ms_grid_point),"\n")
    }
    
    if (config$fixed_paras) {
      config$sim_para <- MixSim::MixSim(
        BarOmega=config$ms_grid_point$BO, 
        MaxOmega=config$ms_grid_point$MO, 
        eps=config$ms_grid_point$eps, 
        sph=config$sph,
        K=config$ms_grid_point$TrueG, 
        p=config$ms_grid_point$DD, 
        PiLow=config$ms_grid_point$PiLow)
    }
    
    for (j in 1:nrow(config$ds_grid)) {

      config$ds_grid_point <- config$ds_grid[j,]
      if (config$verbose) {
        cat("ds_grid_point:",
            paste0(names(config$ds_grid_point),"=",config$ds_grid_point),"\n")
      }
      
      for (k in 1:config$ds_draws) {
        
        if (config$verbose) {
          cat("ds_draw",k,"of",config$ds_draws,"\n")
        }
        results <- c(one_procedure(config), results)
        
      }
    }
  }
  
  results
}


simulation <- function(config) {
  
  if (config$parallel) {
    
    results <- do.call(foreach::foreach, config$foreach_args) %dopar% {
      one_mixture_grid(config)
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
      results <- c(one_mixture_grid(config), results)
      save_results(structure(results, config=config))
    }
  }
  
  results
}

# Summary functions -------------------------------------------------------

first_accept  <- function(res) {
  lapply(res, function(x){
    if (any(x < -99, na.rm=T)) {return(NA)}
    apply(x, MARGIN=2, FUN=function(x){which(x >= 0.05)[1]})
  })
}

# Reshape and combine results for first accept, in long format
combine_first_accepts <- function(results) {
  grid <- results %>% 
    lapply(function(x){attr(x,"grid_pos")}) %>% 
    do.call(rbind, .)
  res <- results %>% 
    first_accept %>% 
    do.call(rbind, .) %>% 
    data.frame() %>% 
    set_names(c("p1","p2","p3"))
  res <- cbind(res,grid)
  res
}

# Add a quantile interval for each grid point
add_quantiles <- function(results) {
  res <- results %>% 
    group_by(across(!starts_with("p"))) %>% 
    mutate(p1_qL = quantile(p1, probs = 0.025)
           ,p1_qU = quantile(p1, probs = 0.975)
           ,p2_qL = quantile(p2, probs = 0.025)
           ,p2_qU = quantile(p2, probs = 0.975)
           ,p3_qL = quantile(p3, probs = 0.025)
           ,p3_qU = quantile(p3, probs = 0.975))
  res
}


# Old functions -----------------------------------------------------------

## sim_recursive just jumps into the config, identifies grids,
# and repeats the procedure for each possibly combination of
# grid elements. It isn't designed to work with a grid that
# contains ms_grid and ds_grid

sim_recursive <- function(config, prev_results=list()) {
  
  grids <- which(lapply(config, length) > 1)
  
  if(length(grids) == 0) {
    stop("sim_re cursive needs at least one grid")
  }
  
  next_sim <- ifelse(length(grids) == 1, 
                     one_procedure, 
                     sim_recursive)
  
  next_grid <- config[[grids[1]]]
  
  if (config$verbose) {
    cat("Entered grid", names(config[grids[1]]), "at", date(),
        "with length(prev_results) = ",
        length(prev_results),"\n")
  }
  
  results <- prev_results
  for( x in next_grid ) {
    
    configx <- config
    configx[[grids[1]]] <- x
    results <- next_sim(configx, results)
  }
  
  if (config$verbose) {
    cat("Unwound grid", 
        names(config[grids[1]]), "at", date(),
        "with length(prev_results) = ",
        length(prev_results), "\n")
  }
  
  return(results)
}
