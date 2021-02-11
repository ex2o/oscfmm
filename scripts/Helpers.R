# Utilities ---------------------------------------------------------------

# Peek at a list object with regular structure

peek <- function(obj) {
  str(obj, max.level=1, vec.len=2, list.len=2)
}

load_packages <- Vectorize(function(package) {
  
  if (!require(package, character.only = T)) {
    install.packages(package)
    require(package, character.only = T)
    return(F)
  }
  return(T)
})

combiner <- function(a, ...) {
  c(a, ...)
}

worker <- function() {
  paste(Sys.info()[['nodename']], Sys.getpid(), sep='-')
}

create_config <- function(...) {
  
  config <- list(...)
  
  cat("---- Mixture draw Grid ----\n")
  print(config$ms_grid)
  cat("\n")
  cat("---- Dataset draw Grid ----\n")
  print(config$ds_grid)
  cat("\n")
  
  if (config$use_rcpp) {
    config$g_stepper <- one_g_step_rcpp
    cat("Planning to use Rcpp","\n")
  } else {
    config$g_stepper <- one_g_step
    cat("Not planning to use Rcpp","\n")
  }
  
  config <- prepare_model_names(config)
  
  config <- prepare_parallel(config)
  
  config <- prepare_for_slurm_array(config)
  
  config
}

prepare_parallel <- function(config) {
  
  if (config$parallel) {
    
    config <- prepare_no_cores(config)
    
    config <- prepare_cluster(config)
    
  } else {
    
    if (is.null(config$ms_draws)) {
      config$ms_draws <- 1
      cat("ms_draws=NULL, automatically set to",1,
          "since parallel=FALSE\n")
    }
  }
  
  config
}

prepare_no_cores <- function(config) {
  
  # Detect, modify and check number of cores
  config$no_cores <- parallel::detectCores()
  if (config$free_core) {
    config$no_cores <- config$no_cores - 1
  }
  if (is.null(config$ms_draws)) {
    config$ms_draws <- config$no_cores
    cat("ms_draws=NULL, automatically set to no_cores=",config$ms_draws,"\n")
  }
  if (config$ms_draws %% config$no_cores != 0) {
    warning(
      "ms_draws=", config$ms_draws,
      " is not a multiple of no_cores=", config$no_cores,
      ", consider setting increasing ms_draws for free",
      " or setting ms_draws=NULL to automatically",
      " and set equal to number of cores")
  }
  config
}

prepare_cluster <- function(config) {
  
  # Prepare foreach loop arguments
  config$foreach_args <- list(
    i = 1:config$ms_draws
    ,.inorder = FALSE
    ,.combine = combiner, .init = list()
    ,.export = c("one_mixture_grid"
                ,"one_procedure" 
                ,"one_g_step"
                ,"retry_on_error" 
                ,"mixsim" 
                ,"simdata"
                ,"worker"
                ,"is_there_time"
                ,"kmeans"
                ,"create_rcpp_functions"
                ,"one_g_step_rcpp"
                ,"make_GMM_arma"
                ,"cxxfunction"
                ,"signature")
  )
  
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
  
  config
}

# If this is a SLURM array job, then it is assumed that the shell script 
# contains #SBATCH --array=1-n where n=length(NN). This way, each job
# corresponds to a single value of NN.

prepare_for_slurm_array <- function(config) {
  
  if (config$slurm_array) {
    
    config$slurm_task_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
    
    if (!is.na(config$slurm_task_id)) {
      cat("slurm array job, task id:", config$slurm_task_id,"\n")
    } else {
      stop("slurm array job but no slurm task id was found",
           "check the slurm_array setting in config\n")
    }
    
    slg <- config$slurm_array_grid
    slv <- config$slurm_array_col
    all_v <- config[[slg]][[slv]]
    unique_v <- unique(all_v)
    if (config$slurm_task_id > length(unique_v)) {
      stop("Task id", config$slurm_task_id,
           "greater than length(unique_v)=",
           length(unique_v))
    }
    this_v <- unique_v[config$slurm_task_id]
    config[[slg]] <- config[[slg]][all_v == this_v,]

    cat("this task uses this_v=",this_v,"\n")
    
  } else {
    
    config$slurm_task_id <- ""
  }
  
  config
}

prepare_model_names <- function(config) {
  
  config$modelNames <- if (config$sph) {
    "VII"
    message("spherical model only")
  } else {
    "VVV"
  }
  config
}


# Assuming no more than 100 simultaneous jobs on SLURM
save_results <- function(results,slurm_task_id) {
  
  for (i in 1:1000) {
    file <- paste0("results",i,"_",slurm_task_id,".rds")
    if (!file.exists(file)) {
      saveRDS(results, file = file)
      return("Success")
    }
  }
  error("Unable to continue saving results")
}

# After max_errors, this function gives up and returns an error, unless
# no_error=TRUE, in which case it returns only a message string.
# if this function is used in a foreach loop then config should be passed
# in to allow the cluster to be stopped

# After max_errors, this function gives up and returns an error, unless
# no_error=TRUE, in which case it returns only a message string.
# if this function is used in a foreach loop then config should be passed
# in to allow the cluster to be stopped

# retry_on_error <- function(
#   expr,
#   is_error=function(x) "try-error" %in% class(x),
#   max_errors=10,
#   no_error=F,
#   prev_attempts=0,
#   config=NULL) {
# 
#   attempts <- prev_attempts
#   retval <- try(eval(expr))
#   while (is_error(retval)) {
#     attempts <- attempts + 1
#     if (attempts >= max_errors) {
#       msg <- sprintf("retry: too many retries [[%s]]\n",
#                      retval[1])
#       retval <- msg
#       if(!no_error) {
#         if (!is.null(config)) {
#           stopCluster(config$cl)
#         }
#         stop(msg)
#       }
#       break
#     } else {
#       msg <- sprintf("retry: error in attempt %i/%i [[%s]]\n",
#                      attempts,
#                      max_errors,
#                      retval[1])
#       if (!is.null(config)) {
#         if (config$verbose) {
#           cat(">>>>>>>>>> errored <<<<<<<<<< at", date(),"on",worker(),"\n")
#         }
#       }
#       warning(msg)
#     }
#     retval <- try(eval(expr))
#   }
#   return(structure(retval, attempts = attempts))
# }

# Uncomment below (and comment above) for debugging

retry_on_error <- function(
  expr,
  is_error=function(x) "try-error" %in% class(x),
  max_errors=10,
  no_error=F,
  prev_attempts=0,
  config=NULL) {
  
  attempts <- prev_attempts
  retval <- eval(expr)
  return(structure(retval, attempts = attempts))
}

is_there_time <- function(config) {
  
  elapsed <- difftime(Sys.time(), config$start_time, units = "secs")[[1]]
  if(elapsed > config$max_elapsed) {
    answer <- FALSE
  } else {
    answer <- TRUE
  }
  
  config$there_is_time <- structure(answer, elapsed=elapsed)
  config
}

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

# Sim functions -----------------------------------------------------------

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

# An rcpp alternative to one_g_step. Creating these in here is a work around
# for Rcpp parallel processing problems with exporting to cores in foreach.
# These functions rely on GMM_arma, and they were not found unless they were
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

# one_g_step_rcpp <- function(Data1, Data2, GG, config) {
# 
#   LL <- config$ds_grid_point$LL
# 
#   # Initialize a clustering using k-means
#   K_means_init <- kmeans(Data1$X, centers = GG,
#                          iter.max=100, nstar=100)
# 
#   # Use mclust to convert kmeans to parameters
#   param_init <- mclust::meVVV(Data1$X,
#                    z=mclust::unmap(K_means_init$cluster),
#                    control=mclust::emControl(itmax=1))
# 
#   # Run GMM using arma
#   MC1_null <- GMM_arma(t(Data1$X),
#                  param_init$parameters$pro,
#                  param_init$parameters$mean,
#                  param_init$parameters$variance$sigma,
#                  2000, GG)
# 
#   # Initialize a clustering using k-means
#   K_means_init <- kmeans(Data2$X, centers = GG,
#                          iter.max=100, nstar=100)
# 
#   # Use mclust to convert kmeans to parameters
#   param_init <- mclust::meVVV(Data2$X,
#                    z=mclust::unmap(K_means_init$cluster),
#                    control=mclust::emControl(itmax=1))
# 
#   # Run GMM using arma
#   MC2_null <- GMM_arma(t(Data2$X),
#                        param_init$parameters$pro,
#                        param_init$parameters$mean,
#                        param_init$parameters$variance$sigma,
#                        2000, GG)
# 
#   # Initialize a clustering using k-means
#   K_means_init <- kmeans(Data1$X, centers = GG+LL,
#                          iter.max=100, nstar=100)
# 
#   # Use mclust to convert kmeans to parameters
#   param_init <- mclust::meVVV(Data1$X,
#                    z=mclust::unmap(K_means_init$cluster),
#                    control=mclust::emControl(itmax=1))
# 
#   # Run GMM using arma
#   MC1_alt <- GMM_arma(t(Data1$X),
#                        param_init$parameters$pro,
#                        param_init$parameters$mean,
#                        param_init$parameters$variance$sigma,
#                        2000, GG+LL)
# 
#   # Initialize a clustering using k-means
#   K_means_init <- kmeans(Data2$X, centers = GG+LL,
#                          iter.max=100, nstar=100)
# 
#   # Use mclust to convert kmeans to parameters
#   param_init <- mclust::meVVV(Data2$X,
#                    z=mclust::unmap(K_means_init$cluster),
#                    control=mclust::emControl(itmax=1))
# 
#   # Run GMM using arma
#   MC2_alt <- GMM_arma(t(Data2$X),
#                       param_init$parameters$pro,
#                       param_init$parameters$mean,
#                       param_init$parameters$variance$sigma,
#                       2000, GG+LL)
# 
# 
#   # Calculate test statistics
#   # Calculate test statistics
#   # Use arma code to calculate likelihood at specified 
#   # parameters (i.e., run with zero iterations)
#   Log_lik_MC2_alt_on_D1 <- GMM_arma(t(Data1$X),
#                                     MC2_alt$proportions, 
#                                     MC2_alt$means,
#                                     MC2_alt$covariances,
#                                     0,GG+LL)$`log-likelihood` 
#   
#   # Compute likelihood
#   V1 <- exp(Log_lik_MC2_alt_on_D1 - MC1_null$`log-likelihood`) # Compute
#   
#   Log_lik_MC1_alt_on_D2 <- GMM_arma(t(Data2$X),
#                                     MC1_alt$proportions, 
#                                     MC1_alt$means,
#                                     MC1_alt$covariances,
#                                     0,GG+LL)$`log-likelihood`
#   
#   V2 <- exp(Log_lik_MC1_alt_on_D2 - MC2_null$`log-likelihood`)
#   
#   V_bar <- (V1+V2)/2
# 
#   # Calculate P-values
#   P1 <- min(1/V1,1)
#   P2 <- min(1/V2,1)
#   P_bar <- min(1/V_bar,1)
# 
#   return(c(P1, P2, P_bar))
# }

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

# Summary functions -------------------------------------------------------

# return logical vector indicating if each result errored out
errored <- function(results) {
  
  unlist(lapply(results, function(x){grepl("retry.*",x[1])}))
}

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
add_quantiles <- function(results, alpha=0.1) {
  
  res <- results %>% 
    group_by(across(!starts_with("p"))) %>% 
    mutate(p1_qL = quantile(p1, probs=alpha/2)
          ,p2_qL = quantile(p2, probs=alpha/2)
          ,p3_qL = quantile(p3, probs=alpha/2)
          ,p1_qU = quantile(p1, probs=1-alpha/2)
          ,p2_qU = quantile(p2, probs=1-alpha/2)
          ,p3_qU = quantile(p3, probs=1-alpha/2))
  res
}

# The function takes an input vector of paths to multiple
# result files in rds format, and combines their cleaned
# first_accepts tables

read_clean_and_combine_first_accepts <- function(files) { 
  
  # Remove errors and check error rates
  results <- lapply(results_list, function(x){
    errd <- errored(x)
    cat("Errored=",sum(errd)/length(x),"\n")
    x[!errd]
  })
  
  # combine first accepts into a tables
  res <- lapply(results, combine_first_accepts)
  
  # extract ms_draw ids and add to the tables
  msid <- lapply(results, ms_draw_ids)
  for (i in 1:length(res)) {
    res[[i]] <- data.frame(res[[i]], msid = msid[[i]])
    cat("Number of ms_draws=",
        length(unique(msid[[i]])),"\n")
  }
  
  # Remove NAs and look at NA rates
  res0 <- lapply(res, function(x){
    x0 <- drop_na(x)
    cat("NA ratio=",(1 - nrow(x0)/nrow(x)),"\n")
    x0
  })
  
  # Bind into a single data.frame
  out <- do.call(rbind, res0)
  rownames(out) <- NULL
  
  # Drop the column eps
  out[,!(names(out) %in% "eps")]
}



# Give a unique ID to each unique ms_draw in the results

ms_draw_ids <- function(results) {
  
  # Using just the first element of the Mu matrix as a unique id for an ms_draw
  as.integer((
    unlist(lapply(results, function(x){attr(x,"sim_para")$Mu[1,1]}))
  )*1000000000)
}





