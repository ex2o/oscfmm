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

retry_on_error <- function(
  expr,
  is_error=function(x) "try-error" %in% class(x),
  max_errors=10,
  no_error=F,
  prev_attempts=0,
  config=NULL) {

  attempts <- prev_attempts
  retval <- try(eval(expr))
  while (is_error(retval)) {
    attempts <- attempts + 1
    if (attempts >= max_errors) {
      msg <- sprintf("retry: too many retries [[%s]]\n",
                     retval[1])
      retval <- msg
      if(!no_error) {
        if (!is.null(config)) {
          stopCluster(config$cl)
        }
        stop(msg)
      }
      break
    } else {
      msg <- sprintf("retry: error in attempt %i/%i [[%s]]\n",
                     attempts,
                     max_errors,
                     retval[1])
      if (!is.null(config)) {
        if (config$verbose) {
          cat(">>>>>>>>>> errored <<<<<<<<<< at", date(),"on",worker(),"\n")
        }
      }
      warning(msg)
    }
    retval <- try(eval(expr))
  }
  return(structure(retval, attempts = attempts))
}

# Uncomment below (and comment above) for debugging

# retry_on_error <- function(
#   expr,
#   is_error=function(x) "try-error" %in% class(x),
#   max_errors=10,
#   no_error=F,
#   prev_attempts=0,
#   config=NULL) {
#   
#   attempts <- prev_attempts
#   retval <- eval(expr)
#   return(structure(retval, attempts = attempts))
# }


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