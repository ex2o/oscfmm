library(mclust)
library(MixSim)
library(Rcpp)
library(RcppArmadillo)
library(inline)


source("Big_test_arma_3.R")

# Seed
# set.seed(4200) 


# Create results directory --------------------------------------
mainDir <- "."
temp_results <- "temp_results"
dir.create(file.path(mainDir, temp_results), showWarnings = FALSE)

# A grid of simulation parameters
# params <- expand.grid(NN = c(300,500,1000,2000,5000,10000)
#             ,TrueG = c(5,10)
#             ,LL = c(1,2)
#             ,DD = c(2,4)
#             ,BO = c(0.01,0.05,0.1)
#             ,NumSim = 2
#             ,GGextra = 1)
params <- expand.grid(NN = c(300, 500)
            ,TrueG = c(5)
            ,LL = c(1,2)
            ,DD = c(2)
            ,BO = c(0.01)
            ,NumSim = 2
            ,GGextra = 1)
nrow(params)
head(params, n=16)

make_indices <- function(params) {
  
  sat_id <- Sys.getenv('SLURM_ARRAY_TASK_ID')
  if (sat_id != "") {
    # Dividing the parameters into sets of 16 so that NN, TrueG and LL are
    # balanced, means we have 96/16 = 6 runs, that each run should take
    # about 26 hours. Provide 30 hours for each run in SLURM.
    index_set <- split(1:nrow(params), rep(1:6, each=16))
    
    sat_id <- as.integer(sat_id)
    indices <- index_set[[sat_id]]
  } else {
    indices <- 1:nrow(params) 
  }
  return(indices)
}

indices <- make_indices(params)

params <- params[indices,]

params$GGmax <- params$TrueG + params$GGextra
params$GGextra <- NULL


# Run simulations ---------------------------------------------------------
results_list <- list()
time_list <- list()
npar <- nrow(params)
for (i in 1:npar) {
  par <- params[i,]
  t <- system.time({
    results_i <- do.call(arma3, par)
    results_list[[i]] <- results_i
  })
  time_list[[i]] <- t
  filename <- paste0(temp_results,"/",paste0(names(par),"=",par,"-", collapse = ""),".rds")
  saveRDS(structure(results_i,t=t), filename)
  cat(paste0("completed params[",i,",] in time t=",round(t[3],1),"s. Remaining params=",npar-i,"\n"))
}
saveRDS(results_list, paste0(temp_results, "/",
  "all_results",Sys.getenv('SLURM_ARRAY_TASK_ID'),".rds"))


# Run simulations for loglik only -----------------------------------------

source("AIC_BIC_arma_3.R")

# A grid of simulation parameters for loglik only
# params <- expand.grid(NN = c(1000,2000,5000,10000)
#             ,TrueG = c(5,10)
#             ,LL = c(1,2)
#             ,DD = c(2,4)
#             ,BO = c(0.01,0.05,0.1)
#             ,NumSim = 100
#             ,GGextra = 2)
params <- expand.grid(NN = c(1000, 2000)
                      ,TrueG = c(5)
                      ,LL = c(1,2)
                      ,DD = c(2)
                      ,BO = c(0.01)
                      ,NumSim = 2
                      ,GGextra = 1)
nrow(params)
head(params, n=16)

make_indices <- function(params) {
  
  sat_id <- Sys.getenv('SLURM_ARRAY_TASK_ID')
  if (sat_id != "") {
    # Dividing the parameters into sets of 16 so that NN, TrueG and LL are
    # balanced, means we have 96/16 = 6 runs, that each run should take
    # about 26 hours. Provide 30 hours for each run in SLURM.
    index_set <- split(1:nrow(params), rep(1:6, each=16))
    
    sat_id <- as.integer(sat_id)
    indices <- index_set[[sat_id]]
  } else {
    indices <- 1:nrow(params) 
  }
  return(indices)
}

indices <- make_indices(params)

params <- params[indices,]

params$GGmax <- params$TrueG + params$GGextra
params$GGextra <- NULL

results_list <- list()
time_list <- list()
npar <- nrow(params)
for (i in 1:npar) {
  par <- params[i,]
  t <- system.time({
    results_i <- do.call(loglik_arma3, par)
    results_list[[i]] <- results_i
  })
  time_list[[i]] <- t
  filename <- paste0(temp_results,"/","loglik_",paste0(names(par),"=",par,"-", collapse = ""),".rds")
  saveRDS(structure(results_i,t=t), filename)
  cat(paste0("completed params[",i,",] in time t=",round(t[3],1),"s. Remaining params=",npar-i,"\n"))
}
saveRDS(results_list, paste0(temp_results, "/",
                             "loglik_all_results",Sys.getenv('SLURM_ARRAY_TASK_ID'),".rds"))

# estimate_full_time <- function(t) {
#   # Df <- 2
#   # Bf <- 3
#   # Lf <- 2
#    NSf <- 50
#    #return(t[1]*Df*Bf*Lf*NSf)
#    return(t[1]*NSf)
# }
# 
# full_time <- sum(unlist(lapply(time_list, estimate_full_time))/3600)







