library(mclust)
library(MixSim)
library(Rcpp)
library(RcppArmadillo)
library(inline)

source("Big_test_arma_3.R")

params <- read.csv("to-catch.csv")
params <- params[,-1]
nrow(params)


make_indices <- function(params) {
  
  # An 8 element slurm array for 16 grid points has 2 grid points
  # per array element.
  index_set <- split(1:nrow(params), rep(1:8, each=2))
  
  sat_id <- Sys.getenv('SLURM_ARRAY_TASK_ID')
  if (sat_id != "") {
    sat_id <- as.integer(sat_id)
    indices <- index_set[[sat_id]]
  } else {
    indices <- 1:nrow(params) 
  }
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
    results_i <- do.call(arma3, par)
    results_list[[i]] <- results_i
  })
  time_list[[i]] <- t
  filename <- paste0(paste0(names(par),"=",par,"-", collapse = ""),".rds")
  saveRDS(structure(results_i,t=t), filename)
}
saveRDS(results_list, paste0(
  "all_results",Sys.getenv('SLURM_ARRAY_TASK_ID'),".rds"))
