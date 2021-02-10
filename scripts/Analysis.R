####
# It could be worth considering RcppArmadillo replacement for mclust
# https://github.com/hiendn/StoEMMIX/blob/master/R/EMalgorithmforGaussianMixtureARMA_exc.R
####

# Simulation --------------------------------------------------------------

source("Helpers.R")
packages <- c("mclust", 
              "MixSim", 
              "doParallel", 
              "parallel", 
              "foreach",
              "magrittr",
              "ggplot2",
              "dplyr",
              "tidyr",
              "Rcpp",
              "RcppArmadillo",
              "inline")
load_packages(packages)

### Choose config script:
### * test_config.R - debugging / exploration
### * main_config.R - for full simulation
source("Test_config.R") 
#source("Main_config.R")

# Print the recorded start_time from config
cat("Start time:",substr(capture.output(print(config$start_time)),6,29),"\n")

# Perform simulations
results <- NULL
t <- system.time({
  results <- simulation(config)
})

# Save results
save_results(results, config$slurm_task_id)

# Print diagnostics
cat("Number of ms_draws errored out: ", sum(errored(results)),"\n")
cat("User time = ", t[1],"s\n")
cat("---- structure of output ----\n")
peek(results)