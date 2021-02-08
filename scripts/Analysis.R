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
              "tidyr")
load_packages(packages)

### Choose config script:
### * test_config.R - debugging / exploration
### * main_config.R - for full simulation
source("Test_config.R") 
#source("Main_config.R")

# Perform simulations
results <- NULL
t <- system.time({
  results <- simulation(config)
})

# Save results
save_results(results, config$slurm_task_id)

# Print diagnostics
cat("Number of ms_draws errored out: ", sum(errored(results)),"\n")
cat("User time = ", t[1],"s")
cat("---- structure of output ----\n")
str(results, max.level=1, vec.len=2, list.len=2)