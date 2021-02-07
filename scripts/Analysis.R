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
              "dplyr")
load_packages(packages)

### Choose config script:
### * test_config.R - debugging / exploration
### * main_config.R - for full simulation
source("test_config.R") 
#source("main_config.R")

# one_procedure is useful for any simulation with no grids (just repetitions)
t <- system.time({
  repeat_sim(config, one_procedure)
})

# sim_recursive suits simulations where the procedure is 
# much more computationally intensive than switching between grid points
t <- system.time({
  results <- repeat_sim(config, sim_recursive)
})
cat("User time = ", t[1],"s")

save_results(results)