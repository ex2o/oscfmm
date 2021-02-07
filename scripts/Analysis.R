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

t <- system.time({
  results <- simulation(config)
})

cat("User time = ", t[1],"s")

save_results(results)