source("Helpers.R")
packages <- c("mclust", 
              "MixSim", 
              "doParallel", 
              "parallel", 
              "foreach")
load_packages(packages)

NN <- 2000 # n_1 (n_2) or n/2 sample size
DD <- 2 # Dimension of mixture model
BO <- 0.2 # Difference in mixtures for MixSim function
TrueG <- 4 # Generative number of components
LL <- 1 # additional components in alternative
Gmax <- TrueG + 1 # max number of null components in procedure
parallel <- T # If you experience trouble try setting parallel = FALSE
stop <- T # Whether to stop on the first failure to reject
free_core <- T # Whether to leave a free core for the OS
save <- T  # whether to save results
name <- "" # A name string for saving
reps <- 2 # Number of repetitions of the sim function for repeat_sim
sim <- sim_NN # The simulation function to repeat
params <- list( # parameters to pass to sim
  NN_range = c(seq(100,1000,by=100), seq(1500,3000,by=500), 4000, 5000)
 ,DD = DD
 ,BO = BO
 ,TrueG = TrueG
 ,LL = LL
 ,Gmax = Gmax
 ,stop = stop
 ,free_core = free_core)


results <- repeat_sim(sim, params, reps, save, name, parallel)

saveRDS(results, file = "results.rds")

lapply(results, first_reject)


