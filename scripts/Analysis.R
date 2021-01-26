source("Helpers.R")
packages <- c("mclust", 
              "MixSim", 
              "doParallel", 
              "parallel", 
              "foreach")
load_packages(packages)

config <- create_config(
  sim = sim_NN     # the simulation function to repeat
 # ,NN = 2000        # n_1 (n_2) or n/2 sample size
 ,NN_range = c(seq(100,1000,by=100), seq(1500,3000,by=500), 4000, 5000)
 ,DD = 2           # dimension of mixture model
 ,BO = 0.2         # difference in mixtures for MixSim function
 ,TrueG = 4        # generative number of components
 ,LL = 1           # additional components in alternative
 ,Gextra = 5       # force stop procedure when GG = TrueG + Gextra
 ,parallel = T     # if you experience trouble try setting parallel = FALSE
 ,stop = T         # whether to stop on the first failure to reject
 ,free_core = T    # Whether to leave a free core for the OS
 ,save = T         # whether to save results
 ,reps = 1         # number of repetitions of the sim function for repeat_sim
 ,name = ""        # A name string for saving
)

results <- repeat_sim(config)

saveRDS(results, file = "results.rds")

lapply(results, first_reject)


