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

# The simulation is chosen based on whatever variables in the config
# are of length > 1.
## **** Note **** 
# If config is created with no grids 
# then a quick visualisation of the
# mixture is produced.
## **************
config <- create_config(
  # n_1 (n_2) or n/2 sample size      
  NN = 1500 #seq(500,1500,by=500)
  # dimension of mixture model
  ,DD = c(2,4)   
  # value of desired average overlap for MixSim function (can be NULL)
  ,BO = c(0.001,0.01,0.05)
  # value of desired maximum overlap for MixSim function (can be NULL)
  #,MO = c(0.001,0.01,0.05)
  # generative number of components
  ,TrueG = c(5,10)  
  # lower bound for mixing proportions
  ,PiLow = 0.1   
  # additional components in alternative
  ,LL = c(1,2)      
  # covariance matrix structure (FALSE = non-spherical, TRUE = spherical).
  ,sph = F
  # force stop procedure when GG = TrueG + Gextra
  ,Gextra = 5  
  # error bound for overlap computation (default).
  ,eps = 1e-03
  # if you experience trouble try setting parallel = FALSE
  ,parallel = F     
  # whether to stop on the first failure to reject
  ,stop_on_accept = T  
  # Whether to leave a free core for the OS
  ,free_core = T  
  # whether to save results
  ,save = T 
  # number of repetitions of the sim function for repeat_sim
  ,reps = 2  
  # A name string for saving
  ,name = ""
  # whether to print diagnostic messages
  ,verbose= T
)

#config <- test_config

# sim_recursive suits simulations where the procedure is 
# much more computationally intensive than switching between grid points
t <- system.time({
  results <- repeat_sim(config, sim_recursive)
})
cat("User time = ", t[1],"s")

save_results(results)

results <- readRDS("results1.rds")

# Reshape and combine results for first accept, in long format
grid <- results %>% 
  lapply(function(x){attr(x,"grid_pos")}) %>% 
  do.call(rbind, .)
res <- results %>% 
  first_accept %>% 
  do.call(rbind, .) %>% 
  data.frame() %>% 
  set_names(c("p1","p2","p3"))
res <- cbind(res,grid)


# Plot 1 ------------------------------------------------------------------
#
# First accept vs N with error bars for 95% quantile intervals.
# Red horizontal line indicates TrueG (ideal first accept).
# Notes: - Power is proportion that land on TrueG.
#        - A general notion of power can include up to g-k for tolerance k.

tg <- 4 # Pick which TrueG to plot
res_tg <- filter(res, TrueG == tg)

# Add a quantile interval to the plot for each NN
res_tg_q <- res_tg %>% 
  group_by(NN) %>% 
  mutate(p1_qL = quantile(p1, probs = 0.025)
        ,p1_qU = quantile(p1, probs = 0.975)
        ,p2_qL = quantile(p2, probs = 0.025)
        ,p2_qU = quantile(p2, probs = 0.975)
        ,p3_qL = quantile(p3, probs = 0.025)
        ,p3_qU = quantile(p3, probs = 0.975))

pdf(file="sim_NN.pdf",width=5,height=4)
ggplot(data = res_tg_q, aes(x = NN, y = p3)) +
  #geom_point() +
  geom_jitter(col = "green", height = 0.1, width = 5) +
  geom_errorbar(aes(ymin = p3_qL, ymax = p3_qU), 
                col = "blue", size = 1.1) +
  geom_hline(yintercept = 4, col = "red", 
             linetype="dashed", size = 1.1) +
  scale_y_continuous(breaks = 1:6) +
  theme_minimal()
dev.off()


# Plot 2 ------------------------------------------------------------------




# Test config -------------------------------------------------------------

test_config <- create_config(
  # n_1 (n_2) or n/2 sample size      
  NN = c(500,500) #c(seq(500, 1000, by=100), seq(1500,3000,by=500), 4000, 5000)
  # dimension of mixture model
  ,DD = 2 #c(2,4)   
  # value of desired average overlap for MixSim function (can be NULL)
  ,BO = 0.001 #c(0.001,0.01,0.05)
  # value of desired maximum overlap for MixSim function (can be NULL)
  #,MO = 0.15
  # generative number of components
  ,TrueG = 5 #c(5,10)  
  # lower bound for mixing proportions
  ,PiLow = 0.1   
  # additional components in alternative
  ,LL = 1 #c(1,2)      
  # force stop procedure when GG = TrueG + Gextra
  ,Gextra = 5  
  # error bound for overlap computation (default).
  ,eps = 1e-03
  # if you experience trouble try setting parallel = FALSE
  ,parallel = F     
  # whether to stop on the first failure to reject
  ,stop_on_accept = T  
  # Whether to leave a free core for the OS
  ,free_core = T  
  # whether to save results
  ,save = T 
  # number of repetitions of the sim function for repeat_sim
  ,reps = 3  
  # A name string for saving
  ,name = ""
  # whether to print diagnostic messages
  ,verbose= T
)
