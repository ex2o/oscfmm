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
# If config is created with no grids then a quick visualisation of the
# mixture is produced.
config <- create_config(
  # n_1 (n_2) or n/2 sample size      
  NN = c(100,200,300) #c(seq(100,1000,by=100), seq(1500,3000,by=500), 4000, 5000)
  # dimension of mixture model
  ,DD = 2   
  # difference in mixtures for MixSim function
  ,BO = 0.05 #seq(0.025, 0.1, by=0.025)
  # generative number of components
  ,TrueG = 4:5 #seq(2, 8, by=2)     
  # lower bound for mixing proportions
  ,PiLow = 0.1   
  # additional components in alternative
  ,LL = 1      
  # force stop procedure when GG = TrueG + Gextra
  ,Gextra = 5   
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

results <- repeat_sim(config, sim_recursive)

save_results(results)
results <- readRDS("results1.rds")

# Reshape and combine results for first accept, in long format
grid <- results %>% 
  lapply(function(x){attr(x,"config")}) %>% 
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


