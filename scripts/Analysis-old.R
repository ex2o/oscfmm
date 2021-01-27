####
# It could be worth considering RcppArmadillo replacement for mclust
# https://github.com/hiendn/StoEMMIX/blob/master/R/EMalgorithmforGaussianMixtureARMA_exc.R
####

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

# The simulation is chosen based on whatever variable in the config
# is of length > 1 (provided the sim has been implemented).
config <- create_config(
 # ,NN = 2000        # n_1 (n_2) or n/2 sample size
  NN = 100 #c(100,200,300) #c(seq(100,1000,by=100), seq(1500,3000,by=500), 4000, 5000)
 ,DD = 2           # dimension of mixture model
 ,BO = 0.05 #c(0.05,0.2) # difference in mixtures for MixSim function
 ,TrueG = c(4,5)   # generative number of components
 ,PiLow = 0.1      # lower bound for mixing proportions
 ,LL = 1           # additional components in alternative
 ,Gextra = 5       # force stop procedure when GG = TrueG + Gextra
 ,parallel = F     # if you experience trouble try setting parallel = FALSE
 ,stop = T         # whether to stop on the first failure to reject
 ,free_core = T    # Whether to leave a free core for the OS
 ,save = T         # whether to save results
 ,reps = 1         # number of repetitions of the sim function for repeat_sim
 ,name = ""        # A name string for saving
)

results <- repeat_sim(config)
save_results(results)
#results <- readRDS("results.rds")

# Reshape and combine results in long format
fr <- results %>%
  lapply(first_reject) %>% 
  lapply(function(x){
    df <- data.frame(do.call(rbind,x))
    df$NN <- config$NN_range; df
  }) %>% 
  do.call(rbind, .) %>% 
  (function(x){
    x$run <- factor(rep(
      1:length(results),each=length(config$NN_range))) 
    x
  })


# Add a quantile interval to the plot for each NN
q <- fr %>% 
  group_by(NN) %>% 
  mutate(X1_qL = quantile(X1, probs = 0.025),
         X1_qU = quantile(X1, probs = 0.975),
         X2_qL = quantile(X2, probs = 0.025),
         X2_qU = quantile(X2, probs = 0.975),
         X3_qL = quantile(X3, probs = 0.025),
         X3_qU = quantile(X3, probs = 0.975),)

# Plot
pdf(file="sim_NN.pdf",width=5,height=4)
ggplot(data = fr, aes(x = NN, y = X3)) +
  #geom_point() +
  geom_jitter(col = "green") +
  geom_errorbar(data = q, aes(ymin = X3_qL, ymax = X3_qU), 
                col = "blue", size = 1.1) +
  geom_hline(yintercept = 4, col = "red", 
             linetype="dashed", size = 1.1) +
  scale_y_continuous(breaks = 1:6) +
  theme_minimal()
dev.off()