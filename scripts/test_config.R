# The simulation is chosen based on whatever variables in the config
# are of length > 1.
# **** Note **** 
# If config is created with no grids 
# then a quick visualisation of the
# mixture is produced.
# **************

config <- create_config(
  # n_1 (n_2) or n/2 sample size      
  NN = 1500
  # dimension of mixture model
  ,DD = 2  
  # value of desired average overlap for MixSim function (can be NULL)
  ,BO = c(0.01,0.05)
  # value of desired maximum overlap for MixSim function (can be NULL)
  #,MO = c(0.001,0.01,0.05)
  # generative number of components
  ,TrueG = 5  
  # lower bound for mixing proportions
  ,PiLow = 0.1   
  # additional components in alternative
  ,LL = 1      
  # covariance matrix structure (FALSE = non-spherical, TRUE = spherical).
  ,sph = F
  # force stop procedure when GG = TrueG + Gextra
  ,Gextra = 5  
  # whether to use fixed MixSim params wherever possible
  ,fixed_paras = T
  # error bound for overlap computation (default).
  ,eps = 1e-03
  # if you experience trouble try setting parallel = FALSE
  ,parallel = T     
  # whether to stop on the first failure to reject
  ,stop_on_accept = T  
  # Whether to leave a free core for the OS
  ,free_core = T  
  # whether to save results
  ,save = T 
  # number of repetitions of the sim function for repeat_sim
  ,reps = 100  
  # A name string for saving
  ,name = ""
  # whether to print diagnostic messages
  ,verbose= T
)
