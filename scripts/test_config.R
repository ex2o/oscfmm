# ********************* MixSim parameters *********************
# New mixture parameters will be drawn using MixSim::MixSim 
# for each row of ms_grid

# Basic grid positions:
ms_grid <- expand.grid(
  
  # dimension of mixture model
  DD=2  
  
  # value of desired average overlap for MixSim function (can be NULL)
  ,BO=c(0.01,0.05)
  
  # value of desired maximum overlap for MixSim function (can be NULL)
  #,MO=c(0.001,0.01,0.05)
  
  # generative number of components
  ,TrueG=5 #c(5,10)  
  
  # error bound for overlap computation (default).
  ,eps=1e-03
  
)

# Parameters defined in terms of basic grid positions:
ms_grid <- mutate(ms_grid
                  
                  # lower bound for mixing proportions
                  ,PiLow = 1/(2*TrueG)
                  
)


# ***************** Data Simulation grid ********************
# A new data set will be drawn using MixSim::simdataset 
# for each row of ds_grid, and this will be repeated 
# for each row of ms_grid

ds_grid <- expand.grid(
  
  # n_1 (n_2) or n/2 sample size      
  NN=c(1000,5000,10000)
  
  # additional components in alternative
  ,LL=c(1,2) 
  
)


# ***************** Once-only configurations ****************
# Each config parameter is unchanged for the entire simulation

config <- create_config(
  
  # number of draws of random parameters from MixSim (NULL is automatic)
  ms_draws=NULL
  
  # number of draws of random datasets from simdataset
  ,ds_draws=100
  
  # whether to parallelise (at the level of ms_draws only)
  ,parallel=T  
  
  # whether this is a SLURM array job
  ,slurm_array=T
  
  # the name of the grid where slurm_array_col is located
  ,slurm_array_grid="ds_grid"
  
  # the column of slurm_array_grid that will be used for the slurm_array job
  ,slurm_array_col="NN"
  
  # max time (in sec) that this job is allowed to run for before wrapping up
  ,max_elapsed=100
  
  # MixSim parameter grid
  ,ms_grid=ms_grid
  
  # Dataset Simulation grid
  ,ds_grid=ds_grid
  
  # force stop procedure when GG = TrueG + Gextra
  ,Gextra=5  
  
  # whether to use fixed MixSim params wherever possible
  ,fixed_paras=T
  
  # whether to stop on the first failure to reject
  ,stop_on_accept=T  
  
  # Whether to leave a free core for the OS
  ,free_core=F  
  
  # whether to print diagnostic messages
  ,verbose=T
  
  # covariance matrix structure (FALSE = non-spherical, TRUE = spherical).
  ,sph=F
  
  # the approximate time that the script started running
  ,start_time=Sys.time()
  
)

rm(list = c("ms_grid", "ds_grid"))
