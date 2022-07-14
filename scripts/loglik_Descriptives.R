library(dplyr)
source("Descriptive_functions.R")
source("Utility_functions.R")
packages <- c("magrittr",
              "ggplot2",
              "dplyr",
              "tidyr",
              "xtable")
load_packages(packages)

dir <- "./temp_results/" # set to the directory of the simulation results files

# Formatting --------------------------------------------------------------

# load the result files
files <- paste0(dir, c("all_results.rds", "loglik_all_results.rds"))
results_list <- sapply(files, readRDS)

full_results <- unlist(results_list[[1]], recursive=F) # NN <= 500
loglik_results <- unlist(results_list[[2]], recursive=F) # NN >= 1000

# EXAMPLE: look at the data structure
attr(full_results[[1]], 'Log_lik_full')
attr(loglik_results[[1]], 'Log_lik_full')
sapply(full_results, min_AIC) %>% table
sapply(full_results, min_BIC) %>% table
lapply(loglik_results, \(x){attr(x, 'Log_lik_full')})

# EXAMPLE: extract all the datasets
datasets <- c(extract_datasets(full_results),
              extract_datasets(loglik_results))
str(datasets)
length(datasets)


# compute AIC and BIC
res_IC <- rbind(combine_AIC_BIC_results(full_results),
                combine_AIC_BIC_results(loglik_results))
saveRDS(res_IC, file=paste0(dir,"formatted_IC_results.rds"))

# Produce result tables ---------------------------------------------------

res_IC <- readRDS(paste0(dir,"formatted_IC_results.rds"))

# Group and summarise
resg_IC <- group_by(res_IC, TrueG, BO, DD, NN, LL)
msg_IC <- summarize(resg_IC,
                    mean_comp_AIC = mean(minAIC),
                    mean_comp_BIC = mean(minBIC),
                    prop_corr_AIC = mean(minAIC == TrueG),
                    prop_corr_BIC = mean(minBIC == TrueG),)

# Look at result tables
msg_IC %>% filter(BO==0.01, TrueG == 5) %>% View
msg_IC %>% filter(BO==0.05, TrueG == 5)
msg_IC %>% filter(BO==0.1, TrueG == 5)

msg_IC %>% filter(BO==0.01, TrueG == 10) %>% View
msg_IC %>% filter(BO==0.05, TrueG == 10)
msg_IC %>% filter(BO==0.1, TrueG == 10)

# LaTeX tables
xtable(msg_IC %>% filter(BO==0.01, TrueG == 5))
xtable(msg_IC %>% filter(BO==0.05, TrueG == 5))
xtable(msg_IC %>% filter(BO==0.1, TrueG == 5))

xtable(msg_IC %>% filter(BO==0.01, TrueG == 10))
xtable(msg_IC %>% filter(BO==0.05, TrueG == 10))
xtable(msg_IC %>% filter(BO==0.1, TrueG == 10))

xtable(msg_IC)


