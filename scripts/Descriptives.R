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

### Uncomment if you have just run the simulations yourself: -----------

files <- paste0(dir, c("all_results.rds", "loglik_all_results.rds"))

## Uncomment to load from temp result files instead 
# names <- list.files(dir)
# names <- names[names != "all_results.rds"]
# names <- names[names != "formatted_results.Rds"]
# files <- paste0(dir, names)

# load the result files
results_list <- sapply(files, readRDS)

# Peek at the result list structure
peek(results_list)

full_results <- unlist(results_list[[1]], recursive=F)
loglik_results <- unlist(results_list[[2]], recursive=F)

peek(full_results)
peek(loglik_results)

res <- combine_first_accepts(full_results)
head(res)
saveRDS(res, file = paste0(dir,"formatted_results.Rds"))


# Get AIC and BIC ---------------------------------------------------------

# Example vector of log likelihoods on rbind(Data1$X, Data2$X), for AIC and BIC
attr(full_results[[1]], 'Log_lik_full')
attr(loglik_results[[1]], 'Log_lik_full')

sapply(full_results, min_AIC) %>% table
sapply(full_results, min_BIC) %>% table
lapply(loglik_results, \(x){attr(x, 'Log_lik_full')})

res_IC <- rbind(combine_AIC_BIC_results(full_results),
                combine_AIC_BIC_results(loglik_results))

saveRDS(res_IC, file=paste0(dir,"formatted_IC_results.rds"))

# Checking formatted results ----------------------------------------------

## Load the formatted / simplified results  
res <- readRDS(paste0(dir,"formatted_results.Rds"))
res_IC <- readRDS(paste0(dir,"formatted_IC_results.rds"))

# Checking results --------------------------------------------------------

names(res)

res %>% filter(TrueG == 10, ) %>% .[["p3"]] %>% hist(main = "First DNR, TrueG = 10")

res %>% filter(TrueG == 10, DD == 2, BO == 0.01) %>% .[["p3"]] %>% hist(main = "First DNR, TrueG = 10, DD = 2, BO = 0.01")

res %>% filter(TrueG == 10, DD == 2, BO == 0.1) %>% .[["p3"]] %>% hist(main = "First DNR, TrueG = 10, DD = 2, BO = 0.1")

res %>% filter(TrueG == 5) %>% .[["p3"]] %>% hist(main = "First DNR, TrueG = 5")

res %>% filter(TrueG == 5, DD == 2, BO == 0.01) %>% .[["p3"]] %>% hist(main = "First DNR, TrueG = 5, DD = 2, BO = 0.01")

res %>% filter(TrueG == 5, DD == 2, BO == 0.1) %>% .[["p3"]] %>% hist(main = "First DNR, TrueG = 5, DD = 2, BO = 0.1")

# Look at significance levels
res$LL <- as.integer(res$LL)
res$TrueG <- as.integer(res$TrueG)
res$NN <- as.integer(res$NN)
res$DD <- as.integer(res$DD)

resg <- group_by(res, TrueG, BO, DD, NN, LL)
msg <- summarize(resg
                 ,alpha_p1 = mean(p1 > TrueG) 
                 ,alpha_p2 = mean(p2 > TrueG)
                 ,alpha_p3 = mean(p3 > TrueG)
                 ,n = n())
as.data.frame(msg)

# Check if there were any false positives
res %>% summarise(sum(p1 > TrueG), sum(p2 > TrueG), sum(p3 > TrueG))



# Format results for paper ------------------------------------------------

# Look at power by:
# tabulate the mean number of components, 
# and the proportion of times it selected the true number.

resg <- group_by(res, TrueG, BO, DD, NN, LL)
msg <- summarize(resg,
                 mean_comp1 = mean(p1),
                 mean_comp2 = mean(p2),
                 mean_comp3 = mean(p3),
                 prop_corr1 = mean(p1 == TrueG),
                 prop_corr2 = mean(p2 == TrueG),
                 prop_corr3 = mean(p3 == TrueG))
View(msg)

msg %>% filter(BO==0.01, TrueG == 5)
msg %>% filter(BO==0.05, TrueG == 5)
msg %>% filter(BO==0.1, TrueG == 5)

msg %>% filter(BO==0.01, TrueG == 10)
msg %>% filter(BO==0.05, TrueG == 10)
msg %>% filter(BO==0.1, TrueG == 10)

xtable(msg %>% filter(BO==0.01, TrueG == 5))
xtable(msg %>% filter(BO==0.05, TrueG == 5))
xtable(msg %>% filter(BO==0.1, TrueG == 5))

xtable(msg %>% filter(BO==0.01, TrueG == 10))
xtable(msg %>% filter(BO==0.05, TrueG == 10))
xtable(msg %>% filter(BO==0.1, TrueG == 10))

xtable(msg)

# AIC and BIC results
resg_IC <- group_by(res_IC, TrueG, BO, DD, NN, LL)
msg_IC <- summarize(resg_IC,
                    mean_comp_AIC = mean(minAIC),
                    mean_comp_BIC = mean(minBIC),
                    prop_corr_AIC = mean(minAIC == TrueG),
                    prop_corr_BIC = mean(minBIC == TrueG),)


msg_IC %>% filter(BO==0.01, TrueG == 5) %>% View
msg_IC %>% filter(BO==0.05, TrueG == 5)
msg_IC %>% filter(BO==0.1, TrueG == 5)

msg_IC %>% filter(BO==0.01, TrueG == 10) %>% View
msg_IC %>% filter(BO==0.05, TrueG == 10)
msg_IC %>% filter(BO==0.1, TrueG == 10)

xtable(msg_IC %>% filter(BO==0.01, TrueG == 5))
xtable(msg_IC %>% filter(BO==0.05, TrueG == 5))
xtable(msg_IC %>% filter(BO==0.1, TrueG == 5))

xtable(msg_IC %>% filter(BO==0.01, TrueG == 10))
xtable(msg_IC %>% filter(BO==0.05, TrueG == 10))
xtable(msg_IC %>% filter(BO==0.1, TrueG == 10))

xtable(msg_IC)

# DECIDING ON NEW DATA TO ADD --------------------------------------
# 
# # Table of sample sizes in each parameter combination
# t <- table(dplyr::select(res, !starts_with("p")))
# tdf <- string2numeric(as.data.frame(t, stringsAsFactors = F))
# tdf_low <- subset(tdf, Freq < 100)
# grid <- dplyr::select(tdf_low, !starts_with("p"), -Freq)
# grid$NumSim <- 100
# grid$GGextra <- 1
# grid
# write.csv(grid, "to-catch.csv")
