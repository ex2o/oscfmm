source("Descriptive_functions.R")
source("Utility_functions.R")
packages <- c("magrittr",
              "ggplot2",
              "dplyr",
              "tidyr")
load_packages(packages)

#results <- readRDS("../results/local/sim_recursive/results_9529835.rds")
#results <- readRDS("../results/local/sim_timed/results1_1.rds")

#results <- readRDS("results1_.rds")

dir <- "../../results/local/new_method/"
names <- list.files(dir)
files <- paste0(dir, names)

# load the result files
results_list <- sapply(files, readRDS)

# Peek at the result list structure
peek(results_list)

#undebug(read_clean_and_combine_first_accepts)

results_list <- remove_errors(results_list)
res <- clean_and_combine_first_accepts(results_list)

#write.csv(res, file = "preliminary_results.csv")
head(res)


# DECIDING ON WHICH NEW DATA TO ADD --------------------------------------

# Table of sample sizes in each parameter combination
t <- table(dplyr::select(res, !starts_with("p"), -msid)); t

# Get parameter combinations that need more data
# and export them for the HPC
tdf <- string2numeric(as.data.frame(t, stringsAsFactors = F))
tdf_low <- subset(tdf, (NN < 1e5 & Freq < 1000) | (NN == 1e5 & Freq < 200))
ms_grid <- dplyr::select(tdf_low, c(DD, BO, TrueG))
ds_grid <- dplyr::select(tdf_low, c(NN, LL))
ms_grid$eps <- 0.001
ms_grid
ds_grid
# saveRDS(ms_grid, "ms_grid.rds")
# saveRDS(ds_grid, "ds_grid.rds")

min_draws <- 16
min_test_runs <- 10
num_hypers <- 144

# Count the number of unique sets of hyperparameters
count_unique_hypers(res, target_num_hypers = num_hypers)

# Check number of test runs per mixture draw
count_test_runs_per_draw(res, min_test_runs = min_test_runs)

# Check number of mixture draws per hyperparameter set
count_draws_per_hyper(res, min_draws = min_draws)

# If we remove the mixture draws with less than min_test_runs test runs,
# then we need this many more test runs
discard_draws_and_compute_required_runs(
  res, min_test_runs = min_test_runs, min_draws = min_draws, 
  target_num_hypers = num_hypers)

# If, on the other hand, we add test runs to the existing mixture draws
# so that every mixture draw has at least min_test_runs, but only after
# restricting our attention to the top k=min_draws draws with the most test runs,
# then we only need to make this many more test runs:
restk <- group_by(res, NN, LL, DD, BO, TrueG, PiLow) %>% 
  filter(is_in_top_k(msid, k = min_draws))
compute_test_runs_to_add(restk, min_test_runs = min_test_runs)


# LOOKING AT SIGNIFICANCE LEVELS ------------------------------------------

# Check the alpha_p3 (significance) for each msid
resg <- group_by(res, msid, TrueG)
msg <- summarize(resg
                 ,alpha_p1 = mean(p1 > TrueG) 
                 ,alpha_p2 = mean(p2 > TrueG)
                 ,alpha_p3 = mean(p3 > TrueG)
                 ,n = n()
)
msg
plot(msg$n, msg$alpha_p1, main = "alpha_p3 by mixture draw sample size")
hist(msg$alpha_p3, main = "alpha_p3 for each mixture draw", breaks=200)


# obtain msids above a certain alpha_3 cutoff point
good_msids <- filter(msg, alpha_p3 < 0.2)$msid

# Check that there are no bad msids
check <- length(good_msids) == nrow(msg)
if (!check) warning("There are ", nrow(msg)-length(good_msids)," bad msids ",
                    "ratio = ", 1-length(good_msids)/ nrow(msg))

# uncomment the RHS to filter only good msids 
resf <- filter(res, msid %in% good_msids)

# Group by unique parameter combination for summary stats
resg <- group_by(resf, NN, LL, DD, BO, TrueG, PiLow)
resgs <- summarize(resg
         ,n = n() 
         ,mu_p1 = mean(p1) 
         ,mu_p2 = mean(p2) 
         ,mu_p3 = mean(p3)
         ,alpha_p1 = mean(p1 > TrueG) 
         ,alpha_p2 = mean(p2 > TrueG)
         ,alpha_p3 = mean(p3 > TrueG)
         ,ms_draws = length(unique(msid))
) 
cat("There are",nrow(resgs),"hyperparameter sets")

#%>% filter(LL == 2) %>% 
#  select(!starts_with("mu"), -n, -TrueG, -ms_draws)

cat("These are the alphas over 0.05:")
resgs$alpha_p3[resgs$alpha_p3 > 0.05]

# Subset the results
res1 <- resf %>% filter(
  BO == 0.01
 ,LL == 1
 ,TrueG == 5
 ,DD == 2
 ,msid %in% good_msids
)

# Look at subset sample size
nrow(res1)

# Look at subset size ratio
nrow(res1)/nrow(resf)

# EXPLORING THE BAD MSIDS -------------------------------------------------

bad_msids <- (group_by(res, msid, TrueG) %>% 
  summarize(alpha_p3 = mean(p3 > TrueG)) %>% 
  filter(alpha_p3 >= 0.5))$msid
length(bad_msids)

# extract the mixture parameters from the results list if they belong
# to the bad_msids.

all_sim_paras <- unlist(lapply(results_list, function(xl) {
  lapply(xl, function(x){
    attr(x,"sim_para")  
  })
}), recursive=F)

all_msids <- unlist(lapply(all_sim_paras, make_msid), use.names=F)

is_bsp <- all_msids %in% bad_msids
bad_sim_paras <- unique(all_sim_paras[is_bsp])
good_sim_paras <- unique(all_sim_paras[!is_bsp])

length(bad_sim_paras) / (length(bad_sim_paras) + length(good_sim_paras))

get_dets <- function(sim_para) {
  ns <- dim(sim_para$S)[3]
  dets <- vector(mode="numeric",length=ns)
  for( i in 1:ns ) {
    S <- sim_para$S[,,i]
    dets[i] <- determinant(S)$modulus
  }
  dets
}

bad_dets <- unlist(lapply(bad_sim_paras, get_dets))
good_dets <- unlist(lapply(good_sim_paras, get_dets))

plot(good_dets)
points(bad_dets, col="red")

hist(good_dets)
hist(bad_dets)

determinant(spb$S[,,1])
MixSim::pdplot(spb$Pi, spb$Mu, spb$S)

spg <- good_sim_paras[[1]]
spg$S
MixSim::pdplot(spg$Pi, spg$Mu, spg$S)


# Plot 1 ------------------------------------------------------------------
#
# First accept vs N with error bars for 95% quantile intervals.
# Red horizontal line indicates TrueG (ideal first accept).
# Notes: - Power is proportion that land on TrueG.
#        - A general notion of power can include up to g-k for tolerance k.

# Add a quantile interval for each grid point
resq <- add_quantiles(res1)

pconfig <- list(
   point_col      = "green"
  ,point_size     = 1.2
  ,jitter_height  = 0.1
  ,jitter_width   = 0.1
  ,errorbar_size  = 1.1
  ,errorbar_col   = "blue"
  ,hline_size     = 1.1
  ,hline_col      = "red" 
  ,hline_linetype = "dashed"
  ,y_scale_breaks = 1:6
)

pdf(file="first_accept_vs_N.pdf",width=5,height=4)
ggplot(data = resq, aes(x = NN, y = p3)) +
  geom_jitter(col = pconfig$point_col, 
              height = pconfig$jitter_height, 
              width = pconfig$jitter_width,
              size = pconfig$point_size) +
  geom_errorbar(aes(ymin = p3_qL, ymax = p3_qU), 
                col = pconfig$errorbar_col, 
                size = pconfig$errorbar_size) +
  geom_hline(yintercept = unique(resq$TrueG), 
             col = pconfig$hline_col, 
             linetype=pconfig$hline_linetype, 
             size = pconfig$hline_size) +
  scale_y_continuous(breaks = pconfig$y_scale_breaks) +
  theme_minimal()
dev.off()