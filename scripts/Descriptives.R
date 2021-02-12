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

files <- paste0("../results/local/sim_rcpp/results1_",c(1,2,3,4,5,6),".rds")

# load the result files
results_list <- sapply(files, readRDS)

# Peek at the result list structure
peek(results_list)

#undebug(read_clean_and_combine_first_accepts)

res <- read_clean_and_combine_first_accepts(results_list)

head(res)

# Table of sample sizes in each parameter combination
table(select(res, !starts_with("p"), -msid))

# Check the alpha_p3 (significance) for each msid
resg <- group_by(res, msid, TrueG)
msg <- summarize(resg
                 ,alpha_p1 = mean(p1 > 5) 
                 ,alpha_p2 = mean(p2 > 5)
                 ,alpha_p3 = mean(p3 > 5)
                 ,n = n()
)
msg
plot(msg$n, msg$alpha_p1, main = "alpha_p3 by mixture draw sample size")
hist(msg$alpha_p3, main = "alpha_p3 for each mixture draw", breaks=200)

# obtain msids above a certain alpha_3 cutoff point
good_msids <- filter(msg, alpha_p3 < 0.2)$msid

# Check that there are no bad msids
check <- length(good_msids) == nrow(msg)
if (!check) warning("There are bad msids")

# uncomment the RHS to filter only good msids 
resf <- filter(res, msid %in% good_msids)

# Group by unique parameter combination for summary stats
resg <- group_by(resf, across(!starts_with(c("p","msid"))))
summarize(resg
         ,n = n() 
         ,mu_p1 = mean(p1) 
         ,mu_p2 = mean(p2) 
         ,mu_p3 = mean(p3)
         ,alpha_p1 = mean(p1 > 5) 
         ,alpha_p2 = mean(p2 > 5)
         ,alpha_p3 = mean(p3 > 5)
         ,ms_draws = length(unique(msid))
) %>% filter(LL == 2) %>% 
  select(!starts_with("mu"), -n, -TrueG, -ms_draws)

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