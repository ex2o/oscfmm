results <- readRDS("../results/local/sim_recursive/results_9529835.rds")


res <- combine_first_accepts(results)

# Look at NA frequency
length(res$p1)
sum(is.na(res$p1[res$BO == 0.001]))

# Remove errors and look at error rate
res0 <- res %>% drop_na()
(1 - nrow(res0)/nrow(res))

# Subset the results
res1 <- res0 %>% filter(
  BO != 0.001
  ,LL == 1
  ,TrueG == 5
  ,DD == 2
  ,NN == 500)

# Look at fail rate
sum(res1$p3 > 5)/nrow(res0)

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