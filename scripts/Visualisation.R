results <- readRDS("../results/local/sim_recursive/results_9529835.rds")

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

# Look at NA frequency
length(res$p1)
sum(is.na(res$p1[res$BO == 0.001]))

# Subset the results
res <- res[!is.na(res$p1) & res$BO != 0.001 & 
             res$LL == 1 & res$TrueG == 5 & res$DD == 2,]

sum(res$p3[res$NN == 500] > 5)/length(res$p3[res$NN == 500])
sum(res$p3[res$NN == 1000] > 5)/length(res$p3[res$NN == 1000])
sum(res$p3[res$NN == 1500] > 5)/length(res$p3[res$NN == 1500])


res1 <- res[res$TrueG == 5 & res$BO == 0.01 & res$L == 1 & res$NN == 1500 & res$DD == 2,]
nrow(res1)
sum(res1$p3 > 5, na.rm = T)

# Plot 1 ------------------------------------------------------------------
#
# First accept vs N with error bars for 95% quantile intervals.
# Red horizontal line indicates TrueG (ideal first accept).
# Notes: - Power is proportion that land on TrueG.
#        - A general notion of power can include up to g-k for tolerance k.

tg <- 5 # Pick which TrueG to plot
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
  geom_hline(yintercept = tg, col = "red", 
             linetype="dashed", size = 1.1) +
  scale_y_continuous(breaks = 1:6) +
  theme_minimal()
dev.off()