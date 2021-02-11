# Summary functions -------------------------------------------------------

# return logical vector indicating if each result errored out
errored <- function(results) {
  
  unlist(lapply(results, function(x){grepl("retry.*",x[1])}))
}

first_accept  <- function(res) {
  
  lapply(res, function(x){
    if (any(x < -99, na.rm=T)) {return(NA)}
    apply(x, MARGIN=2, FUN=function(x){which(x >= 0.05)[1]})
  })
}

# Reshape and combine results for first accept, in long format
combine_first_accepts <- function(results) {
  
  grid <- results %>% 
    lapply(function(x){attr(x,"grid_pos")}) %>% 
    do.call(rbind, .)
  res <- results %>% 
    first_accept %>% 
    do.call(rbind, .) %>% 
    data.frame() %>% 
    set_names(c("p1","p2","p3"))
  res <- cbind(res,grid)
  res
}

# Add a quantile interval for each grid point
add_quantiles <- function(results, alpha=0.1) {
  
  res <- results %>% 
    group_by(across(!starts_with("p"))) %>% 
    mutate(p1_qL = quantile(p1, probs=alpha/2)
           ,p2_qL = quantile(p2, probs=alpha/2)
           ,p3_qL = quantile(p3, probs=alpha/2)
           ,p1_qU = quantile(p1, probs=1-alpha/2)
           ,p2_qU = quantile(p2, probs=1-alpha/2)
           ,p3_qU = quantile(p3, probs=1-alpha/2))
  res
}

# The function takes an input vector of paths to multiple
# result files in rds format, and combines their cleaned
# first_accepts tables

read_clean_and_combine_first_accepts <- function(results_list) { 
  
  # Remove errors and check error rates
  results <- lapply(results_list, function(x){
    errd <- errored(x)
    cat("Errored=",sum(errd)/length(x),"\n")
    x[!errd]
  })
  
  # combine first accepts into a tables
  res <- lapply(results, combine_first_accepts)
  
  # extract ms_draw ids and add to the tables
  msid <- lapply(results, ms_draw_ids)
  for (i in 1:length(res)) {
    res[[i]] <- data.frame(res[[i]], msid = msid[[i]])
    cat("Number of ms_draws=",
        length(unique(msid[[i]])),"\n")
  }
  
  # Remove NAs and look at NA rates
  res0 <- lapply(res, function(x){
    x0 <- drop_na(x)
    cat("NA ratio=",(1 - nrow(x0)/nrow(x)),"\n")
    x0
  })
  
  # Bind into a single data.frame
  out <- do.call(rbind, res0)
  rownames(out) <- NULL
  
  # Drop the column eps
  out[,!(names(out) %in% "eps")]
}



# Give a unique ID to each unique ms_draw in the results

ms_draw_ids <- function(results) {
  
  # Using just the first element of the Mu matrix as a unique id for an ms_draw
  as.integer((
    unlist(lapply(results, function(x){attr(x,"sim_para")$Mu[1,1]}))
  )*1000000000)
}