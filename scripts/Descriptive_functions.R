# Summary functions -------------------------------------------------------

# logical vector of whether each msid in a list of non-unique msids 
# appears in the top k number of times in that list
is_in_top_k <- function(msids, k) {
  msidt <- table(msids)
  best_msids <- as.integer(names(sort(msidt, dec = T)[1:k]))
  msids %in% best_msids
}

compute_test_runs_to_add <- function(res, min_test_runs = 10) {
  resr <- count_test_runs_per_draw(res, min_test_runs = min_test_runs)
  resr_low <- resr$n[resr$n < min_test_runs]
  return(length(resr_low)*min_test_runs - sum(resr_low))
}

discard_draws_and_compute_required_runs <- function(
  res, min_test_runs = 10, 
  min_draws = 16, 
  target_num_hypers = 144) {
  
  resd <- discard_draws_with_few_test_runs(res, few = min_test_runs)
  # Check number of mixture draws per hyperparameter set
  rdh <- count_draws_per_hyper(resd, min_draws = min_draws)
  rdh_low <- rdh[rdh$n < min_draws,]
  target <- nrow(rdh_low)*min_draws 
  target <- target + (target_num_hypers-count_unique_hypers(resd))*min_draws
  actual <- sum(rdh_low$n)
  message("REQUIRE: ", target-actual, 
      " new mixture draws with ",min_test_runs,
      " test runs each ")
  return((target-actual)*min_test_runs)
}

count_unique_hypers <- function(res, target_num_hypers = 144) {
   ruh <- group_by(res, NN, LL, DD, BO, TrueG, PiLow) %>% 
     summarise()
   n_hypers <- nrow(ruh)
   n_missing_hypers <- target_num_hypers - n_hypers
   if (n_missing_hypers > 0) {
     message("There are ", n_missing_hypers, " missing hypers") 
   }
   return(n_hypers)
}

# Check number of test runs per mixture draw
count_test_runs_per_draw <- function(res, min_test_runs = 10) {
  rpd <-  group_by(res, msid) %>% 
    summarise(n=n())
  message("ratio of draws with < ",min_test_runs," test runs = ", 
          round(nrow(rpd[rpd$n < min_test_runs,])/nrow(rpd),2))
  return(return(rpd))
}

# Check number of mixture draws per hyperparameter set
count_draws_per_hyper <- function(res, min_draws = 16) {
  dph <- group_by(res, NN, LL, DD,BO, TrueG, PiLow) %>% 
    summarise(n=length(unique(msid)))
  message("ratio of hypers with < ", min_draws," draws = ", 
          round(nrow(dph[dph$n < min_draws,])/nrow(dph),2))
  return(dph)
}

# Remove the mixture draws with less than 10 test runs
discard_draws_with_few_test_runs <- function(res, few = 10) {
  rpd <- count_test_runs_per_draw(res, min_test_runs = few)
  bad_msids <- rpd[rpd$n < few,]$msid
  resd <- res[!(res$msid %in% bad_msids),]
  count_unique_hypers(resd)
  return(resd)
}



# return logical vector indicating if each result errored out
errored <- function(results) {
  
  unlist(lapply(results, function(x){grepl("retry.*",x[1])}))
}

first_accept  <- function(res) {
  
  lapply(res, function(x){
    if (any(x < -99, na.rm=T)) {return(NA)}
    apply(x, MARGIN=2, FUN=function(x){
      w <- which(x >= 0.05)
      if (length(w) == 0) {return(NA)}
      w[1]})
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

remove_errors <- function(results_list) {
  # Remove errors and check error rates
  results <- lapply(results_list, function(x){
    errd <- errored(x)
    cat("Errored=",sum(errd)/length(x),"\n")
    x[!errd]
  })
  results
}



# The function takes an input vector of paths to multiple
# result files in rds format, and combines their cleaned
# first_accepts tables

clean_and_combine_first_accepts <- function(results_list) { 
  
  # combine first accepts into a tables
  res <- lapply(results_list, combine_first_accepts)
  
  # extract ms_draw ids and add to the tables
  msid <- lapply(results_list, ms_draw_ids)
   
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
  unlist(lapply(results, function(x){make_msid(attr(x,"sim_para"))}))
}

# ms_draw_ids <- function(results) {
# 
#   # Using just the first element of the Mu matrix as a unique id for an ms_draw
#   as.integer((
#     unlist(lapply(results, function(x){attr(x,"sim_para")$Mu[1,1]}))
#   )*1000000000)
# }

make_msid <- function(sim_para) {
  as.integer(sim_para$Mu[1,1]*1000000000)
}
