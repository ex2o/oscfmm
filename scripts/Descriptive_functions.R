grid_position <- function(result) {
  config <- attr(result, "config")
  grid_pos <- c(config$NN, config$DD, 
                config$BO, config$TrueG, 
                config$LL)
  names(grid_pos) <- c("NN","DD","BO","TrueG","LL")
  return(grid_pos)
}

first_accept  <- function(res) {
  
  if (any(res < -99, na.rm=T)) {return(NA)}
  apply(res, MARGIN=2, FUN=function(x){
    w <- which(x >= 0.05)
    if (length(w) == 0) {return(NA)}
    w[1]})
}

combine_first_accepts <- function(results_list) {
  gp <- lapply(results_list, grid_position) %>% 
    do.call(rbind, .)
  fa <- lapply(results_list, first_accept) %>% 
    do.call(rbind, .)
  colnames(fa) <- c("p1","p2","p3")
  gpfa <- cbind(fa,gp)
  return(as_tibble(gpfa))
}

min_AIC <- function(x) {
  DD <- attr(x, 'config')$DD
  GG <- attr(x, 'config')$GG
  logL <- attr(x, 'Log_lik_full')
  k <- (1 + 2*DD + choose(DD,2))*GG # Each component density has 1 + 2d + choose(d,2) parameters
  which.min(2*(k - logL))
}

min_BIC <- function(x) {
  DD <- attr(x, 'config')$DD
  GG <- attr(x, 'config')$GG
  logN <- log(2*attr(x, 'config')$NN)
  logL <- attr(x, 'Log_lik_full')
  k <- (1 + 2*DD + choose(DD,2))*GG
  which.min(k*logN - 2*logL)
}

combine_AIC_BIC_results <- function(results_list) {
  gp <- lapply(results_list, grid_position) %>% 
    do.call(rbind, .)
  minAIC <- sapply(results_list, min_AIC)
  minBIC <- sapply(results_list, min_BIC)
  as_tibble(cbind(gp, minAIC, minBIC))
}
