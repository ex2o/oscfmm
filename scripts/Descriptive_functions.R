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
