source("Descriptive_functions.R")
source("Utility_functions.R")
packages <- c("magrittr",
              "ggplot2",
              "dplyr",
              "tidyr")
load_packages(packages)


dir <- "../results/local/new_method/"
names <- list.files(dir)
files <- paste0(dir, names)

# load the result files
results_list <- sapply(files, readRDS)

# Peek at the result list structure
peek(results_list)
str(results_list[[1]])

res <- combine_first_accepts(results_list)
head(res)


# DECIDING ON NEW DATA TO ADD --------------------------------------

# Table of sample sizes in each parameter combination
t <- table(dplyr::select(res, !starts_with("p")))
tdf <- string2numeric(as.data.frame(t, stringsAsFactors = F))
tdf_low <- subset(tdf, Freq < 100)
grid <- dplyr::select(tdf_low, !starts_with("p"), -Freq)
grid$NumSim <- 100
grid$GGextra <- 1
grid
write.csv(grid, "to-catch.csv")
