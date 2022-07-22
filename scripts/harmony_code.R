# Load libraries
library(harmony)
library(mclust)
library(fields)

# Load and process data
data(cell_lines)
V <- cell_lines$scaled_pcs
meta_data <- cell_lines$meta_data

# Put data in XX and draw two samples of size n/2
XX <- V[,1:2]
sample_index <- sample(dim(XX)[1])

X1 <- XX[sample_index[1:(length(sample_index)/2)],]
X2 <- XX[sample_index[(length(sample_index)/2+1):(length(sample_index))],]


# Testing G = gg versus G = gg + ll
gg <- 4
ll <- 2
MC_gg <- Mclust(XX,G=gg,model='VVV')
MC_ggll <- Mclust(XX,G=gg+ll,model='VVV')

null_1 <- em(X1,'VVV',MC_gg$parameters)
null_2 <- em(X2,'VVV',MC_gg$parameters)

alter_1 <- em(X1,'VVV',MC_ggll$parameters)
alter_2 <- em(X2,'VVV',MC_ggll$parameters)

a2l1 <- estepVVV(X1,alter_2$parameters)
a1l2 <- estepVVV(X2,alter_1$parameters)

# log-pvalues
holdout_1 <- null_1$loglik- a2l1$loglik
holdout_2 <- null_2$loglik - a1l2$loglik 

# Compute p-values
exp(holdout_1)
exp(holdout_2)
2/(1/exp(holdout_1)+1/exp(holdout_2))

# Compute BIC and AIC
gg <- 8
big_MC <- Mclust(XX,G=gg,model='VVV')
(-2*big_MC$loglik + (gg*(1+2*2+choose(2,2)))*log(dim(XX)[1]))/dim(XX)[1]
(-2*big_MC$loglik + gg*(1+2*2+choose(2,2))*2)/dim(XX)[1]

# Plot data
plot(XX,col=tim.colors(5,alpha=0.6)[as.numeric(factor(paste0(meta_data$dataset, meta_data$cell_type)))],
     pch=as.numeric(factor(paste0(meta_data$dataset, meta_data$cell_type))),lwd=2)
