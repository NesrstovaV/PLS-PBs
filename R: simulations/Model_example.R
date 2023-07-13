##########################################
## Model example on how to obtain PLS PBs
##########################################

### Clean workspace
rm(list = ls())

# 1) Load necessary functions and libraries ------------------------------------
source('fBalChip_PLS.R')   # Main function 
source('fBPMaxOrthNewChip_PLS.R')
source('fBalChipman_PLS.R')   
source('fBPUpChi_PLS.r')      

library(MASS)
library(compositions)


# 2) Generate sample data ------------------------------------------------------
n <- 100              # observations
D <- 15               # parts/variables
Sig <- diag(D-1)      # positive-definite symmetric matrix -> covariance matrix
mu <- c(rep(0, D-1))  # means of variables

set.seed(123)
# ilr coordinates
Z <- mvrnorm(n,mu,Sigma = Sig) 

# Z -> CoDa X
V <- ilrBase(D = D)  # ilrBase() in library(compositions)

X <- as.matrix(as.data.frame(acomp(exp(Z%*%t(V)))))

# Response y:
beta <- runif(D-1,0.1,1)
eps <- rnorm(n)
y <- Z%*%beta+eps


# 3) Calculate PLS PBs -----------------------------------------------------------
PLS_balances <- fBalChip_PLS(X,y,version = "cov")     # version = "cov" -> max. covariance

balances <- PLS_balances$bal
balances  
  



