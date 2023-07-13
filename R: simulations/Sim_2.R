##################################################
### Simulation 2: 4 same-sized blocks structure
##################################################

### Clean workspace
rm(list = ls())


# R version: 4.1.3 (2022-03-10) -- "One Push-Up"

library(compositions)  # to handle CoDa; version 2.0-4
library(pls)           # PLS; version 2.8-0
library(matrixcalc)    # is.positive.definite; version 1.0-6
library(MASS)          # mvrnorm(); version 7.3-55
library(selbal)        # version 0.1.0
library(tictoc)        # sim. time measurement

# For parallel computing
library(parallel)      # version 4.1.3
library(foreach)       # version 1.5.2
library(doSNOW)        # version 1.0.20




### Source functions -----------------------------------------------------------

# PCA-PBs functions by Martin-Fernandez et al. (2018): Advances in Principal Balances for Compositional Data
source('fBalChip.r')  # Main function
source('fBPMaxOrthNewChip.r')
source('fBalChipman.r')
source('fBPUpChi.r')

# our proposal: PLS-PBs

source('fBalChip_PLS.R')   # Main function 
source('fBPMaxOrthNewChip_PLS.R')
source('fBalChipman_PLS.R')   
source('fBPUpChi_PLS.r')         

# Auxiliary functions for simulations
source('Auxiliary_functions.R')

#####


### Data simulation ------------------------------------------------------------
n=250   # number of observatios  
D=100   # number of parts/components

# 4 same-sized blocks of markers -> 20 markers in each
# each block contains different values 
# pattern of signs in the covariance structure: 7-3-5-5
# signs of betas reflect the signs of the sets within blocks: if 7+, 3-, 5+, 5-, then 7 betas+, 3 betas- etc.
# the last marker block is not related to the response at all (betas = 0)

r <- D/10  # number of markers associated to the response
codes1 <- matrix(rep(0, (D-1)*D), ncol=D)  # to define SBP -> construction of pivot coordinates
for(ico in 1:(D-1)){
  codes1[ico,] = c(rep(0, ico-1), 1, rep(-1, D-ico))
}  


# Compute the basis of a clr-plane, to use with isometric log-ratio or planar transform of a (dataset of) compositions.
V1 <- gsi.buildilrBase(t(codes1))  

mu1 = c(rep(0, D-1))

bl.size <- 10   # all blocks will be of same size, but different values inside

Sig2 <- matrix(0, nrow = D-1, ncol = D-1)

### Blocks of markers:
# first block of markers
val1 <- seq(10,1,length.out = 2*bl.size)    # non-diagonal elements (remove "2" in the next step)
val1 <- val1[-1]  # remove a diagonal element (10)

oD1 <- list()
ind <- 1
for(i in (2*bl.size-1):1){
  oD1[[ind]] <- val1[1:i]
  ind <- ind + 1
}

oD1 <- unlist(oD1)
B1.2 <- matrix(0,2*bl.size,2*bl.size)
B1.2[lower.tri(B1.2)] <- oD1
B1.2[upper.tri(B1.2)] <- t(B1.2)[upper.tri(t(B1.2))]
diag(B1.2) <- 10

# second block of markers
val2 <- seq(2,0.05,length.out = 2*bl.size)
val2 <- val2[-1]  # remove diagonal element = 2

oD2 <- list()
ind <- 1
for(i in (2*bl.size-1):1){
  oD2[[ind]] <- val2[1:i]
  ind <- ind + 1
}

oD2 <- unlist(oD2)
B2.2 <- matrix(0,2*bl.size,2*bl.size)
B2.2[lower.tri(B2.2)] <- oD2
B2.2[upper.tri(B2.2)] <- t(B2.2)[upper.tri(t(B2.2))]
diag(B2.2) <- 2


# third block of markers
val3 <- seq(4,0.1,length.out = 2*bl.size)
val3 <- val3[-1]  # remove diagonal element = 4

oD3 <- list()
ind <- 1
for(i in (2*bl.size-1):1){
  oD3[[ind]] <- val3[1:i]
  ind <- ind + 1
}

oD3 <- unlist(oD3)
B3.2 <- matrix(0,2*bl.size,2*bl.size)
B3.2[lower.tri(B3.2)] <- oD3
B3.2[upper.tri(B3.2)] <- t(B3.2)[upper.tri(t(B3.2))]
diag(B3.2) <- 4

# fourth block of markers
val4 <- seq(6,1,length.out = 2*bl.size)
val4 <- val4[-1]   # remove diagonal element = 6

oD4 <- list()
ind <- 1
for(i in (2*bl.size-1):1){
  oD4[[ind]] <- val4[1:i]
  ind <- ind + 1
}

oD4 <- unlist(oD4)
B4.2 <- matrix(0,2*bl.size,2*bl.size)
B4.2[lower.tri(B4.2)] <- oD4
B4.2[upper.tri(B4.2)] <- t(B4.2)[upper.tri(t(B4.2))]
diag(B4.2) <- 6


# create a chess-like block structure (7-3-5-5):
# first create a matrix of signs and use it to "resign" Sig2
signs.matrix <- as.matrix(rbind(matrix(rep(c(rep(1,7),rep(-1,3),rep(1,5),rep(-1,5)),7),byrow = T, nrow = 7),
                                matrix(rep(c(rep(-1,7),rep(1,3),rep(-1,5),rep(1,5)),3),byrow = T, nrow = 3),
                                matrix(rep(c(rep(1,7),rep(-1,3),rep(1,5),rep(-1,5)),5),byrow = T, nrow = 5),
                                matrix(rep(c(rep(-1,7),rep(1,3),rep(-1,5),rep(1,5)),5),byrow = T, nrow = 5)))

# matrix Sig2A completed:
Sig2[1:(2*bl.size),1:(2*bl.size)] <- B1.2*signs.matrix
Sig2[21:(20+2*bl.size),21:(20+2*bl.size)] <- B2.2*signs.matrix
Sig2[41:(40+2*bl.size),41:(40+2*bl.size)] <- B3.2*signs.matrix
Sig2[61:(60+2*bl.size),61:(60+2*bl.size)] <- B4.2*signs.matrix

diag(Sig2)[81:(D-1)] <- 1    # diagonal elements outside the blocks are set to 1

# Check if Sig2 is symmetrical and positive definite
isSymmetric(Sig2)
is.positive.definite(Sig2) 


### Cross-validation -----------------------------------------------------------

S <- 100   # number of simulation runs
K <- 5     # number of CV folds
R <- 1     # number of CV runs

### Initial setting for parallel computing

cl <- makeCluster(detectCores()-2)
registerDoSNOW(cl)

comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}



# AMD EPYC 7232P 8-Core Processor
# Cores: 8
# Logical processors: 16
# Freq: 3.10 GHz

# Approximate times
# sol1: 13.1min
# sol2: 32.3min
# sol3: 1.85min
# sol4: 26.3min



### PCA PBs algorithm

ptm <- proc.time()

sol1a = foreach(i = 1:S, .combine = "comb", .multicombine=TRUE, .verbose = T, .init=list(list()), 
                .packages = c("compositions","pls","MASS"))%dopar%{
                  
                  set.seed(i,"L'Ecuyer-CMRG") 
                  # data simulaton
                  d <- data_simulation(n,D,mu = mu1, codes = codes1, V = V1, r = D/10, Sig=Sig2, typeSig = "Sig2")   
                  X <- d$X
                  y <- d$y
                  nbal <- D-1
                  

                  # calculation
                  rmsep.vec1 <- c()
                  for(m in 1:nbal){
                    rmsep.vec1[m] <- PBs_PCA_CV(X,y,nbalances=m,R=R, K=K, n=n) }  
                  
                  # store results 
                  list(rmsep.vec1)
                }
proc.time()-ptm



### PLS PBs algorithm

ptm <- proc.time()

sol2a = foreach(i = 1:S, .combine = "comb", .multicombine=TRUE, .verbose = T, .init=list(list()), 
                .packages = c("compositions","pls","MASS"))%dopar%{
                  
                  set.seed(i,"L'Ecuyer-CMRG") 
                  # data simulaton
                  d <- data_simulation(n,D,mu = mu1, codes = codes1, V = V1, r = D/10, Sig=Sig2, typeSig = "Sig2")
                  X <- d$X
                  y <- d$y
                  nbal <- D-1

                  
                  # calculation
                  rmsep1 <- c()
                  for(m in 1:nbal){
                    rmsep1[m] <- PBs_PLS_CV(X,y,nbalances=m,R=R, K=K, n=n)
                  }
                  
                  # store results 
                  list(rmsep1)
                }
proc.time()-ptm 



### Standard PLS

ptm <- proc.time()

sol3a = foreach(i = 1:S, .combine = "comb", .multicombine=TRUE, .verbose = T, .init=list(list()), 
                .packages = c("compositions","pls","MASS"))%dopar%{
                  
                  set.seed(i,"L'Ecuyer-CMRG") 
                  # data simulaton
                  d <- data_simulation(n,D,mu = mu1, codes = codes1, V = V1, r = D/10, Sig=Sig2, typeSig = "Sig2")
                  X <- d$X
                  y <- d$y
                  nload <- D
                  
                  # calculation
                  rmsep.vec1 <- c()
                  for(m in 1:nload){
                    rmsep.vec1[m] <- PLS_CV(X,y,nload=m, R=R, K=K, n=n)        
                  }
                  
                  # store results 
                  list(rmsep.vec1)
                }
proc.time()-ptm 



### Selbal

ptm <- proc.time()

sol4a = foreach(i = 1:S, .combine = "comb", .multicombine=TRUE, .verbose = T, .init=list(list()), 
                .packages = c("compositions","pls","MASS"))%dopar%{
                  
                  set.seed(i,"L'Ecuyer-CMRG") 
                  # data simulaton
                  d <- data_simulation(n,D,mu = mu1, codes = codes1, V = V1, r = D/10, Sig=Sig2, typeSig = "Sig2")
                  X <- d$X
                  y <- d$y

                  # calculation
                  rmsep.vec1 <- c()
                  rmsep.vec1 <- Selbal_CV(X,y,R=R,K=K,n=n)                 

                  # store results 
                  list(rmsep.vec1)
                }
proc.time()-ptm 

stopCluster(cl)

### Colect the results
t1a <- ToTable(sol1a,S=S,D=D)  # PCA PBs
t2a <- ToTable(sol2a,S=S,D=D)  # PLS PBs
t3a <- ToTable(sol3a,S=S,D=D)  # PLS
t4a <- ToTableSelbal(sol4a,S=S,D=1) # Selbal



### Ability to correctly identify the biomarkers -------------------------------
# Construct "ideal balances" for comparison 
# -> balances, that consist just of marker variables and nothing else (num, den according to cov matrix)
# -> Fig 6

# Run sim 100 times; just construct 1st PCA_PB, PLS_PB, selbal balance

bal2.PCA <- matrix(NA,nrow=S,ncol=D)
bal2.PLS <- matrix(NA,nrow=S,ncol=D)

set.seed(1554)
for(i in 1:S){
  d <- data_simulation(n,D,mu = mu1, codes = codes1, V = V1, r = D/10,Sig=Sig2,typeSig = "Sig2")  
  X <- d$X
  y <- d$y
  
  # PCA PBs
  modelA <- fBalChip(X)
  bal2.PCA[i,] <- modelA$bal[1,]
  
  # PLS PBs
  modelB <- fBalChip_PLS(X,y,ver="cov")
  bal2.PLS[i,] <- modelB$bal[1,]
  
  print(i)
}


bal2.selbal <- matrix(NA,nrow=S,ncol=D)
set.seed(1554)
tic()
for(i in 1:S){
  d <- data_simulation(n,D,mu = mu1, codes = codes1, V = V1, r = D/10,Sig=Sig2,typeSig = "Sig2")
  X <- d$X
  y <- as.numeric(d$y)
  model <- selbal::selbal(X,y,draw=F)
  bilance <- model$balance.values
  
  # to obtain a logcontrast
  num <- model$numerator  # components in numerator
  den <- model$denominator  # components in denominator
  
  r <- length(num)
  s <- length(den)
  
  bpos <- (1/r)*sqrt((r*s)/(r+s))
  bneg <- (-1/s)*sqrt((r*s)/(r+s))
  
  cols <- colnames(X)
  
  indpos <- match(num,cols)   # indices of variables in the numerator
  indneg <- match(den,cols)   # indices of variables in the denominator
  
  logcontr <- rep(0,ncol(X))   # balance-logcontrast
  logcontr[indpos] <- bpos
  logcontr[indneg] <- bneg
  
  # final logcontrast:
  bal2.selbal[i,] <- logcontr
  print(i)
}
toc()

# Construct "ideal balance" for comparison -> balance, that consist just of marker variables and nothing else (num, den according to cov matrix)

inds1 <- rep(c(rep(T,7),rep(F,3),rep(T,5),rep(F,5)),3) # last block NOT related to the response -> not in ideal balance
posl1 <- seq(1,60,1)
indpos2 <- posl1[inds1]
indneg2 <- posl1[inds1==F]

r2 <- length(indpos2)
s2 <- length(indneg2)

bpos2 <- (1/r2)*sqrt((r2*s2)/(r2+s2))
bneg2 <- (-1/s2)*sqrt((r2*s2)/(r2+s2))

logcontr2 <- rep(0,D)   # balance-logcontrast
logcontr2[indpos2] <- bpos2
logcontr2[indneg2] <- bneg2
logcontr2

















