##################################################
### Simulation 1: one-block structure
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

# 1 block of markers
# within the block, there are multiple smaller sets of positive and negative markers


r <- D/10  # number of markers associated to the response
codes1 <- matrix(rep(0, (D-1)*D), ncol=D)  # to define SBP -> construction of pivot coordinates
for(ico in 1:(D-1)){
  codes1[ico,] = c(rep(0, ico-1), 1, rep(-1, D-ico))
}  


# Compute the basis of a clr-plane, to use with isometric log-ratio or planar transform of a (dataset of) compositions.
V1 <- gsi.buildilrBase(t(codes1))  

mu1 = c(rep(0, D-1))


# A matrix Sig1
Sig1 <- matrix(0, nrow = D-1, ncol = D-1)

posl <- seq(2,0.5,length.out = 2*r)    # non-diagonal elements (remove "2" in the next step)
posl <- posl[-1]  # remove a diagonal element (10)

aux <- list()
ind <- 1
for(i in (2*r-1):1){
  aux[[ind]] <- posl[1:i]
  ind <- ind + 1
}

aux<- unlist(aux)
block1 <- matrix(0,2*r,2*r)
block1[lower.tri(block1)] <- aux
block1[upper.tri(block1)] <- t(block1)[upper.tri(t(block1))]
diag(block1) <- 2

signs.matrix <- as.matrix(rbind(matrix(rep(c(rep(1,7),rep(-1,3),rep(1,5),rep(-1,5)),7),byrow = T, nrow = 7),
                                matrix(rep(c(rep(-1,7),rep(1,3),rep(-1,5),rep(1,5)),3),byrow = T, nrow = 3),
                                matrix(rep(c(rep(1,7),rep(-1,3),rep(1,5),rep(-1,5)),5),byrow = T, nrow = 5),
                                matrix(rep(c(rep(-1,7),rep(1,3),rep(-1,5),rep(1,5)),5),byrow = T, nrow = 5)))

Sig1[1:20,1:20] <- block1*signs.matrix
diag(Sig1)[21:(D-1)] <- 1


# Check if Sig1 is symmetrical and positive definite
isSymmetric(Sig1)
is.positive.definite(Sig1) 


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

sol1 = foreach(i = 1:S, .combine = "comb", .multicombine=TRUE, .verbose = T, .init=list(list()), 
               .packages = c("compositions","pls","MASS"))%dopar%{
                 
                 set.seed(i,"L'Ecuyer-CMRG")  
                 # data simulaton
                 d <- data_simulation(n,D,mu = mu1, codes = codes1, V = V1, r = D/10, Sig=Sig1, typeSig = "Sig1")   
                 X <- d$X
                 y <- d$y
                 nbal <- D-1

                 # calculation
                 rmsep.vec1 <- c()
                 for(m in 1:nbal){
                   rmsep.vec1[m] <- PBs_PCA_CV(X,y,nbalances=m,R=R, K=K, n=n)
                 } 
                 
                 # store results 
                 list(rmsep.vec1)
               }
proc.time()-ptm 


### PLS PBs algorithm
ptm <- proc.time()

sol2 = foreach(i = 1:S, .combine = "comb", .multicombine=TRUE, .verbose = T, .init=list(list()), 
               .packages = c("compositions","pls","MASS"))%dopar%{
                 
                 set.seed(i,"L'Ecuyer-CMRG")  
                 # data simulaton
                 d <- data_simulation(n,D,mu = mu1, codes = codes1, V = V1, r = D/10, Sig=Sig1, typeSig = "Sig1")   
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

sol3 = foreach(i = 1:S, .combine = "comb", .multicombine=TRUE, .verbose = T, .init=list(list()), 
               .packages = c("compositions","pls","MASS"))%dopar%{
                 
                 set.seed(i,"L'Ecuyer-CMRG")  
                 # data simulaton
                 d <- data_simulation(n,D,mu = mu1, codes = codes1, V = V1, r = D/10, Sig=Sig1, typeSig = "Sig1")  
                 X <- d$X
                 y <- d$y
                 nload <- D

                 # calculation
                 rmsep.vec1 <- c()
                 for(m in 1:nload){
                   rmsep.vec1[m] <- PLS_CV(X,y,nload=m,R=R, K=K, n=n)    # pripadne nahradit zpet PLS_CV               
                 }
                 
                 # store results 
                 list(rmsep.vec1)
               }
proc.time()-ptm 



### Selbal
ptm <- proc.time()

sol4 = foreach(i = 1:S, .combine = "comb", .multicombine=TRUE, .verbose = T, .init=list(list()), 
               .packages = c("compositions","pls","MASS"))%dopar%{
                 
                 set.seed(i,"L'Ecuyer-CMRG")  
                 # data simulaton
                 d <- data_simulation(n,D,mu = mu1, codes = codes1, V = V1, r = D/10, Sig=Sig1, typeSig = "Sig1")   # one block
                 X <- d$X
                 y <- d$y
                 #nload <- D
                 

                 # calculation
                 rmsep.vec1 <- c()
                 rmsep.vec1 <- Selbal_CV(X,y,R=R,K=K,n=n)                 
                 
                 # store results 
                 list(rmsep.vec1)
               }
proc.time()-ptm 

stopCluster(cl)

### Colect the results
t1 <- ToTable(sol1,S=S,D=D)         # PCA PBs
t2 <- ToTable(sol2,S=S,D=D)         # PLS PBs
t3 <- ToTable(sol3,S=S,D=D)         # PLS
t4 <- ToTableSelbal(sol4,S=S,D=1)   # Selbal



### Ability to correctly identify the biomarkers -------------------------------
# Construct "ideal balances" for comparison 
# -> balances, that consist just of marker variables and nothing else (num, den according to cov matrix)
# -> Fig 6

# Run sim 100 times; just construct 1st PCA_PB, PLS_PB, selbal balance


bal1.PCA <- matrix(NA,nrow=S,ncol=D)
bal1.PLS <- matrix(NA,nrow=S,ncol=D)

set.seed(1554)
for(i in 1:S){
  d <- data_simulation(n,D,mu = mu1, codes = codes1, V = V1, r = D/10, Sig=Sig1, typeSig = "Sig1")  
  X <- d$X
  y <- d$y
  
  # PCA PBs
  modelA <- fBalChip(X)
  bal1.PCA[i,] <- modelA$bal[1,]
  
  # PLS PBs
  modelB <- fBalChip_PLS(X,y,ver="cov")
  bal1.PLS[i,] <- modelB$bal[1,]
  
  print(i)
}

bal1.selbal <- matrix(NA,nrow=S,ncol=D)
set.seed(1554)
tic()
for(i in 1:S){
  d <- data_simulation(n, D, mu = mu1, codes = codes1, V = V1, r = D/10, Sig=Sig1, typeSig = "Sig1")
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
  bal1.selbal[i,] <- logcontr
  print(i)
}
toc()

# Construct "ideal balance" for comparison -> balance, that consist just of marker variables and nothing else (num, den according to cov matrix)

r1 <- 12
s1 <- 8

bpos1 <- (1/r1)*sqrt((r1*s1)/(r1+s1))
bneg1 <- (-1/s1)*sqrt((r1*s1)/(r1+s1))

indpos1 <- c(1:7,11:15)   # indices of variables in the numerator
indneg1 <- c(8:10,16:20)   # indices of variables in the denominator

logcontr1 <- rep(0,D)   # balance-logcontrast
logcontr1[indpos1] <- bpos1
logcontr1[indneg1] <- bneg1



