#####################################################
### Simulation 3: 4 different-sized blocks structure
#####################################################

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



# blocks of 30 and 10 markers
# all blocks have same range of values
# signs of betas again reflect the covariance structure
# last block of markers not related to the response

r <- D/10  # number of markers associated to the response
codes1 <- matrix(rep(0, (D-1)*D), ncol=D)  # to define SBP -> construction of pivot coordinates
for(ico in 1:(D-1)){
  codes1[ico,] = c(rep(0, ico-1), 1, rep(-1, D-ico))
}  


# Compute the basis of a clr-plane, to use with isometric log-ratio or planar transform of a (dataset of) compositions.
V1 <- gsi.buildilrBase(t(codes1))  

mu1 = c(rep(0, D-1))

mark1 <- 15  # half of a first block of markers (first block =  30 markers)
mark2 <- 5   # half of a second block of markers (second block = 10 markers)
mark3 <- 15  # analogous
mark4 <- 5   # analogous

Sig3 <- matrix(0, nrow = D-1, ncol = D-1)


### Blocks of markers generated separately:
# first block of markers
h1 <- seq(10,1,length.out = 2*mark1)    # non-diagonal elements (remove "2" in the next step)
h1 <- h1[-1]  # remove 2 -> it will be a diagonal element

oDi1 <- list()
ind <- 1
for(i in (2*mark1-1):1){
  oDi1[[ind]] <- h1[1:i]
  ind <- ind + 1
}

oDi1 <- unlist(oDi1)
BL1 <- matrix(0,2*mark1,2*mark1)
BL1[lower.tri(BL1)] <- oDi1
BL1[upper.tri(BL1)] <- t(BL1)[upper.tri(t(BL1))]
diag(BL1) <- 10

BL3 <- BL1   # first and third blocks are identical

# second block of markers
h2 <- seq(10,1,length.out = 2*mark2)
h2 <- h2[-1]

oDi2 <- list()
ind <- 1
for(i in (2*mark2-1):1){
  oDi2[[ind]] <- h2[1:i]
  ind <- ind + 1
}

oDi2 <- unlist(oDi2)
BL2 <- matrix(0,2*mark2,2*mark2)
BL2[lower.tri(BL2)] <- oDi2
BL2[upper.tri(BL2)] <- t(BL2)[upper.tri(t(BL2))]
diag(BL2) <- 10

BL4 <- BL2    # second and fourth blocks are identical

# create a chess-like block structure (7-3-5-5-5-5) and (7-3):
# first create a matrix of signs and use it to "resign" Sig3; sign.matrix2 is for larger blocks
signs.matrix2 <- as.matrix(rbind(matrix(rep(c(rep(1,7),rep(-1,3),rep(1,5),rep(-1,5),rep(1,5),rep(-1,5)),7),byrow = T, nrow = 7),
                                 matrix(rep(c(rep(-1,7),rep(1,3),rep(-1,5),rep(1,5),rep(-1,5),rep(1,5)),3),byrow = T, nrow = 3),
                                 matrix(rep(c(rep(1,7),rep(-1,3),rep(1,5),rep(-1,5),rep(1,5),rep(-1,5)),5),byrow = T, nrow = 5),
                                 matrix(rep(c(rep(-1,7),rep(1,3),rep(-1,5),rep(1,5),rep(-1,5),rep(1,5)),5),byrow = T, nrow = 5),
                                 matrix(rep(c(rep(1,7),rep(-1,3),rep(1,5),rep(-1,5),rep(1,5),rep(-1,5)),5),byrow = T, nrow = 5),
                                 matrix(rep(c(rep(-1,7),rep(1,3),rep(-1,5),rep(1,5),rep(-1,5),rep(1,5)),5),byrow = T, nrow = 5)))
signs.matrix3 <- as.matrix(rbind(matrix(rep(c(rep(1,7),rep(-1,3)),7),byrow = T, nrow = 7),
                                 matrix(rep(c(rep(-1,7),rep(1,3)),3),byrow = T, nrow = 3)))


# matrix Sig3 completed:
Sig3[1:(2*mark1),1:(2*mark1)] <- BL1*signs.matrix2
Sig3[31:(30+2*mark2),31:(30+2*mark2)] <- BL2*signs.matrix3
Sig3[41:(40+2*mark3),41:(40+2*mark3)] <- BL3*signs.matrix2
Sig3[71:(70+2*mark4),71:(70+2*mark4)] <- BL4*signs.matrix3

diag(Sig3)[81:(D-1)] <- 1    # diagonal elements outside the blocks are set to 1

# Check if Sig3 is symmetrical and positive definite
isSymmetric(Sig3)
is.positive.definite(Sig3)


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


### PCA PBs algorithm
ptm <- proc.time()

sol1b = foreach(i = 1:S, .combine = "comb", .multicombine=TRUE, .verbose = T, .init=list(list()), 
                .packages = c("compositions","pls","MASS"))%dopar%{
                  
                  set.seed(i,"L'Ecuyer-CMRG")
                  # data simulaton
                  d <- data_simulation(n, D, mu = mu1, codes = codes1, V = V1, r = D/10, Sig=Sig3, typeSig = "Sig3")   
                  X <- d$X
                  y <- d$y
                  nbal <- D-1

                  # calculation
                  rmsep.vec1 <- c()
                  for(m in 1:nbal){
                    rmsep.vec1[m] <- PBs_PCA_CV(X,y,nbalances=m,R=R, K=K, n=n)} 
                  
                  # store results 
                  list(rmsep.vec1)
                }
proc.time()-ptm 

### PLS PBs algorithm
ptm <- proc.time()

sol2b = foreach(i = 1:S, .combine = "comb", .multicombine=TRUE, .verbose = T, .init=list(list()), 
                .packages = c("compositions","pls","MASS"))%dopar%{
                  
                  set.seed(i,"L'Ecuyer-CMRG")
                  # data simulaton
                  d <- data_simulation(n, D, mu = mu1, codes = codes1, V = V1, r = D/10, Sig=Sig3, typeSig = "Sig3") 
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

sol3b = foreach(i = 1:S, .combine = "comb", .multicombine=TRUE, .verbose = T, .init=list(list()), 
                .packages = c("compositions","pls","MASS"))%dopar%{
                  
                  set.seed(i,"L'Ecuyer-CMRG")
                  # data simulaton
                  d <- data_simulation(n, D, mu = mu1, codes = codes1, V = V1, r = D/10, Sig=Sig3, typeSig = "Sig3")   
                  X <- d$X
                  y <- d$y
                  nload <- D
                  
                  # calculation
                  rmsep.vec1 <- c()
                  for(m in 1:nload){
                    rmsep.vec1[m] <- PLS_CV(X,y,nload=m,R=R, K=K, n=n)            
                  }
                  
                  # store results 
                  list(rmsep.vec1)
                }
proc.time()-ptm 

### Selbal
ptm <- proc.time()

sol4b = foreach(i = 1:S, .combine = "comb", .multicombine=TRUE, .verbose = T, .init=list(list()), 
                .packages = c("compositions","pls","MASS"))%dopar%{
                  
                  set.seed(i,"L'Ecuyer-CMRG")
                  # data simulaton
                  d <- data_simulation(n, D, mu = mu1, codes = codes1, V = V1, r = D/10, Sig=Sig3, typeSig = "Sig3")
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
t1b <- ToTable(sol1b,S=S,D=D)  # PCA PBs
t2b <- ToTable(sol2b,S=S,D=D)  # PLS PBs
t3b <- ToTable(sol3b,S=S,D=D)  # PLS
t4b <- ToTableSelbal(sol4b,S=S,D=1) # Selbal




### Ability to correctly identify the biomarkers -------------------------------
# Construct "ideal balances" for comparison 
# -> balances, that consist just of marker variables and nothing else (num, den according to cov matrix)
# -> Fig 6

# Run sim 100 times; just construct 1st PCA_PB, PLS_PB, selbal balance

bal3.PCA <- matrix(NA,nrow=S,ncol=D)
bal3.PLS <- matrix(NA,nrow=S,ncol=D)

set.seed(1554)
for(i in 1:S){
  d <- data_simulation(n, D, mu = mu1, codes = codes1, V = V1, r = D/10, Sig=Sig3, typeSig = "Sig3") 
  X <- d$X
  y <- d$y
  
  # PCA PBs
  modelA <- fBalChip(X)
  bal3.PCA[i,] <- modelA$bal[1,]
  
  # PLS PBs
  modelB <- fBalChip_PLS(X,y,ver="cov")
  bal3.PLS[i,] <- modelB$bal[1,]
  
  print(i)
}


bal3.selbal <- matrix(NA,nrow=S,ncol=D)
set.seed(1554)
tic()
for(i in 1:S){
  d <- data_simulation(n, D, mu = mu1, codes = codes1, V = V1, r = D/10, Sig=Sig3, typeSig = "Sig3")
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
  bal3.selbal[i,] <- logcontr
  print(i)
}
toc()

# Construct "ideal balance" for comparison -> balance, that consist just of marker variables and nothing else (num, den according to cov matrix)

inds2 <- c(rep(T,7),rep(F,3),rep(T,5),rep(F,5),rep(T,5),rep(F,5),rep(T,7),rep(F,3),rep(T,7),rep(F,3),rep(T,5),rep(F,5),rep(T,5),rep(F,5))
posl2 <- 1:70
indpos3 <- posl2[inds2==T]   # indices of variables in the numerator
indneg3 <- posl2[inds2==F]    # indices of variables in the denominator

r3 <- length(indpos3)
s3 <- length(indneg3)

bpos3 <- (1/r3)*sqrt((r3*s3)/(r3+s3))
bneg3 <- (-1/s3)*sqrt((r3*s3)/(r3+s3))


logcontr3 <- rep(0,D)   # balance-logcontrast
logcontr3[indpos3] <- bpos3
logcontr3[indneg3] <- bneg3
logcontr3



