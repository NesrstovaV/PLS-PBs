###############################################
### Auxiliary functions
###############################################

# 1) Additional functions -----------------------------------

# Norm of a vector
norm_vec <- function(x) sqrt(sum(x^2))

# Split data into groups for CV -> train and test sets
split_data <- function(X,y,Sigma=NULL,V=NULL) {
  
  n <- nrow(X)
  D <- ncol(X)
  
  # split data: 4/5 train, 1/5 test 
  n_train <- (4/5)*n
  n_test <- (1/5)*n
  v <- 1:nrow(X)
  vysl <- split(v, sample(rep(1:2, c(n_train, n_test))))
  
  vtrain <- vysl$`1`   # indices of "train rows"
  vtest <- vysl$`2`  # indices of "test rows"
  
  # assign Xtest, Xtrain and ytest, ytrain
  Xtest <- X[vtest,]  
  ytest <- y[vtest]
  
  Xtrain <- X[vtrain,]
  ytrain <- y[vtrain]
  
  result <- list("Xtrain" = Xtrain, "ytrain" = ytrain, "Xtest" = Xtest, "ytest" = ytest, "ind_train" = vtrain, "ind_test" = vtest)

}


# RMSEP directly with reziduals
my_rmsep_2 <- function(rezidua,n){
  sqrt(sum(rezidua^2)/n)
}


# 2) Function to simulate data with different covariance matrices - different marker patterns ------------

data_simulation <- function(n, D, mu, V, r = D/10, codes, Sig, typeSig = c("Sig1","Sig2","Sig3")){   
  # n = number of observations
  # D = number of variables
  # mu = matrix of zeros
  # V = matrix of the basis
  # r -> markers associated to the response (10+,10-); this is used in the first and second setting (blocks of markers have the same size)
  # codes = SBP
  # -> n, D, mu, V, codes SAME for all three settings
  # Sig = covariance matrix (for the simulation, we choose from 3 options)
  
  # ilr
  Z <- mvrnorm(n, mu, Sigma = Sig)
  
  # coda
  X <- as.matrix(as.data.frame(acomp(exp(Z%*%t(V))))) 
  
  # Betas and indices of markers 
  if(typeSig=="Sig1"){
    ru <- runif(2*r,0.1,1)  # length(ru1) = number of markers
    beta <- rep(0, D-1)
    signs <- c(rep(1,7),rep(-1,3),rep(1,5),rep(-1,5))   # betas: pattern of signs
    beta[c(1:20)] <- ru*signs   # positive betas: first and fourth block (=sets of markers as there is actually just 1 block altogether)
  } 
  else if(typeSig=="Sig2"){
    ru <- runif(2*r*4,0.1,1)  # length(ru1) = number of markers; *4 four blocks
    signs <- c(rep(1,7),rep(-1,3),rep(1,5),rep(-1,5))   # betas: pattern of signs
    beta = rep(0, D-1)
    beta[c(1:60)] <- ru[1:60]*signs  # positive betas: first block (LAST BLOCK NOT RELATED TO THE RESPONSE) 
  } else if(typeSig=="Sig3"){
    ru <- runif((15+5+15+5)*2,0.1,1)  # last marker block not related to the response
    signs <- c(rep(1,7),rep(-1,3),rep(1,5),rep(-1,5),rep(1,5),rep(-1,5),rep(1,7),rep(-1,3),rep(1,7),rep(-1,3),rep(1,5),rep(-1,5),rep(1,5),rep(-1,5))   # betas: pattern of signs
    beta <- rep(0, D-1)
    beta[1:70] <- ru[1:70]*signs   # LAST BLOCK NOT RELATED TO THE RESPONSE
  } 
  
  # response
  eps <- rnorm(n)      
  y <- Z%*%beta+eps    
  
  results <- list("X" = X, "y" = y, "Sig" = Sig, "Z" = Z)
}



# 3) Functions for CV -----------------------

# a) PCA PBs
function_PCA_PBs <- function(Xtrain, ytrain, Xtest, ytest,intercept=TRUE,nbalance=(ncol(Xtrain)-1)){
  # Xtrain, ytrain = training data
  # Xtest, ytest = test data
  # intercept = whether to consider intercept; default = TRUE
  # nbalance = how many balances do we choose for CV, default = all (i.e. D-1)
  
  # clr+center
  train_clr <-log(Xtrain) - rowMeans(log(Xtrain))
  train_clr_C <- scale(train_clr, center=TRUE, scale = FALSE)
  Xtrain_C <- as.data.frame(compositions::clrInv(train_clr_C))   # back from clr
  ytrain_C <- ytrain - mean(ytrain)
  
  test_clr <- log(Xtest) - rowMeans(log(Xtest))
  test_clr_C <- scale(test_clr, center=TRUE, scale = FALSE)
  Xtest_C <- as.data.frame(compositions::clrInv(test_clr_C))   # back from clr
  ytest_C <- ytest - mean(ytest)
  
  
  # TRAIN data:
  train <- fBalChip(Xtrain)        # model "PBs with PCA"
  balance_train <- train$bal[c(1:nbalance),]      # balances (dim = (number of balances)x(number of variables)) -> rows = balances, cols = variables
  
  if(nbalance>1){   
    soucin <- as.matrix(log(Xtrain_C))%*%t(balance_train)   # dim = (number of observations in Xtrain)x(number of balances) -> rows = observations, cols = balances
  } else if(nbalance==1){
    soucin <- as.matrix(log(Xtrain_C))%*%balance_train      # for the case when only 1 balance is chosen to be calculated
  }
  
  Lmodel <- lm(ytrain_C ~ soucin)
  Beta_train <- Lmodel$coefficients                      # Beta from model with train balances 
  
  # TEST data:
  if(nbalance>1){
    balance_test <- as.matrix(log(Xtest_C))%*%t(balance_train)  # Balances for test data; dim of the result is (nrow Xtest)x(# balances)  
  } else if(nbalance==1){
    balance_test <- as.matrix(log(Xtest_C))%*%balance_train     # for the case when only 1 balance is used
  }
  
  
  # Intercept: TRUE/FALSE
  if(intercept==TRUE){
    ones <- c(rep(1,nrow(balance_test)))
    balance_test <- cbind(ones,balance_test)
    y_hat_test <- balance_test%*%Beta_train      # Intercept = TRUE
  } else{
    balance_test <- balance_test
    y_hat_test <- balance_test%*%Beta_train[-1] # Intercept = FALSE
  }
  
  rezid_test <- ytest_C-y_hat_test

  
  vysledek <- list("logcontr_B_train"=balance_train,"balance_train"=soucin,"Beta_train"=Beta_train,
                   "balance_test"=balance_test,"y_hat_test"=y_hat_test,"rezids_test"=rezid_test)  
}


PBs_PCA_CV <- function(X,y,nbalances, Sigma,V,R, K, n) {
  
  # This function returns residuals for each fold in each CV run
  # For each fold, data will be divided into TRAIN and TEST; TRAIN data can be contaminated

  runs <- list()     # final result of simulation
  outputs <- list()  # to save steps of simulation (each run with 5 folds)
  for(i in 1:R){
    for(j in 1:K){   # for-loop to calculate all "folds"
      
      # split data into TRAIN and TEST
      spl <- split_data(X,y,Sigma=Sigma,V=V)  
      Xtrain <- spl$Xtrain
      ytrain <- spl$ytrain
      Xtest <- spl$Xtest
      ytest <- spl$ytest
      
      outputs[[j]] <- function_PCA_PBs(Xtrain,ytrain,Xtest,ytest,intercept=TRUE, nbal = nbalances)
      runs[[i]] <- outputs
    }
  }
  
  # take all residuals: each component of a list contains a vector of all residuals of 1 CV run (i.e. residuals of all 5 folds together in 1 comp. of a list)
  rezid <- list()
  for(i in 1:R){
    rezid[[i]] <- list(runs[[i]][[1]][["rezids_test"]],runs[[i]][[2]][["rezids_test"]],
                       runs[[i]][[3]][["rezids_test"]],runs[[i]][[4]][["rezids_test"]],
                       runs[[i]][[5]][["rezids_test"]])
  }
  
  results <- rezid
  
}  


# b) PLS PBs
function_PLS_PBs <- function(Xtrain, ytrain, Xtest, ytest,intercept=TRUE,nbalance=(ncol(Xtrain)-1)){
  # Xtrain, ytrain = training data
  # Xtest, ytest = test data
  # intercept = whether to consider intercept; default = TRUE
  # nbalance = how many balances do we choose for CV, default = all (i.e. D-1)
  
  # clr+center
  train_clr <-log(Xtrain) - rowMeans(log(Xtrain))
  train_clr_C <- scale(train_clr, center=TRUE, scale = FALSE)
  Xtrain_C <- as.data.frame(compositions::clrInv(train_clr_C))   # back to clr
  ytrain_C <- ytrain - mean(ytrain)
  
  test_clr <- log(Xtest) - rowMeans(log(Xtest))
  test_clr_C <- scale(test_clr, center=TRUE, scale = FALSE)
  Xtest_C <- as.data.frame(compositions::clrInv(test_clr_C))   # back to clr
  ytest_C <- ytest - mean(ytest)
  
  
  # for TRAIN data
  train_pls <- fBalChip_PLS(Xtrain,ytrain,ver="cov")        # model "PBs with PLS" -> max. covariance
  balance_train <- train_pls$bal[c(1:nbalance),]      # bilance (dim = (# balances)x(# variables)) -> rows = balances, columns = variables
  
  if(nbalance>1){   
    soucin <- as.matrix(log(Xtrain_C))%*%t(balance_train)   # dim = (# observations in Xtrain)x(# balances) -> rows = observations, columns = balances
  } else if(nbalance==1){
    soucin <- as.matrix(log(Xtrain_C))%*%balance_train       # for the case when only 1 balance is used
  }
  
  Lmodel <- lm(ytrain_C ~ soucin)
  Beta_train <- Lmodel$coefficients                      # Beta from the model with training data
  
  # for TEST data:
  if(nbalance>1){
    balance_test <- as.matrix(log(Xtest_C))%*%t(balance_train)  # Balances for test data; dim of the result is (nrow Xtest)x(# balances)  
  } else if(nbalance==1){
    balance_test <- as.matrix(log(Xtest_C))%*%balance_train     # for the case when only 1 balance is used
  }
  
  # include INTERCEPT: yes/no
  if(intercept==TRUE){
    ones <- c(rep(1,nrow(balance_test)))
    balance_test <- cbind(ones,balance_test)
    y_hat_test <- balance_test%*%Beta_train      # Intercept = TRUE
  } else{
    balance_test <- balance_test
    y_hat_test <- balance_test%*%Beta_train[-1] # Intercept = FALSE
  }
  
  rezid_test <- ytest_C-y_hat_test
  
  results <- list("logcontr_B_train"=balance_train,"balance_train"=soucin,"Beta_train"=Beta_train,
                   "balance_test"=balance_test,"y_hat_test"=y_hat_test,"rezids_test"=rezid_test)   # "RMSEP_test"=rmsep_test
  
}


PBs_PLS_CV <- function(X,y,nbalances, Sigma,V,R, K, n) {
  
  # This function returns residuals for each fold in each CV run
  # For each fold, data will be divided into TRAIN and TEST
  
  runs <- list()     # final result of simulation
  outputs <- list()  # to save steps of simulation (each run with 5 folds)
  for(i in 1:R){
    for(j in 1:K){   # for-loop to calculate all "folds"
      
      # split data into TRAIN and TEST
      spl <- split_data(X,y,Sigma=Sigma,V=V)  
      Xtrain <- spl$Xtrain
      ytrain <- spl$ytrain
      Xtest <- spl$Xtest
      ytest <- spl$ytest
      
      outputs[[j]] <- function_PLS_PBs(Xtrain,ytrain,Xtest,ytest,intercept=TRUE, nbal = nbalances)
      runs[[i]] <- outputs
    }
  }
  
  # take all residuals: each component of a list contains a vector of all residuals of 1 CV run (i.e. residuals of all 5 folds together in 1 comp. of a list)
  rezid <- list()
  for(i in 1:R){
    rezid[[i]] <- list(runs[[i]][[1]][["rezids_test"]],runs[[i]][[2]][["rezids_test"]],
                       runs[[i]][[3]][["rezids_test"]],runs[[i]][[4]][["rezids_test"]],
                       runs[[i]][[5]][["rezids_test"]])
  }

  results <- rezid
  
}  


# c) standard PLS
function_PLS <- function(Xtrain, ytrain, Xtest, ytest,intercept=TRUE,nload=ncol(Xtrain)){
  # Xtrain, ytrain = training data
  # Xtest, ytest = test data
  # intercept = whether to consider intercept; default = TRUE
  # nload = how many loadings do we choose for CV, default = all (i.e. D)
  
  # for TRAIN data 
  
  clrC_1 <-log(Xtrain) - rowMeans(log(Xtrain))

  
  df_train <- as.data.frame(cbind(ytrain,clrC_1))
  
  train_pls <- pls::plsr(ytrain~., data = df_train, method="simpls", center=T, ncomp = nload)       # model "PBs with PLS"
  
  
  clr_test <-log(Xtest) - rowMeans(log(Xtest))
  y_hat_test <- predict(train_pls, clr_test, ncomp = nload)
  
  
  rezid_test <- ytest-y_hat_test
  
  results <- list("y_hat_test"=y_hat_test,"rezids_test"=rezid_test)
  
  
}


PLS_CV<- function(X,y,nload,R, K, n) {

  # This function returns residuals for each fold in each CV run
  # For each fold, data will be divided into TRAIN and TEST
  
  runs <- list()     # final result of simulation
  outputs <- list()  # to save steps of simulation (each run with 5 folds)
  for(i in 1:R){
    for(j in 1:K){   # for-loop to calculate all "folds"
      
      # split data into TRAIN and TEST
      spl <- split_data(X,y)  
      Xtrain <- spl$Xtrain
      ytrain <- spl$ytrain
      Xtest <- spl$Xtest
      ytest <- spl$ytest
      
      outputs[[j]] <- function_PLS(Xtrain,ytrain,Xtest,ytest,intercept=TRUE, nload = nload)
      runs[[i]] <- outputs
    }
  }
  
  # take all residuals: each component of a list contains a vector of all residuals of 1 CV run (i.e. residuals of all 5 folds together in 1 comp. of a list)
  rezid <- list()
  for(i in 1:R){
    rezid[[i]] <- list(runs[[i]][[1]][["rezids_test"]],runs[[i]][[2]][["rezids_test"]],
                       runs[[i]][[3]][["rezids_test"]],runs[[i]][[4]][["rezids_test"]],
                       runs[[i]][[5]][["rezids_test"]])
  }
  
  results <- rezid
}  



# d) Selbal
function_selbal <- function(Xtrain, ytrain, Xtest, ytest, intercept = TRUE){
  # Xtrain, ytrain = training data
  # Xtest, ytest = test data
  
  # for TRAIN data
  train_selbal <- selbal::selbal(Xtrain,ytrain,draw=F)   # calculate a balance using selbal
  bilance <- train_selbal$balance.values                 # resulting balance, not a logcontrast
  
  # to obtain a logcontrast
  num <- train_selbal$numerator   # components in numerator
  den <- train_selbal$denominator  # components in denominator
  
  r <- length(num)
  s <- length(den)
  
  bpos <- (1/r)*sqrt((r*s)/(r+s))
  bneg <- (-1/s)*sqrt((r*s)/(r+s))
  
  cols <- colnames(Xtrain)
  
  indpos <- match(num,cols)   # indices of variables in the numerator
  indneg <- match(den,cols)   # indices of variables in the denominator
  
  logcontr_train <- rep(0,ncol(Xtrain))   # balance-logcontrast
  logcontr_train[indpos] <- bpos
  logcontr_train[indneg] <- bneg
  
  # Linear model
  soucin <- as.matrix(log(Xtrain))%*%logcontr_train
  
  Lmodel <- lm(ytrain ~ soucin)
  Beta_train <- Lmodel$coefficients                     # Beta from the model with training data
  
  # for TEST data:
  logcontr_test <- as.matrix(log(Xtest))%*%logcontr_train   
  
  # include INTERCEPT: yes/no
  if(intercept==TRUE){
    ones <- c(rep(1,nrow(logcontr_test)))
    logcontr_test <- cbind(ones,logcontr_test)
    y_hat_test <- logcontr_test%*%Beta_train      # Intercept = TRUE
  } else{
    logcontr_test <- logcontr_test
    y_hat_test <- logcontr_test%*%Beta_train[-1] # Intercept = FALSE
  }
  
  rezid_test <- ytest-y_hat_test
  
  results <- list("logcontr_train"=logcontr_train,"soucin"=soucin,"Beta_train"=Beta_train,
                   "logcontr_test"=logcontr_test,"y_hat_test"=y_hat_test,"rezids_test"=rezid_test) 
  
}


Selbal_CV <- function(X,y,R,K,n){
  
  runs <- list()     # final result of simulation
  outputs <- list()  # to save steps of simulation (each run with 5 folds)
  for(i in 1:R){
    for(j in 1:K){   # for-loop to calculate all "folds"
      
      # split data into TRAIN and TEST
      spl <- split_data(X,y)  
      Xtrain <- spl$Xtrain
      ytrain <- spl$ytrain
      Xtest <- spl$Xtest
      ytest <- spl$ytest
      
      outputs[[j]] <- function_selbal(Xtrain,ytrain,Xtest,ytest,intercept=TRUE)
      runs[[i]] <- outputs
    }
  }
  
  # take all residuals: each component of a list contains a vector of all residuals of 1 CV run (i.e. residuals of all 5 folds together in 1 comp. of a list)
  rezid <- list()
  for(i in 1:R){
    rezid[[i]] <- list(runs[[i]][[1]][["rezids_test"]],runs[[i]][[2]][["rezids_test"]],
                       runs[[i]][[3]][["rezids_test"]],runs[[i]][[4]][["rezids_test"]],
                       runs[[i]][[5]][["rezids_test"]])
  }
  
  results <- rezid
}

# 4) Functions to handle the results of CV -------------------------------------


ToTable <- function(vysl,S,D){
  
  means <- vector(mode = "list", length = S)
  sublist <- vector(mode="list", length = D-1)
  
  for(i in 1:S){
    means[[i]] <- sublist
  }
  
  for(i in 1:S){
    for(j in 1:(D-1)){
      means[[i]][[j]] <- lapply(vysl[[1]][[i]][[j]],FUN = function (x) my_rmsep_2(x,n=length(x)))
    }
  }
  
  MEANS <- vector(mode="list", length=S)
  for(i in 1:S){
    MEANS[[i]] <- sublist
  }
  
  for(i in 1:S){
    for(j in 1:(D-1)){
      MEANS[[i]][[j]] <- mean(unlist(means[[i]][[j]]))
    }
  }
  
  result <- matrix(unlist(MEANS),nrow=S,ncol=D-1,byrow = T)
  
}




ToTableSelbal <- function(vysl,S,D=1) {
  means <- vector(mode = "list", length = S)
  sublist <- vector(mode="list", length = 1)
  
  for(i in 1:S){
    means[[i]] <- sublist
  }
  
  for(i in 1:S){
    means[[i]] <- lapply(vysl[[1]][[i]][[1]],FUN = function (x) my_rmsep_2(x,n=length(x)))
  }
  
  MEANS <- vector(mode="list", length=S)
  for(i in 1:S){
    MEANS[[i]] <- sublist
  }
  
  
  for(i in 1:S){
    MEANS[[i]] <- mean(unlist(means[[i]]))
  }
  
  
  result <- matrix(unlist(MEANS),nrow=S,ncol=1,byrow = T)
  
}



