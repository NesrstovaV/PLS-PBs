fBalChip<-function(Xcoda){
  # given  a coda set Xcoda
  # this function calls fBPMaxOrthNewChip for a searching using the algortihm for Constrained PCs
  # and it returns ALL balances with D parts
  # that maximizes the variance
  # The balances are sorted by the percentatge of variance
  #
  # Returns a list: balances and variance
  # Bres= balances
  # Vres= variance of balances
  
  numbal=ncol(Xcoda)-1
  # call the recursive function
  res<-fBPMaxOrthNewChip(Xcoda)
  Bres<-res$bal
  balname<-paste("bal",1:nrow(Bres),sep="")  
  rownames(Bres)<-balname
  colnames(Bres)<-colnames(Xcoda)
  Vres<-res$varbal
  #
  # sort by expl var
  vopt<-res$varbal
  # sort variance
  vsopt<-sort(vopt,decreasing = TRUE,index.return=TRUE)
  #
  # assign variance explained already ordered
  Vres<-vsopt$x
  #
  # assign balances same order
  Bres<-Bres[vsopt$ix,]
  #  
  # return results: balances and variances
  return(list(bal=Bres,varbal=Vres)) 
}