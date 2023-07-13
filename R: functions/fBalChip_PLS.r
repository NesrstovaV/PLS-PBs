fBalChip_PLS<-function(Xcoda, ycoda, version = "cov"){    # ycoda = response variable
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
  res<-fBPMaxOrthNewChip_PLS(Xcoda,ycoda,version=version)   # NOTE: "fBPMaxOrthNewChip_PLS" for PLS PBs
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
  
  #return(list(bal=Bres,varbal=Vres))    # NOTE: original code
  if(version=="cov"){                    # NOTE: for PLS PBs
    return(list(bal=Bres,cov=Vres))
  } else if(version=="cor"){
    return(list(bal=Bres,cor=Vres))
  }
}

## NOTE: for PLs PBs -> balances sorted according to max(cov) or max(cor)