fBPUpChi_PLS<-function(Yp,r3,b,version="cov")   # r3 = response variable
{
  
  # given coda set Yp 
  # and given a balance with some zero
  # return list of PARENT principal balances basis
  # that maximizises the variance
  # searching by the NO-FULL {0, -1,+1} and COMPLETING the
  # SBP using a loop (UP) scheme by the CHIPMAN procedure
  
  npart=ncol(Yp)
  nbal=ncol(Yp)-1

  # to save balances and variances
  Bal=c()
  VarB=c()
  
  usezero<-sum(b==0)
  
  #log-transfo data
  lYp=as.matrix(log(Yp))
  # while it is not the full balance go up
  k=0
  while (usezero>0){
    # new balance
    k<-k+1
    # non-zero in the will be in the denominator
    den<-sum(b!=0)
    # for only one zero we get the full
    if (usezero==1){
          
      b[b!=0]<--sqrt(1/((den+1)*den))
      b[b==0]<-sqrt(den/(den+1))
      
      # VarB<-cbind(VarB,var(as.vector(lYp%*%b)))    # NOTE: original code
      ##-------------------------------------------  # NOTE: update for PLS PBs
      if(version=="cov"){                            
        VarB <- cbind(VarB, abs(cov(r3,lYp%*%b)))
      }else if (version=="cor"){
        VarB <- cbind(VarB, abs(cor(r3,lYp%*%b)))
      }
      ##-------------------------------------------
      Bal<-rbind(Bal,b)
      usezero<-0
    }
    # for more than one zero we explore other {0,+1} combinations 
    else{
      # create the combination by CHIPMAN procedure
      # search the maximum balance
      
      clrC_1 <-log(Yp[,b==0]) - rowMeans(log(Yp[,b==0]))
      clrC <- scale(clrC_1, center=TRUE, scale = FALSE) # NOTE: added line to center data
      
      yC <- r3 - mean(r3)    # NOTE: new line -> center response r3
      # PCs
      #
      pcClr<-pls::plsr(yC~clrC,method="simpls",center=F) 

      #first PC: PC1
      bx<-pcClr$loadings[,1]
      
      #
      # look for change of sign
      if (abs(min(bx))>max(bx)){bx<--bx}
      # force zeros to the other sign
      bx[bx<0]<-0
      # matrix of {0,+1} possibilities
      M<-matrix(0,sum(bx>0),length(bx))
      # sort
      bxsort<-sort(bx,decreasing = TRUE,index.return=TRUE)
      # index
      col<-bxsort$ix
      # create M
      for (i in 1:nrow(M)){
        M[i,col[1:i]]<-abs(bx[col[1:i]])
      }
      # sign
      M<-sign(M)
      
      # create a balance
      balax<-b
      # old non-zero to denominator
      balax[b!=0]<--1
      # search the max variance
      VarSBPx<-matrix(0,1,nrow(M))
      balsx<-c()
      for (i in 1:nrow(M))
      {
        # take one possibility
        balax[b==0]<-M[i,]
        # create the coefficients
        num<-sum(balax==1)
        balax[balax==1]<-sqrt(den/((den+num)*num))
        balax[balax==-1]<--sqrt(num/((den+num)*den))
        balsx=rbind(balsx,balax)
        
        #VarSBPx[i]<-var(as.vector(lYp%*%balax))     # NOTE: original code
        ##-------------------------------------------# NOTE: update for PLS PBs
        if(version=="cov"){                          
          VarSBPx[i]<-abs(cov(r3,lYp%*%balax))
        } else if(version=="cor"){
          VarSBPx[i]<-abs(cor(r3,lYp%*%balax))
        }
        ##-------------------------------------------
      }
      VarB=cbind(VarB,max(VarSBPx))
      Bal=rbind(Bal,balsx[VarSBPx==VarB[k],])
      
      rm(M)
      usezero<-sum(Bal[k,]==0)
      b<-Bal[k,]
      # end else GO UP
    }
    
    # end GO UP while
  }
  # return results
  return(list(bal=Bal,varbal=VarB))
  # end function
}