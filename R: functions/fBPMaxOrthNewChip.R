fBPMaxOrthNewChip<-function(Y,angle=TRUE)
{
  
  # recursion: given a coda set Y
  # return list of principal balances basis
  # that maximizises the variance
  # searching by the NO-FULL and COMPLETING the
  # SBP using loop (UP: fBUpChi.r) and recursive (DOWN) schemes
  # both based on Chipman procedure
  
  numpart=ncol(Y)
  numbal=ncol(Y)-1
  B=c()
  #B=matrix(0,numbal,numpart) to save balances
  V=c()
  #V=matrix(0,1,numbal) to save variances
  
  #first optimal in data set Y
  res<-fBalChipman(Y,angle=angle)
  B<-res$bal
  V<-res$varbal
 # if necessary GO UP to complete
  if (sum(B==0)>0){  
    res<-fBPUpChi(Y,B)
  B=rbind(B,res$bal)
  V=cbind(V,res$varbal)
  }
 # control number of balances added
 if (is.vector(B)) B<-matrix(B,1,length(B)) 
 numbaladd<-nrow(B)-1
   
  ### GO DOWN THE CURRENt LIST AND THE FIRST
  ## first go down from the first optimal balance
  
  usenum<-(B[1,]>0)
  useden<-(B[1,]<0)
 # GO DOWN from numerator of the first optimal balance
  if(sum(usenum)>1){
    resP<-fBPMaxOrthNewChip(Y[,usenum],angle=angle)
    Bx<-matrix(0,length(resP$varbal),numpart)
    Bx[,usenum]<-resP$bal
    B<-rbind(B,Bx)
    V<-cbind(V,resP$varbal)
  }# end if
# GO DOWN from denominator of the first optimal balance
  if(sum(useden)>1){
    resP<-fBPMaxOrthNewChip(Y[,useden],angle=angle)
    Bx<-matrix(0,length(resP$varbal),numpart)
    Bx[,useden]<-resP$bal
    B<-rbind(B,Bx)
    V<-cbind(V,resP$varbal)
  }# end if
  
# REVISIT list of balances added GO UP so as to complete the SBP if necessary GO DOWN by the POSITIVE

   if (numbaladd > 0){
     for (k in 2:(1+numbaladd)){
       usepos=(B[k,]>0)
       if (sum(usepos)>1) {
         resP<-fBPMaxOrthNewChip(Y[,usepos],angle=angle)
         Bx<-matrix(0,length(resP$varbal),numpart)
         Bx[,usepos]<-resP$bal
         B<-rbind(B,Bx)
         V<-cbind(V,resP$varbal)
       }#end if2
      }# end for
     }# end if1

# return results
#
V<-as.matrix(V,1,lenght(V))
#
return(list(bal=B,varbal=V))
  
}