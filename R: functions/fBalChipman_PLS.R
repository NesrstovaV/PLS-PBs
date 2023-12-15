fBalChipman_PLS<-function(C,r2,angle=TRUE, version = "cov"){      # r2 = response variable
# given a coda set C
# This function returns the balance with D parts
# that is close (angle) to the first PC
# If "angle=FALSE" the Max the var of scores
#
# columns
col<-dim(C)[2] 
nbal<-col-1

# clr-transfo
clrC_1 <-log(C) - rowMeans(log(C))
clrC <- scale(clrC_1, center=TRUE, scale = FALSE) # # NOTE: added line to center data

yC <- r2 - mean(r2)    # NOTE: new line -> center response r2


# PCs
#

pcClr<-pls::plsr(yC~clrC,method="simpls",center=F) # NOTE: pls::plsr instead of prcomp; center=F -> already centered data
#first PC: PC1
# pcClr1<-pcClr$rotation[,1]  # rotation = the matrix of variable loadings, ORIGINAL CODE
# NOTE: using results of plsr, name of results remained as "pcClr" to avoid errors when rewriting the code
pcClr1 <- pcClr$loadings[,1]

balsig<-sign(pcClr1)
bal<-matrix(0,nbal,col)
colnames(bal)<-colnames(C)
# balances associated to the PCs
# first bal
bal[1,pcClr1==max(pcClr1)]<-1     # line 1 of "bal" and column of a variable having the highest value of pcClr1
bal[1,pcClr1==min(pcClr1)]<--1    # line 1 of "bal" and column of a variable having the lowest value of pcClr1
numbal=1

# other bal
if (col>2){
numbal=numbal+1
while (numbal<col){    
bal[numbal,]<-bal[numbal-1,]
useonly<-(bal[numbal-1,]==0)
bal[numbal,abs(pcClr1)==max(abs(pcClr1[useonly]))]<-balsig[abs(pcClr1)==max(abs(pcClr1[useonly]))] # on a given row of "bal" fill to the corresponding column (of a variable with the highest abs(pcClr1)) the corresponding sign from "balsig" 
numbal=numbal+1
}#end while
}#end if

# OUTPUT: matrix "(D-1) x D" with -1,0,1

# coefficients & angle
VarSBP<-rep(0,nbal)      # UPDATE: calculating cov (or cor); name kept as VarSBP (so that it is not needed to be overwritten everywhere)
for (f in 1:nbal) {
  den<-sum(bal[f,]==-1)  # how many times there is "-1" in the row "f"
  num<-sum(bal[f,]==1)   # how many times there is "+1" in the row "f"
  bal[f,bal[f,]==1]<-sqrt(den/((den+num)*num))     # calculates values of balances "in cells with 1"
  bal[f,bal[f,]==-1]<--sqrt(num/((den+num)*den))   # calculates values of balances "in cells with -1"
  # variance of the balance:
  #VarSBP[f]<-abs(sum(bal[f,]*pcClr1))      # NOTE: original code
  
  if (version == "cov"){                    # NOTE: for PLS PBs (two possible versions: "cov" or "cor")
    VarSBP[f] <- abs(cov(r2,as.matrix(log(C))%*%t(bal)[,f]))  
  } else if (version=="cor") {
    VarSBP[f] <- abs(cor(r2,as.matrix(log(C))%*%t(bal)[,f]))
  }
    
}
# log-trasnform
lC<-as.matrix(log(C))
#mvar=var(as.vector(lC%*%bal[VarSBP==max(VarSBP),]))   # NOTE: orginal code
##---------------------------------------------------  # NOTE: for PLS PBs (two possible versions: "cov" or "cor")
if (version=="cov"){                                   
  mvar = abs(cov(r2,as.vector(lC%*%bal[VarSBP==max(VarSBP),])))
} else if (version == "cor") {
  mvar = abs(cor(r2,as.vector(lC%*%bal[VarSBP==max(VarSBP),]))) 
}
##---------------------------------------------------
if (!angle) {

# calculate variance in the balance direction
  VarSBP<-rep(0,nbal)
  
  for (i in 1:nbal)
{
  Proj<-as.vector(lC%*%(bal[i,]))   
  #VarSBP[i]<-var(Proj)                                  # NOTE: original code
  ##---------------------------------------------------  # NOTE: for PLS PBs (two possible versions: "cov" or "cor")
  if(version=="cov"){
    VarSBP[i] <- cov(r2,Proj)
  } else if(version=="cor"){
    VarSBP[i] <- cor(r2,Proj)
  }
  ##---------------------------------------------------  
      
}# end for
mvar=max(VarSBP)
}# end if

# return results
return(list(bal=bal[VarSBP==max(VarSBP),],varbal=mvar))


}
