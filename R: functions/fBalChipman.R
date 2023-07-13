fBalChipman<-function(C,angle=TRUE){
# given a coda set C
# This function returns the balance with D parts
# that is close (angle) to the first PC
# If "angle=FALSE" the Max the var of scores
#
# columns
col<-dim(C)[2]
nbal<-col-1
# clr-transfo
clrC<-log(C) - rowMeans(log(C))

# PCs
#
pcClr<-prcomp(clrC)
#first PC: PC1
pcClr1<-pcClr$rotation[,1]

balsig<-sign(pcClr1)
bal<-matrix(0,nbal,col)
colnames(bal)<-colnames(C)
# balances associated to the PCs
# first bal
bal[1,pcClr1==max(pcClr1)]<-1
bal[1,pcClr1==min(pcClr1)]<--1
numbal=1

# other bal
if (col>2){
numbal=numbal+1
while (numbal<col){
bal[numbal,]<-bal[numbal-1,]
useonly<-(bal[numbal-1,]==0)
bal[numbal,abs(pcClr1)==max(abs(pcClr1[useonly]))]<-balsig[abs(pcClr1)==max(abs(pcClr1[useonly]))]
numbal=numbal+1
}#end while
}#end if

# coefficients & angle
VarSBP<-rep(0,nbal)
for (f in 1:nbal) {
  den<-sum(bal[f,]==-1)
  num<-sum(bal[f,]==1)  
  bal[f,bal[f,]==1]<-sqrt(den/((den+num)*num))
  bal[f,bal[f,]==-1]<--sqrt(num/((den+num)*den))
  # variance of the balance
  VarSBP[f]<-abs(sum(bal[f,]*pcClr1))
}
# log-trasnform
lC<-as.matrix(log(C))
mvar=var(as.vector(lC%*%bal[VarSBP==max(VarSBP),]))

if (!angle) {

# calculate variance in the balance direction
  VarSBP<-rep(0,nbal)
  
  for (i in 1:nbal)
{
  Proj<-as.vector(lC%*%(bal[i,]))
  VarSBP[i]<-var(Proj)
}# end for
mvar=max(VarSBP)
}# end if

# return results
return(list(bal=bal[VarSBP==max(VarSBP),],varbal=mvar))


}