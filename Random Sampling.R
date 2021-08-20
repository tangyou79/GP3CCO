setwd("/data")
myY=read.table("CHD_03.txt",head=F)

x <- c(0)


for(line in 1:nrow(myY)){
  temp=0
  for(column in 4:9){ 
    if(myY[line,column]==0){
      temp=temp+1       
    }
    if(temp==6){
      x <- append(x,line)
    }
  }
}


myYout = myY[-x,]

write.table(myYout,paste("CHD_04.txt",sep=''),quote=F,row.names=F,col.names = F)




repA=1
repB=100


setwd("data")
myY=read.table("CHD_04.txt",head=F)
###5 fold cross validation
require(BGLR,quietly = T)
require(ranger,quietly = T)
require(kernlab,quietly = T)

gen=read.table("dog-2.txt")


myY<-myY[myY[,1] %in% gen[,1],]
gen<-gen[gen[,1] %in% myY[,1],]

order=match(myY[,1],gen[,1])
gen=gen[order,]

YYY = myY[,2]

rep=repB
cycles=5
for(i in 1:rep){
myY_R=matrix(NA,nrow=length(YYY), ncol=cycles)
s=data.frame(matrix(sample(1:(floor(length(YYY )*0.2)*5)), ncol=5))
for(r in 1:cycles){
  test=s[,r]
  myY_R[-test,r]=YYY[-test]
  write.table(myY_R,paste("data\\pheno_reference_rep_",i,".txt",sep=''),quote=F,row.names=F,col.names = F)
}
}

repA=1
repB=100
setwd("data")
myY=read.table("CHD_04.txt",head=F)
###5 fold cross validation
require(BGLR,quietly = T)
require(ranger,quietly = T)
require(kernlab,quietly = T)

gen=read.table("dog-2.txt")


myY<-myY[myY[,1] %in% gen[,1],]
gen<-gen[gen[,1] %in% myY[,1],]

order=match(myY[,1],gen[,1])
gen=gen[order,]




#myY <- myY[,2]
gen <- gen[,-1]

gen <- as.matrix(gen)
W = myY[,3:ncol(myY)]

cycles=5





E2 = as.matrix(dist(gen)^2);
G = exp(-(E2/mean(E2)))
eK = eigen(G,symmetric=T)
G = NAM::GRM(gen,TRUE)

method=c("GBLUP")

for(i in repA:repB){
  myY=read.table(paste("data\\pheno_reference_rep_",i,".txt",sep=''),head=F)
  myY_P=matrix(NA,nrow=nrow(myY), ncol=cycles)
  myY_Pr=matrix(NA,nrow=nrow(myY), ncol=cycles)
  for(r in 1:cycles){
    y=as.numeric(myY[,r])
    w=which(is.na(y))
  
   # GBLUP = BGLR(rmExistingFiles = F,y,ETA=list(list(X=W, model="FIXED"),list(X=gen,model='GBLUP')),verbose=F)$yHat
    GBLUP=BGLR( y=y,ETA=list(list(X=W, model="FIXED"),G=list(K=G,model='RKHS')))$yHat

    GBLUPr = GBLUP
    GBLUP = GBLUP[w]
  
  for(j in 1:length(method)){
     
     
      myY_Pr[,r]= GBLUPr
      write.table(myY_Pr,paste("data\\",method[j],"_ALL_pheno_prediction_rep_",i,".txt",sep=''),quote=F,row.names=F,col.names = F)
    
      myGS=eval(as.symbol(method[j]), envir=.GlobalEnv)
      myY_P[w,r]=myGS
      write.table(myY_P,paste("data\\",method[j],"_pheno_prediction_rep_",i,".txt",sep=''),quote=F,row.names=F,col.names = F)
   
    }  
  }
}




method=c("BayesC")

for(i in repA:repB){
  myY=read.table(paste("data\\pheno_reference_rep_",i,".txt",sep=''),head=F)
  myY_P=matrix(NA,nrow=nrow(myY), ncol=cycles)
  myY_Pr=matrix(NA,nrow=nrow(myY), ncol=cycles)
  for(r in 1:cycles){
    y=as.numeric(myY[,r])
    w=which(is.na(y))
  

    BayesC = BGLR(rmExistingFiles = F,y,ETA=list(list(X=W, model="FIXED"),list(X=gen,model='BayesC')),verbose=F)$yHat
    BayesCr = BayesC
    BayesC = BayesC[w]
  
  for(j in 1:length(method)){
     
     
      myY_Pr[,r]= BayesCr
      write.table(myY_Pr,paste("data\\",method[j],"_ALL_pheno_prediction_rep_",i,".txt",sep=''),quote=F,row.names=F,col.names = F)
    
      myGS=eval(as.symbol(method[j]), envir=.GlobalEnv)
      myY_P[w,r]=myGS
      write.table(myY_P,paste("data\\",method[j],"_pheno_prediction_rep_",i,".txt",sep=''),quote=F,row.names=F,col.names = F)
   
    }  
  }
}

