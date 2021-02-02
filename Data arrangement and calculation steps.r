#1.Data Conversion(Use plink to process data)
plink2 --bfile cornell_canine --export A-transpose --out new_file --dog --maf

#2.Dealing with missing values
gen=read.table("/home/lizhuo/dog/dog/data/new_file.traw")#Read data

list <-which(rowSums(is.na(gen)) > 0)#Find the location of missing data

genAll <- gen[-list,]#Delete missing data

write.table(genAll,paste("/home/lizhuo/dog/dog/data/dog4.txt",sep=''),quote=F,row.names=F,col.names = F)#Write the sorted data

#3.Start Calculating
setwd("/home/lizhuo/dog/text")#Set up the workspace
myY=read.table("/home/lizhuo/dog/text/data/PH.txt",head=T)[,2]  #read phenotype file and extract the phenotype vector
rep=100#Set the number of cycles
cycles=5#Randomly divide each column of data into 5 columns for calculation
for(i in 1:rep){#100 copies of phenotypic data, each divided into 5 columns, each column randomly generates 20% missing data
myY_R=matrix(NA,nrow=length(myY), ncol=cycles)
s=data.frame(matrix(sample(1:(floor(length(myY)*0.2)*5)), ncol=5))
for(r in 1:cycles){
  test=s[,r]
  myY_R[-test,r]=myY[-test]
  write.table(myY_R,paste("pheno_reference_rep_",i,".txt",sep=''),quote=F,row.names=F,col.names = F)
  }
}

###7 fold cross validation
require(BGLR,quietly = T)
require(ranger,quietly = T)
require(kernlab,quietly = T)
require(glmnet,quietly = T)
require(NAM,quietly = T)
require(rrBLUP,quietly = T)
require(pls,quietly = T)


gen=read.table("/home/lizhuo/dog/text/data/dog4.txt")##each raw is each individual

#Match phenotype data with genotype data
myY<-myY[myY[,1] %in% gen[,1],]
gen<-gen[gen[,1] %in% myY[,1],]
order=match(myY[,1],gen[,1])
gen=gen[order,]

#Remove id
myY <- myY[,2]
gen <- gen[,-1]

# Kernels
E2 = as.matrix(dist(gen)^2);
G = exp(-(E2/mean(E2)))
eK = eigen(G,symmetric=T)
G = NAM::GRM(gen,TRUE)
##

rep=100
cycles=5

method=c("BayesC","BayesL","BRR","GBLUP","SVM")
for(i in 1:rep){
  myY=read.table(paste("pheno_reference_rep_",i,".txt",sep=''),head=F)
  myY_P=matrix(NA,nrow=nrow(myY), ncol=cycles)
  for(r in 1:cycles){
    y=as.numeric(myY[,r])
    w=which(is.na(y))
 
    # Kernels
    f22b = predict(ksvm(gen,y,type= "eps-svr"),gen)[,1];
    SVM = f22b
   
    #Bayesian
  
  	BayesC = BGLR(rmExistingFiles = F,y,ETA=list(list(X=gen,model='BayesC')),verbose=F)$yHat
  	BayesL = BGLR(rmExistingFiles = F,y,ETA=list(list(X=gen,model='BL')),verbose=F)$yHat
    BRR = BGLR(rmExistingFiles = F,y,ETA=list(list(X=gen,model='BRR')),verbose=F)$yHat
 

    # REML GBLUP
    f10 = mixed.solve(y,K=G)
    GBLUP = c(f10$u[w])

    
    for(j in 1:length(method)){
     
      #Store all predicted data
      myY_Pr[,r]=eval(as.symbol(methodr[j]), envir=.GlobalEnv)
      write.table(myY_Pr,paste(method[j],"_ALL_pheno_prediction_rep_",i,".txt",sep=''),quote=F,row.names=F,col.names = F)

      #Store predicted missing data
      myGS=eval(as.symbol(method[j]), envir=.GlobalEnv)
      myY_P[w,j]=myGS
      write.table(myY_P,paste(method[j],"_pheno_prediction_rep_",i,".txt",sep=''),quote=F,row.names=F,col.names = F)
      
    }  
  }
}
