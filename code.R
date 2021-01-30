

Ypre<-function(traindata,realdata,data,cutoff){


#library  Rdonlp2
library(Rdonlp2)

#take "log" transformation to each dataset.
traindata<-log(traindata)
realdata<-log(realdata)   
data<-log(data)
cutoff<-log(cutoff)  

#N is number of rows of traindata. 
N<-nrow(traindata)
c<-ncol(traindata)

#p is initial value of iteration.
p<-c(0.5,0.5,rep(0,c-2))

#par.l and par.u are the lower and upper bound of parameters respectively
par.l<-c(rep(0,c))
par.u<-c(rep(1,c))

#fn is  objective function.
fn=function(L){
  -(1/N)*sum(1/(1+exp(-(sqrt(N)/3)*(L%*%t(traindata)-cutoff)*(realdata-cutoff))))
}  

#A is coefficient of linear constrains.
A<-matrix(c(rep(1,c)),1,byrow=TRUE)

#lin.l and lin.u are left and right side of linear constrains.
lin.l<-1
lin.u<-1

#run the function.
ret=donlp2(p,fn,par.u=par.u,par.l=par.l,
A,lin.l=lin.l,lin.u=lin.u) 

#get the estimated value of parameters.            
coe<-ret$par


#calculate  prediction  of combined forecast model of test data . 
Y<-exp(coe%*%t(data))
}

