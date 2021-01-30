README
1.Function:
    This code is uesd to predict phenotype value according to combined forecasting method,which is based on prediction of historical reference data according to c kinds of individual predicting methods.
2.Inputting Parameters:
    The code needs four parameters to input:
    (1)traindata,which is a matrix with c columns.It is prediction of historical reference data according to c kinds of individual predicting methods.
    (2)realdata,which is a vector with length equaling to row number of "traindata".It is obeserved value of historical reference data,and it  has the same order with "traindata" in sample number.
    (3)data,which is a matrix with c columns.It is prediction of test data according to c kinds of individual predicting methods.And the column names of "data" need to be the same with "traindata".
    (4)cutoff,which is a threshold value of phenotype.You can get it according to clinical experience.
3.Output:
    You can obtain a prediction vector of test data according to combined forecasting method.
4.Attention:
    R package "Rdonlp2" is required in this code.So you need install the "Rdonlp2" package in R in earlier time.There will be a logarithmic transformation for data in code,so you need make sure parameters for input are positive.
5.Test Data:
    Here we prepare a set of data for testing the code for you:
    Compete code as below in R:

library(GP3CCO)
path="C:\\Users\\jlp\\Desktop\\testdata"
traindata<-read.table(paste(path,"\\traindata.txt",sep=''),head=F)
realdata<-read.table(paste(path,"\\realdata.txt",sep=''),head=F)
data<-read.table(paste(path,"\\data.txt",sep=''),head=F)
y<-Ypre(traindata,realdata,data,0.62)

Then you can get a prediction vector "y". You can check it in set "testdata".
