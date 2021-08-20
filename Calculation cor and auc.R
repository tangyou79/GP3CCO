1.This introduction is used to introduce how to calculate the Pearson correlation and  AUC(the area under the receiver operator characteristic curve).R package "pROC" is required to calculate the AUC .So you need install the "pROC" package in R in earlier time.

2.The relevant code can be written as below:
 
library(pROC)
cor<-cor(x,y);
auc<-roc(z,y)$auc[1];

Where x is a vector with length n,and y is a vector with length n.z is a binary vector with length n.
Inputting these code in R,then you can get both the Pearson correlation of vector x and y and the "AUC" of vector z and y.
