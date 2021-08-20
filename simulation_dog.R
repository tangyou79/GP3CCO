###simulation
setwd("~/Downloads/")
library(data.table)
insert1 <- function(x){
  x=gsub("2","3", x, fixed=TRUE)
  x=gsub("0","4", x, fixed=TRUE)
}
insert2 <- function(x){
  x=gsub("3","0", x, fixed=TRUE)
  x=gsub("4","2", x, fixed=TRUE)
}
breeds=6
#randomly set h2 for different breeds from 0.1 to 0.3
set.seed(47850)
H2=sample(1:4,breeds,replace=T)*0.2
#randomly select 1000 QTN
marker=160469
SNP=sample(1:marker,20)

file_name=c("rottweiler","newfoundland","labrador_retriever","golden_retriever","german_shepherd_dog","english_setter")
for(i in 1:breeds){
  h2=H2[i]
  geno=fread(paste(file_name[i],"_gen.txt",sep=''),head=F)[,-1]
  id=fread(paste(file_name[i],"_gen.txt",sep=''),head=F)[,1]
  geno=t(geno[,..SNP])
  sum_af=apply(geno,1,sum)
  af=sum_af/(ncol(geno)*2)
  index_af=which(af>0.5)
  if(length(index_af)<1)
    index_af=1:nrow(geno)
  geno1=geno[index_af,]
  geno2=geno[-index_af,]
  
  geno1=apply(geno1,2,insert1)
  geno1=apply(geno1,2,insert2)
  geno1=matrix(as.numeric(as.matrix(geno1)),nrow(geno1),ncol(geno1))
  colnames(geno1)=colnames(geno2)
  Z=rbind(geno1,geno2)
  sum_af=apply(Z,1,sum)
  af=sum_af/(ncol(Z)*2)
  index=which(af>0.05)
  Z=Z[index,]
  sum_af=apply(Z,1,sum)
  af=sum_af/(ncol(Z)*2)
  m=nrow(Z)
  n=ncol(Z)
  add=rnorm(m,0,sqrt(1/(2*(af*(1-af)))))
  gene_score=as.matrix(t(Z))%*%as.numeric(add)
  
  var_add=var(gene_score)
  residualvar=(var_add-h2*var_add)/h2
  residual=rnorm(n,0,sqrt(residualvar))
  pheno=gene_score+residual
  write.table(cbind(id,pheno),paste(file_name[i],"_h2=",h2,"_sim_pheno.txt",sep=''),col.names=F,row.names=F,quote=F)
}
