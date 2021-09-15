rm(list=ls())

library(reshape2)
library(plyr)
library(boot)

phase="R2"
datatype="SAV"

#basepath=sprintf("xxx",phase,datatype)
setwd(basepath)

print(basepath)

chaps=c("G","H","K")
fileprefix=as.vector(sapply(seq(1,5), function(x){paste0(x,chaps)}))

nreads=read.table("../ncounts.txt",sep = "\t", row.names = 1, col.names = c("","Counts"))
nreadsmat=matrix(nreads[,1],nrow = 8,ncol = 15)
colnames(nreadsmat)=fileprefix

for(i in fileprefix){
  
  print(sprintf("Processing %s",i))
  
  # importing data
  raw_data=read.table(paste0(i,"_",datatype,".txt"),sep = "\t",header=T,row.names = 1, stringsAsFactors = F)
  counts_data=t(round(apply(raw_data,1,function(x){0.01*x*nreadsmat[,i]+1})))
  
  ## converting to a long format
  counts=melt(counts_data,value.name='Success')
  colnames(counts)[1:2]=c("SAV","Rep")
  counts$Rep=substr(counts$Rep,2,4)
  counts$Treatment=sapply(counts$Rep, function(x){cnd = if(strtoi(substr(x,3,3))<5) "+" else "-"})
  counts$Failures = nreads[counts$Rep,]*10 - counts$Success
  # mapping Treatment variable to replicate
  # counts=merge(merge(counts_success,counts_failure),rep_treatment)
  
  # Fitting GLM
  list_of_models=lapply(unique(counts$SAV),function (S){
    
    glm(cbind(Success,Failures)~Treatment,family=binomial,data=counts[counts$SAV==S,])
    
  })
  
  ## naming model objects by SAV variable
  names(list_of_models)=unique(counts$SAV)
  
  ## getting summary of models: SAV, LR -test and slope coefficient for Treatment predictor
  
  results=do.call('rbind',lapply(names(list_of_models), function (S) {
                    amodel=list_of_models[[S]]
                    data.frame(SAV=S,pval=anova(amodel,test='LRT')[2,"Pr(>Chi)"], beta1=coef(amodel)[2],row.names=NULL)
                  })
  )
  ## adjust p-value using Bonferoni correction/ FDR
  results$pval.bonf=p.adjust(results$pval,method='bonf')
  
  # add average frequencies to raw data
  raw_data=cbind(raw_data, t(apply(raw_data,1,function(x) rbind(mean(x[1:4]),mean(x[5:8]),mean(x),max(x) ) ) ))
  colnames(raw_data)[9:12]=c("Avg.P","Avg.W","Avg.All","Max.all")
  raw_data$SAV=rownames(raw_data)
  
  # Merge raw data with refined results
  refined=merge(results[results$pval.bonf<0.05,], raw_data[raw_data$Max.all>=35,], by="SAV")
  
  refined=refined[with(refined,order(-refined$beta1,refined$pval.bonf)),]
  refined$pval=-log10(refined$pval)
  refined$pval.bonf=-log10(refined$pval.bonf)
  colnames(refined)=sub("pval","-log10.Pval",colnames(refined))
  write.table(format(refined,digits=3),file=sprintf("../GLMfit/%s_%s_GLM.txt",datatype,i),row.names = F,sep="\t",quote = F)
  
  rm(counts,raw_data)
  print(sprintf("Finished %s",i))
}
