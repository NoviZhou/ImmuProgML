####免疫治疗反应 #####差异表达分析为
load("/home/data/t060414/data/network-ML/data/ICB/data_clinical_ICB_list_30.RData")
ann<-read.table("/home/data/t060414/data/network-ML/data/ICB/SKCM_response.txt",header=F,sep="\t")
load("/home/data/t060414/data/network-ML/data/ICB/data_exprs_ICB_list_30.RData")
library(limma)
response_DE<-list()

for(i in ann$V1[3] ){
  mydata<-data_exprs_ICB_list_30[[i]]
  pdata<-data_clinical_ICB_list_30[[i]]
  unique(pdata$benefit_used)
  #pdata<- pdata[-which( pdata$benefit_used=="None"),]
  sample<-intersect(names(mydata),pdata$geo_accession)
  mydata<-mydata[,sample]
  group <- as.matrix(pdata$response_use)
  rownames(group)<- pdata$geo_accession
  group<-group[sample,]
  
  design <- model.matrix(~0+group)
  colnames(design)=levels(factor(group))
  rownames(design)=rownames(group)
  contrast.matrix<-makeContrasts(paste(c("Benefit","NonBenefit"),collapse = "-"),levels = design)
  myMeanFun<-function(x){
    tapply(as.double(x),group,mean)  
  }
  
  ###calculate cancer and normal mean TPM
  meangroup <- t(apply(t(mydata),2,myMeanFun))
  
  fit <- lmFit(log2(mydata+1),design)
  
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  fit2 <- eBayes(fit2)  
  tempOutput = topTable(fit2, coef=1, n=Inf)
  nrDEG = na.omit(tempOutput)
  
  allresult <- nrDEG
  allresultfile <- data.frame(GeneSymbol=rownames(allresult),allresult)
  response_DE[[i]]<- allresultfile
}

#y<-row.names(mydata)
#gene<-unlist(lapply(y,function(y) strsplit(as.character(y)," /// ")[[1]][1]))
#gene<-unlist(lapply(y,function(y) strsplit(as.character(y),";")[[1]][1]))
#row.names(mydata)<- gene 
for(i in ann$V1[c(1,2,4:9)] ){
  mydata<-data_exprs_ICB_list_30[[i]]
 
  
  pdata<-data_clinical_ICB_list_30[[i]]
  unique(pdata$benefit_used)
  names(mydata)
  #pdata<- pdata[-which( pdata$benefit_used=="None"),]
  sample<-intersect(names(mydata),pdata$sample_use)
  mydata<-mydata[,sample]
  group <- as.matrix(pdata$response_use)
  rownames(group)<- pdata$sample_use
  group<-group[sample,]
  
  design <- model.matrix(~0+group)
  colnames(design)=levels(factor(group))
  rownames(design)=rownames(group)
  contrast.matrix<-makeContrasts(paste(c("Benefit","NonBenefit"),collapse = "-"),levels = design)
  myMeanFun<-function(x){
    tapply(as.double(x),group,mean)  
  }
  
  ###calculate cancer and normal mean TPM
  meangroup <- t(apply(t(mydata),2,myMeanFun))
  
  fit <- lmFit(log2(mydata+1),design)
  
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  fit2 <- eBayes(fit2)  
  tempOutput = topTable(fit2, coef=1, n=Inf)
  nrDEG = na.omit(tempOutput)
  
  allresult <- nrDEG
  allresultfile <- data.frame(GeneSymbol=rownames(allresult),allresult)
  response_DE[[i]]<- allresultfile
}

save(response_DE,file="/home/data/t060414/data/network-ML/result/2.Immunotherapy_DE.Rdata")
###
