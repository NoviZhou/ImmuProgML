####machine learning 
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
ann<-read.table("/home/data/t060414/data/network-ML/data/ICB/SKCM_response.txt",header=F,sep="\t")
load("/home/data/t060414/data/network-ML/data/ICB/data_exprs_ICB_list_30.RData")
load("/home/data/t060414/data/network-ML/data/ICB/data_clinical_ICB_list_30.RData")
####select gene list 
load("/home/data/t060414/data/network-ML/result/2.Stage_GSscore.Rdata")
load("/home/data/t060414/data/network-ML/result/2.pagerank.Rdata")
pagerank<-data.frame(gene=row.names(pagerank),score=pagerank)
pagerank$rankscore<-rank(-pagerank$score)
pagerank$rankscore<-(max(pagerank$rankscore,na.rm=T)-(pagerank$rankscore-1))/max(pagerank$rankscore,na.rm=T)
pagerank_GS.score<-merge(Stage_GSscore,pagerank,by.x="GeneSymbol",by.y="gene",all=T)
pagerank_GS.score$GS.pagerank<-pagerank_GS.score$rankscore*pagerank_GS.score$GSscore

pagerank_GS.score<-pagerank_GS.score[order(pagerank_GS.score$GS.pagerank,decreasing = T),]
pagerank_GS.score<-pagerank_GS.score[1:1283,]#top10%
### data prepare
gene<-pagerank_GS.score$GeneSymbol
for(i in ann$V1 ){
  mydata<-data_exprs_ICB_list_30[[i]]
  gene<-intersect(gene,row.names(mydata))
  
}
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

expr_node<-list()
for(i in ann$V1 [-5]){
  mydata<-data_exprs_ICB_list_30[[i]]
  mydata<- apply( mydata,2,fpkmToTpm)
  mydata<-as.data.frame(mydata[gene,])
  pdata<-data_clinical_ICB_list_30[[i]][,c("sample_use","response_use")]
  names(pdata)[2]<-c("group")
  mydata<-as.data.frame(t(mydata))
  mydata$sample_use<-row.names(mydata)
  expr<-merge(pdata,mydata,all=F)
  expr$group<- as.numeric(ifelse(expr$group=="NonBenefit", 0, 1))#XGBoost算法需要的分类变量为数值型
  expr$group <- as.factor(expr$group)
  if(length(which(is.na(expr$OS.time)))!=0){
    expr<-expr[-which(is.na(expr$OS.time)),]}
  expr_node[[i]]<-expr
  
}
for(i in ann$V1 [5]){
  mydata<-data_exprs_ICB_list_30[[i]]
  mydata<-as.data.frame(mydata[gene,])
  pdata<-data_clinical_ICB_list_30[[i]][,c("sample_use","response_use")]
  names(pdata)[2]<-c("group")
  mydata<-as.data.frame(t(mydata))
  mydata$sample_use<-row.names(mydata)
  expr<-merge(pdata,mydata,all=F)
  expr$group<- as.numeric(ifelse(expr$group=="NonBenefit", 0, 1))#XGBoost算法需要的分类变量为数值型
  expr$group <- as.factor(expr$group)
  if(length(which(is.na(expr$OS.time)))!=0){
    expr<-expr[-which(is.na(expr$OS.time)),]}
  expr_node[[i]]<-expr
  
}

###预测响应
library(randomForest)                       
setwd("/home/data/t060414/data/network-ML/result/2.ML_response")
set.seed(100)
rf_input <- expr_node[["GSE91061"]] #为了节约时间，使用部分数据进行展示
names(rf_input)<-gsub("-",".",names(rf_input))
row.names(rf_input)<-rf_input$sample_use 
rf_input<-rf_input[,-1]
n <- length(names(rf_input))
outTab <- data.frame()
rf_input$group
for (i in 1:(n-1)) {
  mtry_fit <- randomForest(group~., data = rf_input, mtry = i)
  err <- mean(mtry_fit$err.rate)
  print(err)
  outTab <- rbind(outTab,err)
  colnames(outTab) <- "err"
}
# 随机设置决策树的大小，以找到模型内误差基本稳定点所对应的决策树数目
# 选取randomforest –mtry节点值，对应误差最小的节点值为41
which.min(outTab$err)
ntree_fit <- randomForest(group~., data = rf_input, mtry = which.min(outTab$err), ntree = 2000)

plot(ntree_fit)

# 之后选择ntree值，ntree指定随机森林所包含的决策树数目，默认为500；
# 在1500左右时，模型内误差基本稳定，故取ntree = 1500
rf <- randomForest(group ~., data = rf_input, mtry = 41, ntree = 2000, importance = T)
rf
importance <- importance(x=rf)
print(importance)
write.table(importance,"output_importance_RF.txt", col.names = T, row.names = T, sep = "\t", quote = F)
threshold <- 0.05  # 设置阈值
selected_features_RF <- row.names(importance[importance[,4] > threshold,])
pdf("B_RF.pdf", width = 5, height = 5)
varImpPlot(rf)
dev.off()
#########xgboost
library(xgboost)
library(Matrix)
#install.packages("Ckmeans.1d.dp")
#如果此处报错，例如
#Error in install.packages : ERROR: failed to lock directory ‘C:\Users\96251\Documents\R\win-library\4.1’ for modifying
#Try removing ‘C:\Users\96251\Documents\R\win-library\4.1/00LOCK’
#只需要去上述所指定的路径下（例子为C:\Users\96251\Documents\R\win-library\4.1/）删除00LOCK文件夹即可

matrix <- sparse.model.matrix(group ~ .-1, data = rf_input)
label <-as.numeric(as.character(rf_input$group))
fin <- list(data=matrix,label=label) 
dmatrix <- xgb.DMatrix(data = fin$data, label = fin$label) 
# 模型训练
xgb <- xgboost(data = dmatrix,max_depth=6, eta=0.5,  
               objective='binary:logistic', nround=25)
# 重要重要性排序 
xgb.importance <- xgb.importance(matrix@Dimnames[[2]], model = xgb)  
head(xgb.importance)
write.table(xgb.importance, "output_importance_XGBoost.txt", col.names = T, row.names = T, sep = "\t", quote = F)

pdf("XGBoost.pdf", width = 5, height = 5)
xgb.ggplot.importance(xgb.importance)
dev.off()
####Boruta算法
library(Boruta)
library(mlbench)
# 为了节约时间，使用部分数据进行展示
Boruta_input <- rf_input

boruta <- Boruta(group ~ ., data = Boruta_input, doTrace = 2, maxRuns = 500)
print(boruta)

#可视化选择的结果
pdf("A_Boruta.pdf",width = 5, height = 5)
plot(boruta, las = 2, cex.axis = 0.7)
dev.off()
#可视化具体的选择过程
pdf("B_Boruta.pdf", width = 5, height = 5)
plotImpHistory(boruta)
dev.off()
#采用TentativeRoughFix的方法来解析不确定的特征
bor <- TentativeRoughFix(boruta)
print(bor)
#显示最终特征选择的结果
attStats(boruta)
selected_features <- getSelectedAttributes(boruta)
write.table(attStats(boruta), "output_importance_Boruta.txt", col.names = T, row.names = T, sep = "\t", quote = F)

####LASSO
library(tidyverse)
library(glmnet)
source('/home/data/t060414/data/network-ML/data/machine test/msvmRFE.R')   #文件夹内自带
library(VennDiagram)
library(sigFeature)
library(e1071)
library(caret)
library(randomForest)
x <- as.matrix(rf_input[,-1]  )
y<-rf_input[,1]  #把分组信息换成01
fit = glmnet(x, y, family = "binomial", alpha = 1, lambda = NULL)
plot(fit, xvar = "dev", label = TRUE)
cvfit = cv.glmnet(x, y, 
                  nfold=10, #例文描述：10-fold cross-validation
                  family = "binomial", type.measure = "class")
plot(cvfit)
cvfit$lambda.min #查看最佳lambda
# 获取LASSO选出来的特征
myCoefs <- coef(cvfit, s="lambda.min");
lasso_fea <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
(lasso_fea <- lasso_fea[-1])
# 把lasso找到的特征保存到文件
setwd("/home/data/t060414/data/network-ML/result/2.ML_response")
write.csv(lasso_fea,"output_importance_lasso.csv")

a<-intersect(row.names(importance),xgb.importance$Feature)
b<-intersect(a,lasso_fea)
c<-intersect(b,selected_features)####SVM-REF-Algorithm
##SVM
input <- rf_input

#采用五折交叉验证 (k-fold crossValidation）
svmRFE(input, k = 5, halve.above = 100) #分割数据，分配随机数
nfold = 5
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
results = lapply(folds, svmRFE.wrap, input, k=5, halve.above=100) #特征选择
top.features = WriteFeatures(results, input, save=F) #查看主要变量
head(top.features)
#把SVM-REF找到的特征保存到文件
write.csv(top.features,"output_importance_svm.csv")
# 运行时间主要取决于选择变量的个数，一般的电脑还是不要选择太多变量

featsweep = lapply(1:300, FeatSweep.wrap, results, input) #300个变量
save(featsweep,file = "featsweep.RData")
# 画图
no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

#dev.new(width=4, height=4, bg='white')
pdf("B_svm-error.pdf",width = 5,height = 5)
PlotErrors(errors, no.info=no.info) #查看错误率
dev.off()
pdf("B_svm-accuracy.pdf",width = 5,height = 5)
Plotaccuracy(1-errors,no.info=no.info) #查看准确率
dev.off()
# 图中红色圆圈所在的位置，即错误率最低点
which.min(errors) 
SVM_fea0ture<-top.features[1:which.min(errors), "FeatureName"]
a<-intersect(selected_features_RF ,SVM_feature)
b<-intersect(a,selected_features)
c<-intersect(b,xgb.importance$Feature)
d<-intersect(c,lasso_fea)

#all<-c(selected_features_RF ,SVM_feature,xgb.importance$Feature,selected_features,lasso_fea)
#Total<-as.data.frame(table(all))
#Total<-Total[Total$Freq>4,]
save(importance,SVM_feature,xgb.importance,lasso_fea,selected_features,file="/home/data/t060414/data/network-ML/result/2.ML_response/ML_response.Rdata")
