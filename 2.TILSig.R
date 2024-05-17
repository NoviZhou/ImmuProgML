###TILSig
install.packages("/home/data/t060414/data/network-ML/code/DealGPL570_0.0.1.tar.gz", repos = NULL, type = "source")
library(DealGPL570)
library(stringr)
library(sva)
library(tibble)
library(dplyr)
library(tidyverse)
library(GenomicFeatures)
library(limma)
library(survival)
library(rtracklayer)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
# 使用DealGPL570处理16个表达谱数据，包含114个细胞系，19种免疫细胞类型
setwd("/home/data/t060414/data/network-ML/data/TILSig")
GSE13906 <- DealGPL570(file = "GSE13906_RAW.tar",type = "geneSymbol")
GSE13906 <- column_to_rownames(GSE13906,"symbol")

GSE23371 <- DealGPL570(file = "GSE23371_RAW.tar",type = "geneSymbol")
GSE23371 <- column_to_rownames(GSE23371,"symbol")

GSE25320 <- DealGPL570(file = "GSE25320_RAW.tar",type = "geneSymbol")
GSE25320 <- column_to_rownames(GSE25320,"symbol")

GSE27291 <- DealGPL570(file = "GSE27291_RAW.tar",type = "geneSymbol")
GSE27291 <- column_to_rownames(GSE27291,"symbol")

GSE27838 <- DealGPL570(file = "GSE27838_RAW.tar",type = "geneSymbol")
GSE27838 <- column_to_rownames(GSE27838,"symbol")

GSE28490 <- DealGPL570(file = "GSE28490_RAW.tar",type = "geneSymbol")
GSE28490 <- column_to_rownames(GSE28490,"symbol")

GSE28698 <- DealGPL570(file = "GSE28698_RAW.tar",type = "geneSymbol")
GSE28698 <- column_to_rownames(GSE28698,"symbol")

GSE28726 <- DealGPL570(file = "GSE28726_RAW.tar",type = "geneSymbol")
GSE28726 <- column_to_rownames(GSE28726,"symbol")

GSE37750 <- DealGPL570(file = "GSE37750_RAW.tar",type = "geneSymbol")
GSE37750 <- column_to_rownames(GSE37750,"symbol")

GSE39889 <- DealGPL570(file = "GSE39889_RAW.tar",type = "geneSymbol")
GSE39889 <- column_to_rownames(GSE39889,"symbol")

GSE42058 <- DealGPL570(file = "GSE42058_RAW.tar",type = "geneSymbol")
GSE42058 <- column_to_rownames(GSE42058,"symbol")

GSE49910 <- DealGPL570(file = "GSE49910_RAW.tar",type = "geneSymbol")
GSE49910 <- column_to_rownames(GSE49910,"symbol")

GSE51540 <- DealGPL570(file = "GSE51540_RAW.tar",type = "geneSymbol")
GSE51540 <- column_to_rownames(GSE51540,"symbol")

GSE59237 <- DealGPL570(file = "GSE59237_RAW.tar",type = "geneSymbol")
GSE59237 <- column_to_rownames(GSE59237,"symbol")

GSE6863 <- DealGPL570(file = "GSE6863_RAW.tar",type = "geneSymbol")
GSE6863 <- column_to_rownames(GSE6863,"symbol")

GSE8059 <- DealGPL570(file = "GSE8059_RAW.tar",type = "geneSymbol")
GSE8059 <- column_to_rownames(GSE8059,"symbol")
#修改列名，仅保留GSE开头的GSE样本编号，去除.CEL后缀
for (i in ls(pattern = "^GSE")) {
  rt <- get(i)
  colnames(rt) <- str_replace_all(colnames(rt),".CEL","")
  assign(i,rt)
}
#从PMID: 28052254文献的补充材料中（文件名:1-s2.0-S2211124716317090-mmc3.xlsx）提取每个样本对应的细胞类型，整理后的文件命名为Immune_cell_line.txt
#读取注释信息
ann <- read.table("Immune_cell_line.txt",header=T,sep="\t",quote="",check.name=F)
#看一下注释信息
head(ann)
#提取16个数据集的注释信息（按GSE号提取，粗略提取）
ann <- ann[which(ann$`Data set`%in%ls(pattern = "^GSE")),]
table(ann$`Data set`)
#提取16个数据集中仅在注释中有的免疫细胞系（按GSM号提取，详细提取）
for (j in ls(pattern = "^GSE")) {
  rt <- get(j)
  rt <- rt[,which(colnames(rt)%in%ann[which(ann$`Data set`==j),"Sample ID"])]
  assign(j,rt)
}
#将数据合并
data <- cbind(GSE13906, GSE23371, GSE25320, GSE27291, GSE27838, GSE28490, GSE28698, 
              GSE28726, GSE37750, GSE39889, GSE42058, GSE49910, GSE51540, 
              GSE59237, GSE6863, GSE8059)
#去除批次效应
batch <- data.frame(batch = rep(c("GSE13906", "GSE23371", "GSE25320", "GSE27291", "GSE27838", "GSE28490", "GSE28698", 
                                  "GSE28726", "GSE37750", "GSE39889", "GSE42058", "GSE49910", "GSE51540", 
                                  "GSE59237", "GSE6863", "GSE8059"),
                                times=c(2,3,4,8,8,5,2,14,8,4,4,26,9,10,3,4)))
modcombat <- model.matrix(~1,data = batch)
data_combat <- as.data.frame(ComBat(dat=as.matrix(data),batch=batch$batch,mod=modcombat))

batch <- data.frame(batch = rep(c("GSE13906", "GSE23371", "GSE25320", "GSE27291", "GSE27838", "GSE28490", "GSE28698", 
                                  "GSE28726", "GSE37750", "GSE39889", "GSE42058", "GSE49910", "GSE51540", 
                                  "GSE59237", "GSE6863", "GSE8059"),
                                times=c(2,3,4,8,8,5,2,14,8,4,4,26,9,10,3,4)))
modcombat <- model.matrix(~1,data = batch)
data_combat <- as.data.frame(ComBat(dat=as.matrix(data),batch=batch$batch,mod=modcombat))

#将去除批次效应后的整合数据输出

write.table(data_combat,"data_combat.txt",quote = F,row.names = T,col.names = T,sep = "\t")
save(data_combat,file="data_combat.Rdata")
## 提取network表达谱数据
load("/home/data/t060414/data/network-ML/result/2.pagerank.Rdata")
#####挑选 top节点####
node=V(subgraph)$name

##最后得到lncRNA表达谱数据
data_node<- data_combat[which(rownames(data_combat)%in%node),]
save(data_node,file="exp_node.Rdata")

# 对每个免疫细胞系的lncRNA进行排序


# 提取lncRNA表达谱数据的免疫细胞系注释信息
ann_114 <- ann[which(ann$`Sample ID`%in%colnames(data_node)),]
table(ann_114$Population)
# 调整表达谱顺序，使得注释信息和表达谱数据的样本顺序一致
loc <- match(ann_114$`Sample ID`,colnames(data_node))
data_node <- data_node[,loc]

# 对于一个lncRNA对应多个相同免疫细胞系，我们取lncRNA在多个相同免疫细胞系表达值的均值来作为后续的排序数值
mean_node<- as.data.frame(t(data_node)) # 将数据转置，行为114个免疫细胞系，列为lncRNA
mean_node$CellType <- ann_114$Population
CellType <- ann_114$Population
mean_node<- aggregate(mean_node,by=list(CellType),FUN = mean) # 取lncRNA在多个相同免疫细胞系中的表达值均值作为该lncRNA在该种免疫细胞系中的表达值
mean_node$CellType <- NULL # 行为19种免疫细胞系（114个细胞系按免疫细胞类型去重后的结果），列为lncRNA
# 对每个lncRNA在所有细胞系的表达谱数据进行排序
mean_node<- column_to_rownames(mean_node,"Group.1")
rank_node<- as.data.frame(t(mean_node))#行为lncRNA，列为19种免疫细胞系
# 本步骤目的是提取出19种免疫细胞系中每种细胞系的lncRNA排序结果
top_outTab <- data.frame()#创建一个空数据框，用于存储每个lncRNA的排序情况
for (m in c(seq(0.01,1,0.01))) {# 按1%为分割间距，将排序结果分为前%1到100%，根据个人需要提取前%多少的基因，原文的提取的5%
  outTab <- data.frame(ID=1:(nrow(rank_node)*m))
  # 对每种免疫细胞系中的lncRNA表达值进行降序排列
  for (i in 1:ncol(rank_node)) {
    top <- rownames(rank_node[order(rank_node[,i],decreasing = T),])[1:(nrow(rank_node)*m)]
    outTab[,i] <- top
  }
  colnames(outTab) <- colnames(rank_node)
  # 这里提供两种方法去提取19种免疫细胞系表达谱的共有lncRNA，分别是交集和并集
  inter <- Reduce(intersect,  list(v1 = outTab[,1],
                                   v2 = outTab[,2],
                                   v3 = outTab[,3],
                                   v4 = outTab[,4],
                                   v5 = outTab[,5],
                                   v6 = outTab[,6],
                                   v7 = outTab[,7],
                                   v8 = outTab[,8],
                                   v9 = outTab[,9],
                                   v10 = outTab[,10],
                                   v11 = outTab[,11],
                                   v12 = outTab[,12],
                                   v13 = outTab[,13],
                                   v14 = outTab[,14],
                                   v15 = outTab[,15],
                                   v16 = outTab[,16],
                                   v17 = outTab[,17],
                                   v18 = outTab[,18],
                                   v19 = outTab[,19]))
  union <- Reduce(union,  list(v1 = outTab[,1],
                               v2 = outTab[,2],
                               v3 = outTab[,3],
                               v4 = outTab[,4],
                               v5 = outTab[,5],
                               v6 = outTab[,6],
                               v7 = outTab[,7],
                               v8 = outTab[,8],
                               v9 = outTab[,9],
                               v10 = outTab[,10],
                               v11 = outTab[,11],
                               v12 = outTab[,12],
                               v13 = outTab[,13],
                               v14 = outTab[,14],
                               v15 = outTab[,15],
                               v16 = outTab[,16],
                               v17 = outTab[,17],
                               v18 = outTab[,18],
                               v19 = outTab[,19]))
  union <- data.frame(top_index=rep(m,length(union)),gene=union) # 取并集用此代码
  #inter <- data.frame(top_index=rep(m,length(inter)),gene=inter) # 取交集用此代码
  top_outTab <- rbind(union,top_outTab) # 取并集用此代码
  #top_outTab <- rbind(inter,top_outTab) # 取交集用此代码
  
}
# 第一列为前百分之多少，阈值为1%到100%，原文选取的5%，第二列是lncRNA
head(top_outTab)

# 计算TSI，提取ihk-lnc

top_TSI_outTab <- data.frame()
for (n in unique(top_outTab$top_index)){
  # 对于前1%到100%的lncRNA均进行TSI计算，大家可以根据自己想设定的阈值进行提取
  TSI <- data.frame(TSI=top_outTab[which(top_outTab$top_index==n),"gene"])
  rownames(TSI) <- TSI$TSI
  # 这步就是实现例子中的步骤
  for (j in top_outTab[which(top_outTab$top_index==n),"gene"]) {
    y=0
    for (m in 1:ncol(rank_node)) {
      x <- 1-(rank_node[j,m]-min(rank_node[j,]))/(max(rank_node[j,])-min(rank_node[j,]))
      y <- sum(y,x)
    }
    TSI[j,1] <- y/(ncol(rank_node)-1)
  }
  TSI$gene <- rownames(TSI)
  TSI <- merge(top_outTab[which(top_outTab$top_index==n),],TSI,by="gene")
  top_TSI_outTab <- rbind(TSI,top_TSI_outTab)
}

# 第一列为lncRNA，第二列为lncRNA的排序情况，第三列为每个lncRNA的TSI值
head(top_TSI_outTab)
# 输出TSI结果
save(top_TSI_outTab,file="output_TSI_all_result.Rdata")
write.table(top_TSI_outTab,"output_TSI_all_result.txt",col.names = T,row.names = F,sep = "\t",quote = F)

a<-unique(top_TSI_outTab[top_TSI_outTab$TSI<0.4,])
a<-unique(a$gene)

