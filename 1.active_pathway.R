###Active Pathway 识别通路
#####
library(ActivePathways)
mut<-read.csv("/home/data/t060414/data/network-ML/result/logistic_univariate_stage_mutations_all.csv")
names(mut)[10]<-"mut.p"
load("/home/data/t060414/data/network-ML/result/Cor_stage_score.Rdata")
methy<-Cor_age_score[,c(3,2)]
names(methy)[2]<-"methy.p"
load("/home/data/t060414/data/network-ML/result/CNV_stage_multicox.Rdata")
CNV<-result_df_all
names(CNV)[1]<-"gene"
load("/home/data/t060414/data/network-ML/result/mRNA_expression_stage.Rdata")
mRNA<-top.table
mRNA$gene<-row.names(mRNA)
library(clusterProfiler)
library(org.Hs.eg.db)
ens2s = bitr(CNV$gene, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")
CNV<-merge(ens2s,CNV,by.x="ENSEMBL",by.y="gene",all=F)
names(CNV)[10]<-"CNV.P"
tmp<-merge(mut[,c(1,10)],CNV[,c(2,10)],by.x="gene",by.y="SYMBOL",all=T)
tmp<-merge(tmp,methy,by.x="gene",by.y="gene",all=T)


mRNA$ENSEMBL<-row.names(mRNA)
mRNA<-merge(ens2s,mRNA,all=F)
tmp<-merge(tmp,mRNA[,c(2,6)],by.x="gene",by.y="SYMBOL",all=T)
names(tmp)[5]<-"mRNA.p"
merge_p<-tmp[-which(duplicated(tmp$gene)),]
row.names(merge_p)<-merge_p$gene
merge_p<-merge_p[,-1]

#tmp1<-intersect(intersect(which(merge_p$mut.p>0.05),which(merge_p$CNV.P>0.05)),which(merge_p$methy.p>0.05))


#merge_p<-merge_p[-tmp1,]
merge_p <- as.matrix(merge_p)
merge_p[is.na(merge_p)] <- 1 
gmt.file <- "/home/data/t060414/data/network-ML/code/Msigdb/c2.all.v2023.1.Hs.symbols.gmt"
pathway_enrich<-ActivePathways(merge_p, gmt.file, significant = 0.05,geneset.filter = c(25, 500))

gmt.file <- "/home/data/t060414/data/network-ML/code/Msigdb/c5.all.v2023.1.Hs.symbols.gmt"
GO_enrich<-ActivePathways(merge_p, gmt.file, significant = 0.05,geneset.filter = c(25,500))
####
gmt.file <- "/home/data/t060414/data/network-ML/code/Msigdb/h.all.v2023.1.Hs.symbols.gmt"
Hallmark_enrich<-ActivePathways(merge_p, gmt.file, significant = 0.05,geneset.filter = c(25,500))

save(pathway_enrich,GO_enrich,Hallmark_enrich,file="/home/data/t060414/data/network-ML/result/ActivePathway.Rdata")
function_list<-list()

c5_geneset2 <- clusterProfiler::read.gmt( "/home/data/t060414/data/network-ML/code/Msigdb/c5.all.v2023.1.Hs.symbols.gmt")#返回的是表
c5_geneset2 <-c5_geneset2 [c5_geneset2$term %in%GO_enrich$term.id,]




for(i in unique(c5_geneset2$term)){

 gene<- c5_geneset2[which(c5_geneset2$term==i),2]
 function_list[[i]]<-gene
}

c2_geneset2 <- clusterProfiler::read.gmt("/home/data/t060414/data/network-ML/code/Msigdb/c2.all.v2023.1.Hs.symbols.gmt")#返回的是表
c2_geneset2 <-c2_geneset2 [c2_geneset2$term %in% pathway_enrich$term.id,]


for(i in unique(c2_geneset2$term)){
  
  gene<- c2_geneset2[which(c2_geneset2$term==i),2]
  function_list[[i]]<-gene
}

Hallmark_geneset2 <- clusterProfiler::read.gmt("/home/data/t060414/data/network-ML/code/Msigdb/h.all.v2023.1.Hs.symbols.gmt")
Hallmark_geneset2 <- Hallmark_geneset2[Hallmark_geneset2$term %in% Hallmark_enrich$term.id,]

for(i in unique(Hallmark_geneset2$term)){
  
  gene<- Hallmark_geneset2[which(Hallmark_geneset2$term==i),2]
  function_list[[i]]<-gene
}
multiomics_pathway<-function_list
save(multiomics_pathway,file="/home/data/t060414/data/network-ML/result/multiomics_pathway.Rdata")

###ensg##
multiomics_pathway_ensg<-list()
for (i in names(multiomics_pathway)) {
  gene<- multiomics_pathway[[i]]
  eg <- bitr(gene,fromType ='SYMBOL',
             toType = c('ENSEMBL'),
             OrgDb='org.Hs.eg.db')
  multiomics_pathway_ensg[[i]]<-eg$ENSEMBL
}
save(multiomics_pathway_ensg,file="/home/data/t060414/data/network-ML/result/multiomics_pathway_ensg.Rdata")
####
