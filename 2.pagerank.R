load("/home/data/t060414/data/network-ML/result/1.network/hanjunwei_graph/graph.rda")
library(igraph)
load("/home/data/t060414/data/network-ML/result/DNB_pathway.Rdata")
load("/home/data/t060414/data/network-ML/result/multiomics_pathway.Rdata")
pathway<-intersect(pathway_TCGA,pathway_GSE98394)
gene<-c()
for(i in pathway){
  gene<-c(gene,multiomics_pathway[[i]])
}
gene<-unique(gene)
gene<-intersect(gene,V(graph)$name)
####构建子网
subgraph <- induced_subgraph(graph, gene)

nodes<-data.frame(node=V(subgraph)$name,
                  weight=0)
###CTpathway: a CrossTalk-based pathway enrichment analysis method for cancer research 引用这个节点加权的方法
load("/home/data/t060414/data/network-ML/result/2.Immunotherapy_DE.Rdata")
ann<-read.table("/home/data/t060414/data/network-ML/data/ICB/SKCM_response.txt",header=F,sep="\t")
for(i in ann$V1[c(1:9)] ){
  DE<-response_DE[[i]][,c(1,2,5)]
  
  DE[,2]<-2^DE[,2]
  for (j in 1:dim(DE)[1]) {
    if(is.na(DE[j,2])){
      
    }else{
    if(DE[j,2]>=1){
      DE[j,4]<-(1-DE[j,3])*(1-DE[j,2]^(-1/2))
    }else{
      DE[j,4]<-(1-DE[j,3])*(1-DE[j,2]^(1/2))
     }
    }
  }
  names(DE)[2:4]<-c(paste0(i,"_FC"),paste0(i,"_P"),paste0(i,"_score"))
  
  nodes<-merge(nodes,DE,by.x = "node", by.y = "GeneSymbol",all.x=T)
}

nodes$weight<-rowMeans(as.matrix(nodes[,c(5,8,11,14,17,20,23,26,29)]),na.rm=T)

nodes$weight[which(nodes$weight=="NaN")]<-0

load("/home/data/t060414/data/network-ML/data/TILSig/output_TSI_all_result.Rdata")
TIL<-unique(top_TSI_outTab[,c(1,3)])
nodes_new<-merge(nodes[,1:2],TIL,by.x="node",by.y="gene",all.x=T)
nodes_new[which(is.na(nodes_new$TSI)),3]=1

nodes_new$new_weight<-(1-as.numeric(nodes_new$TSI))*nodes_new$weight
# 输出标准化并限制范围后的数据
nodes_new$zscore<- as.data.frame(scale(nodes_new$new_weight, center = min(nodes_new$new_weight), scale = max(nodes_new$new_weight)+0.001 - min(nodes_new$new_weight)))




###已知靶点全设置

IO<-read.table("/home/data/t060414/data/network-ML/data/IO.txt",header=F)

IO_inter<-intersect(IO[,1],V(subgraph)$name)
nodes_new[nodes_new$node %in%  IO_inter,5]<-1

V(subgraph)$weight <- nodes_new$zscore

###计算边的权重  皮尔森相关系数

load("/home/data/t060414/data/network-ML/result/2.cor_edge_weigt.Rdata")
df_edges <- as_data_frame(subgraph, what = "edges")[,-3]
names(df_edges)<-c("row" ,"column" )
edges<-merge(df_edges,edge_all,all=F)
edges$cor_score_use[which(edges$cor_score_use=="NaN")]<-0
E(subgraph)$weight <-abs(edges$cor_score_use)
# 计算PageRank
pr <- page.rank(subgraph)
pagerank<-as.matrix(pr$vector)
save(pagerank,subgraph,nodes,file="/home/data/t060414/data/network-ML/result/2.pagerank.Rdata")
# 输出PageRank值
print(pr$vector)
###CTpathway: a CrossTalk-based pathway enrichment analysis method for cancer research 引用这个网络的方法
