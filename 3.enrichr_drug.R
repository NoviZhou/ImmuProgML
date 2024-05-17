####Enrichr 结果####输入了前100个gene
CMap<-read.table("/home/data/t060414/data/network-ML/result/4.drug/Old_CMAP_up_table.txt",sep="\t",header=T)
xsum<-read.csv("/home/data/t060414/data/network-ML/result/4.drug/CMAP.drug.csv")#CMap 新的输入结果
LINC<-read.table("/home/data/t060414/data/network-ML/result/4.drug/LINCS_L1000_Chem_Pert_Consensus_Sigs_table.txt",sep="\t",header=T)
LINC<-LINC[LINC$Adjusted.P.value<0.001,]
IDG<-read.table("/home/data/t060414/data/network-ML/result/4.drug/IDG_Drug_Targets_2022_table.txt",sep="\t",header=T)
IDG<-IDG[IDG$Adjusted.P.value<0.001,]
CMap$Term<- tolower(CMap$Term)
LINC$Term<- tolower(LINC$Term)
IDG$Term<- tolower(IDG$Term)
drug<-intersect(IDG$Term,LINC$Term)
xsum$id<- tolower(xsum$id)
xsum<-xsum[xsum$score<0,]
#intersect(intersect(IDG$Term,LINC$Term),CMap$Term)
#intersect(intersect(IDG$Term,LINC$Term),xsum$id)
# drug
#[1] "erlotinib"   "sorafenib"   "docetaxel"   "midostaurin" "ixabepilone" "cabazitaxel"
#[7] "vandetanib"  "vinorelbine" "vincristine" "gefitinib"   "sunitinib"   "bosutinib"  
#[13]"vinblastine" "nelfinavir"  "pazopanib"
gene<-c("CXCL8","ANXA1","TNNT1","PLK1","NMU","SFN","LAMC2","MMP9","TUBA4A","TUBB2B","TUBB3")
pagerank_GS.score<-pagerank_GS.score[order(pagerank_GS.score$GS.pagerank,decreasing = T),]
input<-pagerank_GS.score[1:100,1]#输入enrichr的基因
load("/home/data/t060414/data/network-ML/result/2.ML.Rdata")
d<-c("SH3KBP1","PRNP","SV2B","PSMB9")##response 的signature
intersect(input,lasso_fea)
IDG_select<-IDG[IDG$Term %in% drug,]
LINC_select<-LINC[LINC$Term %in% drug,]
network<-c()
for(i in 1:dim(IDG_select)[1]){
  drug_tmp<-IDG_select[i,1]
  gene_tmp<-strsplit(IDG_select[i,9],";")[[1]]
  network_tmp<-as.data.frame(gene_tmp)
  network_tmp[,2]<-drug_tmp
  network_tmp[,3]<-"IDG"
  network<-rbind(network,network_tmp)
}
for(i in 1:dim(LINC_select)[1]){
  drug_tmp<-LINC_select[i,1]
  gene_tmp<-strsplit(LINC_select[i,9],";")[[1]]
  network_tmp<-as.data.frame(gene_tmp)
  network_tmp[,2]<-drug_tmp
  network_tmp[,3]<-"LINC"
  network<-rbind(network,network_tmp)
}
network<-unique(network)
names(network)<-c("gene","drug","from")
node1<-data.frame(node=unique(network$gene),attribute=1)
node2<-data.frame(node=unique(network$drug),attribute=2)
node<-rbind(node1,node2)
gene_laiyuan<-network[,c(1,3)]
gene_score<-pagerank_GS.score[pagerank_GS.score$GeneSymbol %in% node1$node,c(1,12)]
drug_score<-node2
names(drug_score)<-names(gene_score)
gene_score<-rbind(gene_score,drug_score)
write.csv(network,"/home/data/t060414/data/network-ML/result/4.drug/edges1.csv",row.names = F,quote=F)
write.csv(node,"/home/data/t060414/data/network-ML/result/4.drug/nodes1.csv",row.names = F,quote=F)
write.csv(gene_score,"/home/data/t060414/data/network-ML/result/4.drug/gene_score.csv",row.names = F,quote=F)
#erlotinib PMID 29579331 厄洛替尼治疗复发性或转移性皮肤鳞状细胞癌
#midostaurin PMID 16969355 多激酶抑制剂 midostaurin (PKC412A) 在转移性黑色素瘤中缺乏活性
#vandetanib PMID 29517106 Vandetanib 在体外和体内对未分化甲状腺癌具有抗肿瘤活性
#gefitinib PMID 21738104 A Phase II Study of Gefitinib in Patients with Metastatic Melanoma
#nelfinavir PMID 17283158 HIV protease inhibitor nelfinavir inhibits growth of human melanoma cells by induction of cell cycle arrest
#sorafenib  PMID 34691254 索拉非尼治疗具有 c-Kit 畸变的转移性黑色素瘤可减少肿瘤生长并促进生存
#vinorelbine 长春瑞滨完成了治疗转移性黑色素瘤的 2 期试验 https://go.drugbank.com/drugs/DB00361/clinical_trials?conditions=DBCOND0030058&phase=2&purpose=treatment&status=completed
#ixabepilone PMID 20098694 伊沙匹隆对于未接受过化疗（之前未经治疗）或之前接受过治疗的转移性黑色素瘤患者没有有意义的活性
#docetaxel PMID 8523052 多西紫杉醇对晚期恶性黑色素瘤具有活性。
#cabazitaxel PMID 26020806 Phase I dose-escalation study of cabazitaxel administered in combination with gemcitabine in patients with metastatic or unresectable advanced solid malignancies
#vinblastine PMID 8922197 采用顺铂、长春花碱和达卡巴嗪 (CVD) 联合化疗以及使用白细胞介素 2 和干扰素 α 的生物疗法治疗转移性黑色素瘤
#vincristine
