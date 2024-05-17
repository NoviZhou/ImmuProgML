rm(list=ls())
###DNB pathway
#load("/home/data/t060414/data/network-ML/result/pathway_cor_stage.Rdata")
load("/home/data/t060414/data/network-ML/result/stage_score.Rdata")
load("/home/data/t060414/data/multiomics-age-related/data/TCGA-SKCM/data_clinical.RData")
# Load data
#stage_score_TCGA<-stage_score_TCGA[Cor_stage_score_all$pathway,]
gene_expr <- as.data.frame(stage_score_TCGA)
names(gene_expr)<-substr(names(gene_expr),1,12)
clinical<-unique(data.frame(patient=data_clinical$bcr_patient_barcode,
                            pathologic_stage=data_clinical$stage_event_pathologic_stage))

clinical$pathologic_stage<-gsub("Stage IV",4,clinical$pathologic_stage)
clinical$pathologic_stage<-gsub("Stage 0",0,clinical$pathologic_stage)
clinical$pathologic_stage<-gsub("Stage IIIC",3,clinical$pathologic_stage)
clinical$pathologic_stage<-gsub("Stage IIIB",3,clinical$pathologic_stage)
clinical$pathologic_stage<-gsub("Stage IIIA",3,clinical$pathologic_stage)
clinical$pathologic_stage<-gsub("Stage III",3,clinical$pathologic_stage)
clinical$pathologic_stage<-gsub("Stage IIC",2,clinical$pathologic_stage)
clinical$pathologic_stage<-gsub("Stage IIB",2,clinical$pathologic_stage)
clinical$pathologic_stage<-gsub("Stage IIA",2,clinical$pathologic_stage)
clinical$pathologic_stage<-gsub("Stage II",2,clinical$pathologic_stage)

clinical$pathologic_stage<-gsub("Stage IA",1,clinical$pathologic_stage)
clinical$pathologic_stage<-gsub("Stage IB",1,clinical$pathologic_stage)
clinical$pathologic_stage<-gsub("Stage I",1,clinical$pathologic_stage)

clinical$pathologic_stage<-gsub("I/II NOS",1,clinical$pathologic_stage)
#clinical<-clinical[-which(clinical$pathologic_stage=="I/II NOS"),]
clinical<-clinical[-which(clinical$pathologic_stage==""),]
row.names(clinical)<-clinical$patient
sample<-intersect(names(gene_expr),row.names(clinical))
clinical<-clinical[sample,]
gene_expr<-gene_expr[,sample]
library(DNB)
# New DNB object
# all-zero genes at any time point will be removed
dnb <- new_DNB(data = gene_expr,
               time = as.factor(clinical$pathologic_stage))

# Calulate correlation and coefficient of variation
dnb <- cal_cor(dnb)
dnb <- cal_cv(dnb)

# Search DNB
dnb <- search_candidates(dnb, min_size =20, max_size =Inf)
dnb <- cal_final(dnb)

# Get results
dnb_genes <- get_DNB_genes(dnb)
candidates <- get_candidates(dnb)
plot_DNB(candidates)
final <- get_final(dnb)
pdf("/home/data/t060414/data/network-ML/result/1.DNB.TCGA.pdf",width = 8,height =6)
plot_DNB(final)
dev.off()
gene<-get_DNB_genes(dnb)
dnb_TCGA<-dnb
pathway_TCGA<-gene

####GEO数据
#melanoma_pdata_list[["GSE19234"]] #转移性  "IIIB" "IIIC" "IV"   "III nA"
#melanoma_pdata_list[["GSE98394"]] #未治疗
#melanoma_pdata_list[["GSE22153"]] #转移 (只有转移) "IV"  "III"
#melanoma_pdata_list[["GSE22154"]] #转移 (only 转移)  "IV"
#melanoma_pdata_list[["GSE54467"]] #转移
load("/home/data/t060414/data/network-ML/result/stage_score.Rdata")
load("/home/data/t060414/data/ActivePathway/data/melanoma_pdata_list.Rdata")
#GSE98394
gene_expr <- as.data.frame(stage_score$GSE98394)
clinical<-melanoma_pdata_list$GSE98394[,c("GEO_ID","Stage_ajcc")]
names(clinical)<-c("patient","pathologic_stage")
unique(clinical$pathologic_stage)

clinical$pathologic_stage<-gsub("IIIc",8,clinical$pathologic_stage)
clinical$pathologic_stage<-gsub("IIIb",7,clinical$pathologic_stage)
clinical$pathologic_stage<-gsub("IIIa",6,clinical$pathologic_stage)
clinical$pathologic_stage<-gsub("III",3,clinical$pathologic_stage)
clinical$pathologic_stage<-gsub("IIc",5,clinical$pathologic_stage)
clinical$pathologic_stage<-gsub("IIb",4,clinical$pathologic_stage)
clinical$pathologic_stage<-gsub("IIa",3,clinical$pathologic_stage)
clinical$pathologic_stage<-gsub("Ia",1,clinical$pathologic_stage)
clinical$pathologic_stage<-gsub("Ib",2,clinical$pathologic_stage)
clinical<-clinical[-which(clinical$pathologic_stage=="na"),]
row.names(clinical)<-clinical$patient

sample<-intersect(names(gene_expr),row.names(clinical))
clinical<-clinical[sample,]
gene_expr<-gene_expr[,sample]
library(DNB)
# New DNB object
# all-zero genes at any time point will be removed
dnb <- new_DNB(data = gene_expr,
               time = as.factor(clinical$pathologic_stage))

# Calulate correlation and coefficient of variation
dnb <- cal_cor(dnb)
dnb <- cal_cv(dnb)

# Search DNB
dnb <- search_candidates(dnb, min_size =20, max_size =Inf)
dnb <- cal_final(dnb)

# Get results
dnb_genes <- get_DNB_genes(dnb)
candidates <- get_candidates(dnb)
plot_DNB(candidates)
final <- get_final(dnb)
pdf("/home/data/t060414/data/network-ML/result/1.DNB.GSE98394.pdf",width = 8,height =6)
plot_DNB(final)
dev.off()
gene<-get_DNB_genes(dnb)
dnb_GSE98394<-dnb
pathway_GSE98394<-gene

###
#GSE190113
load("/home/data/t060414/data/ActivePathway/data/GSE190113_pdata.Rdata")
gene_expr <- as.data.frame(stage_score_GSE190113)
clinical<-GSE190113_pdata[,c("Sample_title","Stage")]
names(clinical)<-c("patient","pathologic_stage")
unique(clinical$pathologic_stage)

clinical$pathologic_stage<-gsub("stage: 3",3,clinical$pathologic_stage)
clinical$pathologic_stage<-gsub("stage: 2",2,clinical$pathologic_stage)
clinical$pathologic_stage<-gsub( "stage: 4" ,4,clinical$pathologic_stage)
clinical$pathologic_stage<-gsub("stage: 1",1,clinical$pathologic_stage)
clinical$patient<-gsub("-",".",clinical$patient)
#clinical<-clinical[-which(clinical$pathologic_stage==" Stage Unknown"),]
row.names(clinical)<-clinical$patient
clinical$pathologic_stage<-as.numeric(clinical$pathologic_stage)

sample<-intersect(names(gene_expr),row.names(clinical))
clinical<-clinical[sample,]
gene_expr<-gene_expr[,sample]
library(DNB)
# New DNB object
# all-zero genes at any time point will be removed
dnb <- new_DNB(data = gene_expr,
               time = as.factor(clinical$pathologic_stage))

# Calulate correlation and coefficient of variation
dnb <- cal_cor(dnb)
dnb <- cal_cv(dnb)

# Search DNB
dnb <- search_candidates(dnb, min_size =40, max_size =Inf)
dnb <- cal_final(dnb)

# Get results
dnb_genes <- get_DNB_genes(dnb)
candidates <- get_candidates(dnb)
plot_DNB(candidates)
final <- get_final(dnb)
pdf("/home/data/t060414/data/network-ML/result/1.DNB.GSE190113.pdf",width = 8,height =6)
plot_DNB(final)
dev.off()
gene<-get_DNB_genes(dnb)
dnb_GSE190113<-dnb
pathway_GSE190113<-gene

pathway_select<-unique(c(intersect(pathway_GSE190113,pathway_GSE98394),intersect(pathway_GSE190113,pathway_TCGA),intersect(pathway_TCGA,pathway_GSE98394)))
save(pathway_select,pathway_GSE190113,dnb_GSE190113,pathway_GSE98394,dnb_GSE98394,pathway_TCGA,dnb_TCGA,file="/home/data/t060414/data/network-ML/result/DNB_pathway.Rdata")
####未用GSE54467
gene_expr <- as.data.frame(stage_score$GSE54467)
clinical<-melanoma_pdata_list$GSE54467[,c("GEO_ID","Stage_clinical")]
names(clinical)<-c("patient","pathologic_stage")
unique(clinical$pathologic_stage)

clinical$pathologic_stage<-gsub(" Stage III",3,clinical$pathologic_stage)
clinical$pathologic_stage<-gsub("Stage II",2,clinical$pathologic_stage)
clinical$pathologic_stage<-gsub( " Stage II" ,2,clinical$pathologic_stage)
clinical$pathologic_stage<-gsub(" Stage I",1,clinical$pathologic_stage)

clinical<-clinical[-which(clinical$pathologic_stage==" Stage Unknown"),]
row.names(clinical)<-clinical$patient
clinical$pathologic_stage<-as.numeric(clinical$pathologic_stage)
sample<-intersect(names(gene_expr),row.names(clinical))
clinical<-clinical[sample,]
gene_expr<-gene_expr[,sample]
library(DNB)
# New DNB object
# all-zero genes at any time point will be removed
dnb <- new_DNB(data = gene_expr,
               time = as.factor(clinical$pathologic_stage))

# Calulate correlation and coefficient of variation
dnb <- cal_cor(dnb)
dnb <- cal_cv(dnb)

# Search DNB
dnb <- search_candidates(dnb, min_size =40, max_size =Inf)
dnb <- cal_final(dnb)

# Get results
dnb_genes <- get_DNB_genes(dnb)
candidates <- get_candidates(dnb)
plot_DNB(candidates)
final <- get_final(dnb)
plot_DNB(final)
gene<-get_DNB_genes(dnb)
