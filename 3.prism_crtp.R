# 确定潜在药物靶标 结果不好
##TCGA 3期
load("/home/data/t060414/data/multiomics-age-related/data/TCGA-SKCM/data_exprs_FPKM.RData")
load("/home/data/t060414/data/multiomics-age-related/data/TCGA-SKCM/data_clinical.RData")
# Load data
#stage_score_TCGA<-stage_score_TCGA[Cor_stage_score_all$pathway,]
gene_expr <- as.data.frame(data_exprs_FPKM)
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
clinical$pathologic_stage<-as.numeric(clinical$pathologic_stage)
row.names(clinical)<-clinical$patient
sample<-intersect(names(gene_expr),row.names(clinical))
clinical<-clinical[sample,]
gene_expr<-gene_expr[,sample]
mydata<-gene_expr
clinical$group <- ifelse(as.numeric(clinical$pathologic_stage)<3,"early","lately")
ctrp.pred.auc<-read.csv("/home/data/t060414/data/network-ML/result/4.drug/SKCM_ctrp.pred.auc.csv",row.names = 1)
names(ctrp.pred.auc)<-substr(names(ctrp.pred.auc),1,12)
names(ctrp.pred.auc)<-gsub("[.]","-",names(ctrp.pred.auc))
prism.pred.auc<-read.csv("/home/data/t060414/data/network-ML/result/4.drug/SKCM_prism.pred.auc.csv",row.names = 1)
names(prism.pred.auc)<-substr(names(prism.pred.auc),1,12)
names(prism.pred.auc)<-gsub("[.]","-",names(prism.pred.auc))
#

top.pps <-clinical[clinical$group=="early", 1]
bot.pps <-clinical[clinical$group=="lately", 1]


## 1.差异药敏分析

ctrp.log2fc <- c()
for (i in 1:nrow(ctrp.pred.auc)) {
  #display.progress(index = i,totalN = nrow(ctrp.pred.auc))
  d <- rownames(ctrp.pred.auc)[i]
  a <- mean(as.numeric(ctrp.pred.auc[d,intersect(top.pps,names(ctrp.pred.auc))])) # 上十分位数的AUC均值
  b <- mean(as.numeric(ctrp.pred.auc[d,intersect(bot.pps,names(ctrp.pred.auc))])) # 下十分位数的AUC均值
  fc <- b/a
  log2fc <- log2(fc); names(log2fc) <- d
  ctrp.log2fc <- c(ctrp.log2fc,log2fc)
}
candidate.ctrp <- ctrp.log2fc[abs(ctrp.log2fc)> 0.1] # 这里我调整了阈值，控制结果数目

prism.log2fc <- c()
for (i in 1:nrow(prism.pred.auc)) {
 # display.progress(index = i,totalN = nrow(prism.pred.auc))
  d <- rownames(prism.pred.auc)[i]
  a <- mean(as.numeric(prism.pred.auc[d,intersect(top.pps,names(prism.pred.auc))])) # 上十分位数的AUC均值
  b <- mean(as.numeric(prism.pred.auc[d,intersect(bot.pps,names(prism.pred.auc))])) # 下十分位数的AUC均值
  fc <- b/a
  log2fc <- log2(fc); names(log2fc) <- d
  prism.log2fc <- c(prism.log2fc,log2fc)
}
candidate.prism <- prism.log2fc[abs(prism.log2fc) > 0.2] # 这里我调整了阈值，控制结果数目


## 2.Spearman相关性分析，用于绘制左图

```{r}
ctrp.cor <- ctrp.cor.p <- c()
for (i in 1:nrow(ctrp.pred.auc)) {
  display.progress(index = i,totalN = nrow(ctrp.pred.auc))
  d <- rownames(ctrp.pred.auc)[i]
  a <- as.numeric(ctrp.pred.auc[d,rownames(Sinfo)]) 
  b <- as.numeric(Sinfo$PPS)
  r <- cor.test(a,b,method = "spearman")$estimate; names(r) <- d
  p <- cor.test(a,b,method = "spearman")$p.value; names(p) <- d
  ctrp.cor <- c(ctrp.cor,r)
  ctrp.cor.p <- c(ctrp.cor.p,p)
}
candidate.ctrp2 <- ctrp.cor[ctrp.cor < -0.4]  # 这里我调整了阈值，控制结果数目
ctrp.candidate <- intersect(names(candidate.ctrp),names(candidate.ctrp2))

prism.cor <- prism.cor.p <- c()
for (i in 1:nrow(prism.pred.auc)) {
  display.progress(index = i,totalN = nrow(prism.pred.auc))
  d <- rownames(prism.pred.auc)[i]
  a <- as.numeric(prism.pred.auc[d,rownames(Sinfo)]) 
  b <- as.numeric(Sinfo$PPS)
  r <- cor.test(a,b,method = "spearman")$estimate; names(r) <- d
  p <- cor.test(a,b,method = "spearman")$p.value; names(p) <- d
  prism.cor <- c(prism.cor,r)
  prism.cor.p <- c(prism.cor.p,p)
}
candidate.prism2 <- prism.cor[prism.cor < -0.35]  
prism.candidate <- intersect(names(candidate.prism),names(candidate.prism2))
```

# 开始画图

## 1. 左侧相关性图

```{r}
# 设置颜色
darkblue <- "#0772B9"
lightblue <- "#48C8EF"

cor.data <- data.frame(drug = ctrp.candidate,
                       r = ctrp.cor[ctrp.candidate],
                       p = -log10(ctrp.cor.p[ctrp.candidate]))
p1 <- ggplot(data = cor.data,aes(r,forcats::fct_reorder(drug,r,.desc = T))) +
  geom_segment(aes(xend=0,yend=drug),linetype = 2) +
  geom_point(aes(size=p),col = darkblue) +
  scale_size_continuous(range =c(2,8)) +
  scale_x_reverse(breaks = c(0, -0.3, -0.5),
                  expand = expansion(mult = c(0.01,.1))) + #左右留空
  theme_classic() +
  labs(x = "Correlation coefficient", y = "", size = bquote("-log"[10]~"("~italic(P)~"-value)")) + 
  theme(legend.position = "bottom", 
        axis.line.y = element_blank())

cor.data <- data.frame(drug = prism.candidate,
                       r = prism.cor[prism.candidate],
                       p = -log10(prism.cor.p[prism.candidate]))
cor.data$drug <- sapply(strsplit(cor.data$drug," (",fixed = T), "[",1)

p2 <- ggplot(data = cor.data,aes(r,forcats::fct_reorder(drug,r,.desc = T))) +
  geom_segment(aes(xend=0,yend=drug),linetype = 2) +
  geom_point(aes(size=p),col = darkblue) +
  scale_size_continuous(range =c(2,8)) +
  scale_x_reverse(breaks = c(0, -0.3, -0.5),
                  expand = expansion(mult = c(0.01,.1))) + #左右留空
  theme_classic() +
  labs(x = "Correlation coefficient", y = "", size = bquote("-log"[10]~"("~italic(P)~"-value)")) + 
  theme(legend.position = "bottom", 
        axis.line.y = element_blank())
```

## 2.右侧箱型图

```{r, fig.width=8, fig.height=8}
ctrp.boxdata <- NULL
for (d in ctrp.candidate) {
  a <- as.numeric(ctrp.pred.auc[d,rownames(top.pps)]) 
  b <- as.numeric(ctrp.pred.auc[d,rownames(bot.pps)])
  p <- wilcox.test(a,b)$p.value
  s <- as.character(cut(p,c(0,0.001,0.01,0.05,1),labels = c("***","**","*","")))
  ctrp.boxdata <- rbind.data.frame(ctrp.boxdata,
                                   data.frame(drug = d,
                                              auc = c(a,b),
                                              p = p,
                                              s = s,
                                              group = rep(c("High PPS","Low PPS"),c(nrow(top.pps),nrow(bot.pps))),
                                              stringsAsFactors = F),
                                   stringsAsFactors = F)
}
p3 <- ggplot(ctrp.boxdata, aes(drug, auc, fill=group)) + 
  geom_boxplot(aes(col = group),outlier.shape = NA) + 
  # geom_text(aes(drug, y=min(auc) * 1.1, 
  #               label=paste("p=",formatC(p,format = "e",digits = 1))),
  #           data=ctrp.boxdata, 
  #           inherit.aes=F) + 
  geom_text(aes(drug, y=max(auc)), 
            label=ctrp.boxdata$s,
            data=ctrp.boxdata, 
            inherit.aes=F) + 
  scale_fill_manual(values = c(darkblue, lightblue)) + 
  scale_color_manual(values = c(darkblue, lightblue)) + 
  xlab(NULL) + ylab("Estimated AUC value") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5,vjust = 0.5,size = 10),
        legend.position = "bottom",
        legend.title = element_blank()) 
dat <- ggplot_build(p3)$data[[1]]

p3 <- p3 + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", inherit.aes = F)

prism.boxdata <- NULL
for (d in prism.candidate) {
  a <- as.numeric(prism.pred.auc[d,rownames(top.pps)]) 
  b <- as.numeric(prism.pred.auc[d,rownames(bot.pps)])
  p <- wilcox.test(a,b)$p.value
  s <- as.character(cut(p,c(0,0.001,0.01,0.05,1),labels = c("***","**","*","")))
  prism.boxdata <- rbind.data.frame(prism.boxdata,
                                    data.frame(drug = d,
                                               auc = c(a,b),
                                               p = p,
                                               s = s,
                                               group = rep(c("High PPS","Low PPS"),c(nrow(top.pps),nrow(bot.pps))),
                                               stringsAsFactors = F),
                                    stringsAsFactors = F)
}
prism.boxdata$drug <- sapply(strsplit(prism.boxdata$drug," (",fixed = T), "[",1)

p4 <- ggplot(prism.boxdata, aes(drug, auc, fill=group)) + 
  geom_boxplot(aes(col = group),outlier.shape = NA) + 
  # geom_text(aes(drug, y=min(auc) * 1.1, 
  #               label=paste("p=",formatC(p,format = "e",digits = 1))),
  #           data=prism.boxdata, 
  #           inherit.aes=F) + 
  geom_text(aes(drug, y=max(auc)), 
            label=prism.boxdata$s,
            data=prism.boxdata, 
            inherit.aes=F) + 
  scale_fill_manual(values = c(darkblue, lightblue)) + 
  scale_color_manual(values = c(darkblue, lightblue)) + 
  xlab(NULL) + ylab("Estimated AUC value") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5,vjust = 0.5,size = 10),
        legend.position = "bottom",
        legend.title = element_blank())
dat <- ggplot_build(p4)$data[[1]]

p4 <- p4 + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", inherit.aes = F)
```

## 3. 合并图像

```{r, fig.width=8, fig.height=8}
plot_grid(p1, p3, p2, p4, labels=c("A", "", "B", ""), 
          ncol=2, 
          rel_widths = c(2, 2)) #左右两列的宽度比例
ggsave(filename = "drug target.pdf",width = 8,height = 8)

# 保存镜像
#save.image("drug.RData")
```

# Session Info

```{r}
sessionInfo()
```