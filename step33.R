rm(list = ls())
pacman::p_load(tidyverse,export)
load('../data/IMMUNE_data/ALL_immu_genesets.Rdata')
load('../data/IMMUNE_data/TIDE.Rdata')
load('outdata/LUAD_TCGA_fpkm_symbol_maxgene.Rdata')
load('outdata/step4_group_data.rds')
fpkm[1:4,1:4]
fpkm <- fpkm[,sample_group$sample]
dim(fpkm)
# 1. immugene数据读取-----------------------------------------------------------------------
fpkm <- t(fpkm)
immugene |> duplicated() |> table();str(immugene)
# 94 个immugenes但是只有92个和fpkm匹配上
df_immune92 <- fpkm[, match(immugene,colnames(fpkm),nomatch = 0)] %>% as.data.frame()
is.na(df_immune92) %>% table()
df_immune92$group <- "ClusterA"
df_immune92$group[rownames(df_immune92) %in% allid[280:497]] <- "ClusterB"
table(df_immune92$group)

################ 2. 比较92个immugenes的差异###############

Innateimmune <- df_immune92
p_load(ggpubr,ggthemes)
bodat <- Innateimmune %>% 
  select(group, everything())
bodat2 <- bodat %>% 
  pivot_longer(cols = c(2:ncol(bodat)))
# 比较92个免疫相关分子的差异
x <- ggpubr::compare_means(value ~ group, group.by = "name",data = bodat2)
xnosig <- x[which(x$p > 0.05),] # p值＞0.05的选出来
nrow(x);nrow(xnosig)
# 92个分子，32个没有差异。
xsignif <- x[which(x$p < 0.05),] 

######### 3.61个免疫检查点的基因########

# 一共61个免疫检查点相关基因
immune_checkpoint <- fpkm[, match(x = c(immugeneINHI,immugeneSTIM),colnames(fpkm),nomatch = 0)] %>% as.data.frame()
str(immune_checkpoint)
setdiff(c(immugeneINHI,immugeneSTIM),colnames(immune_checkpoint)) # fpkm只有

Innateimmune <- immune_checkpoint
is.na(Innateimmune) %>% table()
Innateimmune$group <- "ClusterA"
Innateimmune$group[rownames(Innateimmune) %in% allid[280:497]] <- "ClusterB"
table(Innateimmune$group)

bodat <- Innateimmune %>% 
  select(group, everything())
bodat2 <- bodat %>% 
  pivot_longer(cols = c(2:ncol(Innateimmune)))
x <- compare_means(value ~ group, group.by = "name",data = bodat2)
xnosig <- x[which(x$p > 0.001),]
dim(x);dim(xnosig)
#### 选择p值小于0.001的基因做差异分析
bodatsig <- bodat2 %>%  
  filter(!name %in% xnosig$name)
ggboxplot(bodatsig,x="name",y="value",fill = "group",size = 0.1,
          # facet.by = "name",
          palette = c("#00CCFF","#FF3333"))+
  xlab("")+ylab("InfiltrationScore")+
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 0, vjust = 0))+
  # scale_y_continuous(limits = c(-1, 1),breaks = seq(-1,1,0.2))+
  stat_compare_means(aes(group=group),label = "p.signif",hide.ns = T,
                     bracket.size = 10)+
  theme_gray()+
  theme(axis.text = element_text(size = 8), 
        axis.text.x = element_text(size = 8, angle = 45, vjust = 0.6)) 

graph2ppt(file = "output/plots/IMMU.pptx", append = T)

# 免疫检查点基因有61个，但是不是每个都开发了药物。
# 下面几个是有药物的
# IDO-1,TIM-3(HAVCR2), LAG-3, KIR, GITR(TNFRSF18),
# VISTA, 4-1BB(TNFRSF9), OX40(TNFRSF4 ),
# B7-H3(CD276), CD27

cliblcgenes <-  c("PDCD1","CTLA4","CD274","LAG3","IDO1","HAVCR2",
                  "KIR",'TNFRSF18','VISTA','TNFRSF9','TNFRSF4','CD276','CD27')

bodatpd <- bodat2 %>% 
  filter(name %in% cliblcgenes)
ggboxplot(bodatpd,x="name",y="value",fill = "group",size = 0.1,
          palette = c("#00CCFF","#FF3333"))+
  xlab("")+ylab("log2(FPKM+1)")+
  theme(axis.text.x = element_text(size = 8, angle = 0, hjust = 0, vjust = 0))+
  # scale_y_continuous(limits = c(-1, 1),breaks = seq(-1,1,0.2))+
  stat_compare_means(aes(group=group),label = "p.signif",hide.ns = T,
                     bracket.size = 10)+
  theme_gray()+
  theme(axis.text = element_text(size = 8), 
        axis.text.x = element_text(size = 8, angle = 0, vjust = 0.6))
## 有差异的是 CD274,CTLA4,TIM-3(HAVCR2),CD27

ggboxplot(bodatpd[bodatpd$name %in%c('CD274','CTLA4','HAVCR2','CD27'),],
          x="name",y="value",fill = "group",size = 0.1,
          palette = c("#00CCFF","#FF3333"))+
  xlab("")+ylab("log2(FPKM+1)")+
  theme(axis.text.x = element_text(size = 8, angle = 0, hjust = 0, vjust = 0))+
  # scale_y_continuous(limits = c(-1, 1),breaks = seq(-1,1,0.2))+
  stat_compare_means(aes(group=group),label = "p.signif",hide.ns = T,
                     bracket.size = 10)+
  theme_gray()+
  theme(axis.text = element_text(size = 8), 
        axis.text.x = element_text(size = 8, angle = 0, vjust = 0.6))

graph2ppt(file = "output/plots/IMMU.pptx",app = T,width = 6, aspectr =1.5)


############## 4.热图 求分组均数到总中位数的距离###################
Innateimmune <- df_immune92
InnateimmuneA <- Innateimmune %>% filter(group == "ClusterA")
InnateimmuneA$group <- NULL
InnateimmuneA[] <- lapply(InnateimmuneA, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})
MeanA <- as.data.frame(apply(InnateimmuneA,2,mean))

InnateimmuneB <- Innateimmune%>%filter(group == "ClusterB")
InnateimmuneB$group <- NULL
InnateimmuneB[] <- lapply(InnateimmuneB, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})
MeanB <- as.data.frame(apply(InnateimmuneB,2,mean))

Mean <- cbind( MeanA, MeanB)
colnames(Mean)<-c("ClusterA","ClusterB")

Innateimmune$group <- NULL
Innateimmune[] <- lapply(Innateimmune, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})
Median <- as.data.frame(apply(Innateimmune,2,median))
Mean$ClusterA <- Mean$ClusterA - Median$`apply(Innateimmune, 2, median)`
Mean$ClusterB <- Mean$ClusterB - Median$`apply(Innateimmune, 2, median)`

# Mean2 <- cbind( MeanA, MeanB)
# colnames(Mean2)<-c("Mut","Wild")
# pheatmap(Mean2,scale = "column",cluster_rows = F, cluster_cols = F,
#          show_colnames =T,show_rownames = T, 
#          color = colorRampPalette(color.key)(50),
#          cutree_cols = 0,
#          cellwidth = 8,
#          cellheight = 8,)
#          
#z化
#Mean_Z=t(scale(t(Mean)))
#Mean_Z<-as
#Mean_Z<-as.data.frame(Mean_Z)
#Mean_Z[Mean_Z>  2]=2 #限定上限，使表达量大于0.05的等于0.05
#Mean_Z[Mean_Z< -2]= -2 #限定下限，使表达量小于-2的等于-2
#Mean_log<-10*Mean
# dd1 <- Mean
# dd2=apply(dd1,2,rescale )        ##归一化
# dd2=t(dd2)    
# pheatmap(dd2)


#热图绘制
color.key<-c("#3300CC","#3399FF","white","#FF3333","#CC0000")
library(pheatmap)
pheatmap(Mean,cluster_rows = F, cluster_cols = F,
         show_colnames =T,show_rownames = T, 
         color = colorRampPalette(color.key)(50),
         cutree_cols = 0,
         #annotation_colors = color.annotation,
         #annotation_col = annotation_col,
         cellwidth = 50,
         cellheight = 18,
         gaps_row = c(17,32,53),
        filename = "output/plots/免疫浸润热图.pdf",
         # silent = F
)
# graph2pdf(file = "output/plots/mianyijinrun.pdf",width = 20)
# 筛选出有意义的
xsignif$name
# 如果要横着放，就用t()转置一下
pheatmap(t(Mean[rownames(Mean) %in% xsignif$name,]),
         cluster_rows = F, cluster_cols = F,
         show_colnames =T,show_rownames = T, 
         color = colorRampPalette(color.key)(50),
         cutree_cols = 0,
         #annotation_colors = color.annotation,
         #annotation_col = annotation_col,
         cellwidth = 18,
         cellheight = 50,
         gaps_col =  c(11,23,34), 
         filename = "output/plots/免疫浸润热图_有差异_横.pdf"
         )

# 5.immuportal免疫相关基因集- -----------------------------------------------------------
df_immune_portal <- fpkm[, match(x = unique(immportgenes$Symbol),colnames(fpkm),nomatch = 0)] %>% as.data.frame()
Innateimmune <- df_immune_portal
is.na(Innateimmune) %>% table()
Innateimmune$group <- "ClusterA"
Innateimmune$group[rownames(Innateimmune) %in% allid[280:497]] <- "ClusterB"
table(Innateimmune$group)

# 画图数据准备
bodat <- Innateimmune %>% 
  select(group, everything())
bodat2 <- bodat %>% 
  pivot_longer(cols = c(2:ncol(Innateimmune)))
# 求p值
x <- compare_means(value ~ group, group.by = "name",data = bodat2)
xnosig <- x[which(x$p > 0.05),]
xhsig <- x[which(x$p < 0.001),]
xhhsig <- x[which(x$p.adj < 0.001),]
# 加上免疫功能
xhsig$func <- immportgenes$Category[match(x = xhsig$name, immportgenes$Symbol,nomatch = 0) ]
xhhsig$func <- immportgenes$Category[match(x = xhhsig$name, immportgenes$Symbol,nomatch = 0) ]
table(xhsig$func)
table(xhhsig$func)
table(immportgenes$Category)
# 排序，p值从小到大
xhhsig <- xhhsig[order(xhhsig$p.adj, decreasing = F),]

# quhua <- xhhsig %>% 
#   filter(func == "Chemokine_Receptors"|
#            func == "Chemokines") %>% 
#   select(name)
# bodatsig <- bodat2 %>% 
#   filter(name %in% quhua$name)
# 读取差异基因
all_diff <- read.csv(file = "outdata/mitocluster_mrna_all_diff.csv",row.names = 1)
all_diff$group = "not-significant"
#将adj.P.Val<0.05,logFC>0.5的基因设置为显著上调基因
#将adj.P.Val<0.05,logFC<0.5的基因设置为显著下调基因
all_diff$group[which((all_diff$adj.P.Val<0.05) & (all_diff$logFC> 1))] = "up-regulated"
all_diff$group[which((all_diff$adj.P.Val<0.05) & (all_diff$logFC< -1))] = "down-regulated"
#查看上调和下调基因的数目
table(all_diff$group)
all_diff$group<-as.factor(all_diff$group)
all_diff$logP<- -log10(all_diff$adj.P.Val)
### 选差异top10 
top10 <- intersect(xhhsig$name,all_diff$ID[all_diff$group != "not-significant"])[1:10]
top10
# top10 基因功能
xhhsig$func[which(xhhsig$name %in% top10)]
bodattop <- bodat2[bodat2$name %in% top10,]

ggboxplot(bodattop,x="name",y="value",fill = "group",size = 0.1,
          palette = c("#00CCFF","#FF3333"))+
  xlab("")+ylab("log2(FPKM+1)")+
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 0, vjust = 0))+
  stat_compare_means(aes(group=group),label = "p.signif",hide.ns = T,
                     bracket.size = 10)+
  theme_gray()+
  theme(axis.text = element_text(size = 8), 
        axis.text.x = element_text(size = 8, angle = 45, vjust = 0.6))

graph2ppt(file = "output/plots/IMMU.pptx", append= T)

bodatpd <- bodat2 %>% 
  filter(name %in% c("PDCD1","CTLA4","CD274","TGFB1","VEGFA","VEGFB"))
ggboxplot(bodatpd,x="name",y="value",fill = "group",size = 0.1,
          palette = c("#00CCFF","#FF3333"))+
  xlab("")+ylab("log2(FPKM+1)")+
  theme(axis.text.x = element_text(size = 8, angle = 0, hjust = 0, vjust = 0))+
  # scale_y_continuous(limits = c(-1, 1),breaks = seq(-1,1,0.2))+
  stat_compare_means(aes(group=group),label = "p.signif",hide.ns = T,
                     bracket.size = 10)+
  theme_gray()+
  theme(axis.text = element_text(size = 8), 
        axis.text.x = element_text(size = 8, angle = 0, vjust = 0.6))

############### 7.相关性图############

aidat <- read.table(file = "cibersort/immucellAI/ImmuCellAI_icb_result.txt",header = T)
aidat <- aidat %>% 
  select(Th17,nTreg, InfiltrationScore)
TCGAtmb<-read.csv("outdata/TCGA_LUAD_tmv.csv", row.names = 1)

tmb <- TCGAtmb %>% 
  filter(Tumor_Sample_Barcode %in% allid) %>% 
  select(Tumor_Sample_Barcode, total) %>% 
  rename(SampleID = Tumor_Sample_Barcode,
         tmb =  total)
rownames(tmb) <- tmb[,1]  
tmb <- tmb[allid, -1] 
repla <- mean(tmb,na.rm = T)
tmb <- replace_na(data = tmb,replace = repla)

corgene <- c("PDCD1","CTLA4","CD274","TGFB1","VEGFA","VEGFB", 
             "HLA-DQB2","CD40","CX3CL1","ICAM1","ITGB2")
cordat <- mrna[, match(corgene,colnames(mrna),nomatch = 0)]
cordat <- cbind(aidat,tmb,cordat)

cormat <- round(cor(cordat), 2) 

ggcorrplot(cormat, type = "lower",
           outline.col = "white",
           ggtheme = ggplot2::theme_bw(),
)
graph2ppt(file = "output/plots/6免疫基因.pptx",append = T)
### 样式复杂的相关性图
p_load(PerformanceAnalytics)
my_data <- cordat
chart.Correlation(my_data, histogram=TRUE, pch=19)


