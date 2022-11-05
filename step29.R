rm(list = ls())
load("outdata/step10_GSVA24cell.rds")
p_load(limma,ggplot2,pheatmap,export)
############ 2 limma差异分析##########

pheatmap(mat = as.matrix(allgsva),show_colnames = F,cluster_cols = F)
stopifnot(allid497 == colnames(allgsva))


######## 2.1 limma分组信息#######
dt <- allgsva
group_list <- c(rep(0,279), rep(1, 218))
group <- factor(group_list, 
                levels = c(0,1), 
                labels = c("ClusterA","ClusterB"))
group
design <- model.matrix(~0+group)
rownames(design) <- colnames(dt)
colnames(design) <- levels(group)
design
contrast.matrix <- makeContrasts(ClusterA - ClusterB,levels = design)
contrast.matrix
###############  2.2 差异分析###########
rt <- as.matrix(dt)
fit <- lmFit(rt,design)
fit1 <- contrasts.fit(fit, contrast.matrix)
fit1 <- eBayes(fit1)
qqt(fit1$t, df = fit1$df.prior+fit1$df.residual, pch = 16, cex = 0.2)
abline(0,1) #QQ图看正态分布?
all_diff <- topTable(fit1, adjust.method = 'fdr',coef=1,p.value = 1,lfc <- log(1,2),number = Inf,sort.by = 'logFC')
head(all_diff)
# 所有差异加一列ID改成gene名字
all_diff$ID <- rownames(all_diff)
# 加一列表示logP
all_diff$logP<- -log10(all_diff$adj.P.Val)
p_load(ggpubr,ggthemes)


ggscatter(all_diff,x="logFC",y="logP")+ theme_base()

##美化火山图tsv
all_diff$Group = "not-significant"
#将adj.P.Val<0.05,logFC>0.5的基因设置为显著上调基因
#将adj.P.Val<0.05,logFC<0.5的基因设置为显著下调基因
all_diff$Group[which((all_diff$adj.P.Val<0.05) & (all_diff$logFC> 0.5))] = "up-regulated"
all_diff$Group[which((all_diff$adj.P.Val<0.05) & (all_diff$logFC< -0.5))] = "down-regulated"
#查看上调和下调基因的数目
table(all_diff$Group)
all_diff$Group<-as.factor(all_diff$Group)
all_diff$logP<- -log10(all_diff$adj.P.Val)
#火山图，加上颜色
ggscatter(all_diff,x="logFC",y="logP",
          color = "Group",
          palette = c("green","gray","red"),
          size = 1.5)+theme_base()
#再加上辅助线
p <- ggscatter(all_diff,x="logFC",y="logP",
               color = "Group",
               palette = c("green","gray","red"),
               size = 1.5)+theme_base() + 
  geom_hline(yintercept = 1.30,linetype = "dashed") + 
  geom_vline(xintercept = c(-0.5,0.5),linetype = "dashed")
p
graph2ppt(file = "output/plots/GSVA.pptx",append = T)

# 火山图没法看

DEG <- all_diff
DEG <- DEG[,c(1,2,6)] #要求DEG有4列，geneid、FC、regulation、pval
DEG <- DEG[is.finite(DEG$logFC),]#去掉无限小数的行
DEG <- DEG[ DEG$adj.P.Val<0.05,] # 筛选logFC大于1.5的
names(DEG) <- c("genes","fold_change","p_value")
DEG$regulation <- "up"
DEG$regulation[DEG$fold_change<0] <- "down"

library(scales)

dd <- rt # dd表示总的表达矩阵
DEG1<- DEG[order(abs(DEG$fold_change),decreasing = T),]#按照FC排序
dd1=dd[rownames(dd) %in% DEG1$genes,]# 在表达矩阵中选取这些基因
dd2=apply(dd1,1,rescale)         ##归一化
dd2=t(dd2)                     ##转置回来

# 2.3 热图加上临床注释信息 ----------------------------------------------------------
load(file = "outdata/LUAD_cln_clean.rds")
load(file = "outdata/step4_group_data.rds")
phe3 <- phe_f
rownames(phe3) <-  NULL
phe3 <- phe3[match(x = allid497, table = phe3$sample), ]
rownames(phe3) <- phe3$sample
phe3$group <- sample_group$group[match(x = phe3$sample, table = sample_group$sample)]
phe3_anno <- phe3 %>% dplyr::select(group, age, gender, smoke_group,ajcc_T, ajcc_N, radiotherapy)
pheatmap(dd2,cluster_rows = T,cluster_cols = F,show_colnames = F,annotation_col = phe3_anno)
graph2ppt(file = "output/plots/GSVA.pptx",append = T)



############## 3.1 箱线图绘制###############

library(reshape2)
allgsva1 <- t(allgsva)
allgsva1 <- as.data.frame(allgsva1)
allgsva1$group <- sample_group$group[match(x = rownames(allgsva1), 
                                           table = sample_group$sample)]
allgsva1 <- reshape2::melt(allgsva1,id.vars=c("group"))
allgsva1$group <- as.factor(allgsva1$group)

#箱线图绘制
library(tidyverse)
a = allgsva1 %>% 
  group_by(variable) %>% 
  summarise(m = median(value)) %>% 
  arrange(desc(m)) %>% 
  pull(variable)
allgsva1$variable = factor(allgsva1$variable,levels = a)
allgsva1$value <- as.numeric(allgsva1$value)
# 两个组之间的免疫浸润整体做差异分析
compare_means(value ~ group, data = allgsva1)
p <- ggboxplot(allgsva1,x="variable",y="value",fill = "group",size = 0.1,
               palette = c("#00CCFF","#FF3333"))+
  xlab("")+ylab("GSVA_ES")+
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 0, vjust = 0))+
  scale_y_continuous(limits = c(-1, 1),breaks = seq(-1,1,0.2))+
  stat_compare_means(aes(group=group),label = "p.signif",hide.ns = T,
                     bracket.size = 10)+
  theme_gray()
p 
p+theme(axis.text = element_text(size = 8), 
        axis.text.x = element_text(size = 8, angle = 45, vjust = 0.6)) 
graph2ppt(file= "output/plots/GSVA.pptx", append=T,asp=2,width = 20)
## 只要有差异的细胞的图
allgsva1sig <- allgsva1 %>% 
  filter(!variable %in% c("Tgd","Cytotoxic.cells","Th2.cells",
                          "Mast.cells","Neutrophils","Eosinophils",
                          "NK.CD56dim.cells"))
ggboxplot(allgsva1sig,x="variable",y="value",fill = "group",size = 0.1,
          palette = c("#00CCFF","#FF3333"))+
  xlab("")+ylab("GSVA_ES")+
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 0, vjust = 0))+
  scale_y_continuous(limits = c(-1, 1),breaks = seq(-1,1,0.2))+
  stat_compare_means(aes(group=group),label = "p.signif",hide.ns = T,
                     bracket.size = 10)+
  theme_gray()+
  theme(axis.text = element_text(size = 8), 
        axis.text.x = element_text(size = 8, angle = 45, vjust = 0.6)) 
library(export)
graph2ppt(file= "output/plots/GSVA.pptx",append= T,asp= 2,width=20)

# 4. CD8/TREG比值 --------------------------------------------------------------
allgsva2 <- t(allgsva)
allgsva2 <- as.data.frame(allgsva2)
allgsva2$group <- sample_group$group[match(x = rownames(allgsva2), 
                                           table = sample_group$sample)]

ratioofcd8treg <- allgsva2 %>% 
  select(CD8.T.cells,Treg,group) %>% 
  mutate(ratio = abs( CD8.T.cells / Treg ))
# ratioofcd8treg <- pivot_longer(allgsva2,cols = ratio)

ggviolin(x="group",y="ratio",fill = "group",data = ratioofcd8treg,size = 0.1,
         palette = c("#00CCFF","#FF3333"),add = "boxplot", 
         add.params = list(fill = "white"))+
  xlab("")+
  ylab("")+
  scale_y_continuous(limits =c(0, 1.5),breaks = seq(0,1.6,0.2))+ 
  stat_compare_means(aes(group=group),label = "p.format",hide.ns = F,
                     bracket.size = 20,label.x = 1.35,label.y = 1.3,size = 10)+
  theme_minimal()+
  theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5))+
  labs(title = "CD8+ Cell/ Treg")+
  theme(plot.title = element_text(size = 24, color = "black", 
                                  hjust = 0.5, 
                                  angle = 0))


##这里有问题，就是明明数据都是大于0的，但是小提琴图会有＜0的图出现
## 原因是小提琴图中间是有计算过程的
##y轴是概率密度，以数据为基础计算而来。
## 虽然数据都是大于0的，而且有限个数的，但是密度曲线在x轴上会超出数据的上下限。
## 解决方法：使用boxplot

ggboxplot(ratioofcd8treg,x="group",y="ratio",fill = "group",size = 0.1,
          palette = c("#00CCFF","#FF3333"))+
  xlab("")+ylab("")+
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 0, vjust = 0))+
  scale_y_continuous(limits = c(0, 1),breaks = seq(0,1,0.2))+
  stat_compare_means(aes(group=group),label = "p.format",hide.ns = F,
                     bracket.size = 10)+
  theme_minimal()+
  theme(axis.text = element_text(size = 8), 
        axis.text.x = element_text(size = 8, angle = 0, vjust = 0.6))
graph2ppt(file = "output/plots/GSVA.pptx", append = T,asp=1,width = 10)
graph2pdf(file = "output/plots/CD8+与Treg比值")
library(ggplot2)
ggplot(ratioofcd8treg,aes(x=group,y= ratio,fill = group))+
  geom_boxplot(aes(fill=group,),
               notch=TRUE,outlier.colour="red", outlier.shape=7,outlier.size=4)+
  scale_y_continuous(limits = c(0, 1),breaks = seq(0,1,0.2))+
  scale_fill_discrete(c("#00CCFF","#FF3333"))+
  # stat_summary(fun="mean",geom="point",shape=23,size=3,fill="white")+# 加个平均值
  theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.line = element_line(colour = "black"))#使背景为空白并保留坐标轴为黑色
library(ggpubr)
ggbarplot(tab,x="Group",y= "Number",fill="ICB Response", size = 0.1,
          palette = c("#00CCFF","#FF3333"))+
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 0, vjust = 0))+
  theme_cleveland()+
  theme(axis.text = element_text(size = 8), 
        axis.text.x = element_text(size = 8, angle = 45, vjust = 0.6))


ggbarplot(tab,x="Group",y= "Number",fill="ICB Response", size = 0.1,
          palette = c("#00CCFF","#FF3333"),legend = "right")+
  xlab("")+ylab("")+
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 0, vjust = 0))+
  stat_compare_means(aes(group=Group),label = "p.format",hide.ns = F,method = "chisq.test",
                     bracket.size = 10,)+
  theme_cleveland()+
  theme(axis.text = element_text(size = 8), 
        axis.text.x = element_text(size = 8, angle = 0, vjust = 0.6))
graph2ppt(file = "output/plots/immuneAI.pptx", append = T)

kafang <- matrix(data = c(83,241,28,141),nrow = 2,byrow = T)
chisq.test(kafang)
