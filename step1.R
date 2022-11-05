
# 数据准备 --------------------------------------------------------------------
p_load(tidyverse,ggplot2,GSVA)

################ GSVA 分析 ##########
##### 注意，TCGAgs必须是的变量必须是因子######
# 因子很容易转换成任何其他数据类型
# 但是其他的数据类型，比如字符型无法直接转换成因子

fusion <- data.table::fread("data/mitofusiongeneset.gmt",
                   header = F,stringsAsFactors = T,
                   ) %>% .[-1-2]
autophagy <- data.table::fread("data/mitoautophagy.gmt",
                   header = F,stringsAsFactors = T,
                   )
fusion_reg <- data.table::fread("data/mitofusionregulationgeneset.gmt",
                   header = F,stringsAsFactors = T,
                   )

# rownames(TCGAgs) <- TCGAgs[ , 1]
# TCGAgs <- TCGAgs[ , -1]
TCGAgs <- t(TCGAgs)
#######   矩阵的列转化成factor 解决方法1：自编函数########
# colApply <- function(dat, cols = colnames(dat), func = as.factor) { 
#   dat[cols] <- lapply(dat[cols], func) 
#   return(dat)
# } 
# str(colApply(dat = TCGAgs))
# ############## 解决方法2 ： 下面这段代码 
# # 注意：应用lapply的前提条件是df是data.frame， 不能是matrix或array
# df[1:ncol(TCGAgs)] <- lapply(TCGAgs[1:ncol(TCGAgs)], as.factor)

############## 写出TCGAgs.csv文件。以后用来作为gsva的准备数据
write.csv(TCGAgs, "outdata/TCGAgs.csv")
################开始GSVA分析######################
rm(list = ls())
#### 导入准备文件
#test <- unlist(TCGAgs) %>% as.character() %>% unique()
#intersect(test, rownames(mrna)) #50000多的和20000多的基因表达矩阵都只能匹配到438个
# intersect(test, rownames(fpkm))
# rt <- read_csv(file = "outdata/mrna_fpkm_finish.txt")


##上面是我做的一个小实验，想知道PCmrna和总mrna 用来做24免疫细胞的GSVA是否相同
## 结果发现是一样的，50000多的矩阵只多了3个基因而已


library(tidyverse)
fpkm <- read_csv(file = "outdata/mrna_fpkm_id_finish.txt")
fpkm <- column_to_rownames(fpkm, var = colnames(fpkm)[1])
allid <- read.csv("outdata/allid_493.txt") %>% .[,2]
# fpkm <- read.csv(file = "outdata/mrna_fpkm_id_finish.txt",header = T,row.names = 1,check.names = F)
# BiocManager::install("GSVA")
library(GSVA)
library(pheatmap)

TCGAgs <- read.csv("outdata/TCGAgs.csv", stringsAsFactors = T, row.names = 1)
###
TCGAgs<-as.list(TCGAgs)
TCGAgs<-lapply(TCGAgs, function(x) x[!is.na(x)])#去除list中的NA
# 所谓表达矩阵，矩阵矩阵要转换为matrix
fpkm <- as.matrix(fpkm)
# 得到GSVA得分
res_es_gsva <- gsva(fpkm, TCGAgs, 
                    min.sz = 1,
                    max.sz = Inf,
                    mx.diff=T, #mx.diff=FALSE es值是一个双峰的分布T代表近似正态分布
                    verbose=FALSE, 
                    parallel.sz=0,
                    method= "gsva")
# 先做个热图看看
pheatmap(res_es_gsva, cluster_cols = F,cluster_rows = T,show_colnames = F)
res_es_gsva1 <- res_es_gsva # fpkm
allgsva <- res_es_gsva
allgsva <- as.data.frame(t(allgsva))
allgsva <- as.matrix(allgsva)


############ 写出allgsva.csv文件 作为矩阵准备作图##########
write.csv(allgsva,file="outdata/allgsva.csv")
allid <- read.csv("outdata/allid_493.txt") %>% .[,2]
allgsva<-read.csv("outdata/allgsva.csv",row.names = 1) %>% t()
cln_survival <- read.csv("outdata/cln_xena.csv")
pheatmap(mat = as.matrix(allgsva),show_colnames = F,cluster_cols = F)
stopifnot(allid == colnames(allgsva))


######## 分组信息#######
dt <- allgsva
group_list <- c(rep(1,111), rep(0, 382))
group <- factor(group_list, 
                levels = c(0,1), 
                labels = c("Wild","Mut"))
group
design <- model.matrix(~0+group)
rownames(design) <- colnames(dt)
colnames(design) <- levels(group)
design
contrast.matrix <- makeContrasts(Mut - Wild,levels = design)
contrast.matrix
###############  用limma分组做差异分析###########
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

#保存为csv
write.csv(all_diff,file="outdata/gsva_all_diff.csv")

ggscatter(all_diff,x="logFC",y="logP")+ theme_base()

##美化火山图tsv
all_diff <- read.csv("outdata/gsva_all_diff.csv")
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
dev.size("px")


DEG <- read.csv("outdata/gsva_all_diff.csv") 
DEG <- DEG[,c(2,3,7)] #要求DEG有4列，geneid、FC、regulation、pval
DEG <- DEG[is.finite(DEG$logFC),]#去掉无限小数的行
DEG <- DEG[ DEG$adj.P.Val<0.05,] # 筛选logFC大于1.5的
names(DEG) <- c("genes","fold_change","p_value")
DEG$regulation <- "up"
DEG$regulation[DEG$fold_change<0] <- "down"

library(ggplot2)
library(scales)
library(pheatmap)
dd <- rt # dd表示总的表达矩阵
DEG1<- DEG[order(abs(DEG$fold_change),decreasing = T),]#按照FC排序
dd1=dd[rownames(dd) %in% DEG1$genes,]# 在表达矩阵中选取这些基因
dd2=apply(dd1,1,rescale)         ##归一化
dd2=t(dd2)                     ##转置回来

phe3 <- read.csv("outdata/phe_finish.csv",row.names = 1)
rownames(phe3) <-  NULL
phe3 <- phe3[match(x = allid, table = phe3$sample), ]
rownames(phe3) <- phe3$sample
phe3_anno <- phe3 %>% dplyr::select(group, age, gender, smoke_group,ajcc_T, ajcc_N, radiotherapy)
pheatmap(dd2,cluster_rows = T,cluster_cols = F,show_colnames = F,annotation_col = phe3_anno)
### 热图很差，不要了



############## 箱线图绘制###############
allgsva<-read.csv("outdata/allgsva.csv",row.names = 1) %>% t()
cln_survival <- read.csv("outdata/cln_xena.csv")

library(reshape2)
allgsva1 <- t(allgsva)
allgsva1 <- as.data.frame(allgsva1)
allgsva1$group <- cln_survival$group[match(x = rownames(allgsva1), 
                                           table = cln_survival$sample)]
allgsva1 <- reshape2::melt(allgsva1,id.vars=c("group"))
allgsva1$group <- as.factor(allgsva1$group)

#箱线图绘制

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
library(export)
graph2ppt(file= "output/plots/24种细胞GSVAscore")
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
graph2ppt(file= "output/plots/24种细胞GSVAscore",append= T)

############## 小提琴图########
ggviolin(x="variable",y="value",fill = "group",data = allgsva1,size = 0.1,
         palette = c("#00CCFF","#FF3333"),add = "boxplot", 
         add.params = list(fill = "white"))+
  xlab("")+
  ylab("")+
  scale_y_continuous(limits =c(-1.5, 1.5),breaks = seq(-1.5,1.5,0.2))+ 
  stat_compare_means(aes(group=group),label = "p.format",hide.ns = T,
                     bracket.size = 20,label.x = 1.35,label.y = 1800,size = 15)+
  theme_gray()+
  theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5))+
  labs(title = "Mutation Load Plot")+
  theme(plot.title = element_text(size = 55, color = "black", 
                                  hjust = 0.5, 
                                  angle = 0))



# CD8/TREG比值 --------------------------------------------------------------
allgsva2 <- t(allgsva)
allgsva2 <- as.data.frame(allgsva2)
allgsva2$group <- cln_survival$group[match(x = rownames(allgsva2), 
                                           table = cln_survival$sample)]

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
  theme_cleveland()+
  theme(axis.text = element_text(size = 8), 
        axis.text.x = element_text(size = 8, angle = 0, vjust = 0.6))
graph2ppt(file = "output/plots/24种细胞GSVAscore.pptx", append = T)
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
