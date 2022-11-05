rm(list = ls())
load(file = "outdata/LUAD_TCGA_ALLfpkm.rds")
load(file = "outdata/step4_group_data.rds")
p_load(limma,edgeR,ggplot2,ggthemes,ggpubr,ggrepel,export,tidyverse)
# lncRNA注释 ------------------------------------------------------------------
# R包rtracklayer
AnnoData = rtracklayer::import('../maf分组/data/gencode.v22.long_noncoding_RNAs.gtf')
index = which(AnnoData$type == 'gene')
Target = data.frame(Ensembl_ID = AnnoData$gene_id[index], Symbol = AnnoData$gene_name[index], Biotype = AnnoData$gene_type[index])
Target$Ensembl_ID = gsub('\\..*', '', Target$Ensembl_ID) # 正则表达式

# 把data.frame转换成matrix,因为这样就允许行名重复了
rt <- as.matrix(fpkm) 
dim(rt)# 60483个基因
rt[1:3,1:3]
# 去掉不表达的基因
rt <- rt[rowSums(rt) > 0, ]    # 57943个基因  
# 基因ID去掉点号后面的版本标识，准备转换ID。
rownames(rt) <- substring(text = rownames(rt), first = 1, last = 15)
# 确认没重复的基因ENSEMBL ID
stopifnot(
  {
    duplicated(rownames(rt))== F
  }
)
# 查看行和列数
dim(rt) #57943 497(mut279,wild)
common = intersect(Target$Ensembl_ID, rownames(rt))  #####15900中的15373个lncrna
lnc =rt[common,sample_group$sample]
rownames(lnc) <- Target$Symbol[match(x = rownames(lnc),table = Target$Ensembl_ID)]
# 重复的基因名字进行合并且取均值
dim(lnc) # 15373   497
lnc <- avereps(lnc)
dim(lnc) #15364   497
lnc[1:5,1:2]

# 差异分析 --------------------------------------------------------------------

{
  # sample_group$group 
  # group_list <- c(rep(0,279), rep(1, 218))
  # group <- factor(group_list, 
  #                 levels = c(0,1), 
  #                 labels = c("ClusterA","ClusterB"))
  group <- sample_group$group
  group
  design <- model.matrix(~0+group)
  rownames(design) <- colnames(lnc)
  colnames(design) <- levels(group)
  design
  contrast.matrix <- makeContrasts(ClusterB - ClusterA,levels = design)
  contrast.matrix
}

#######进行limma 差异分析######
{
  fit <- lmFit(lnc,design)
  fit1 <- contrasts.fit(fit, contrast.matrix)
  fit1 <- eBayes(fit1)
  qqt(fit1$t, df = fit1$df.prior+fit1$df.residual, pch = 16, cex = 0.2)
  abline(0,1) #QQ图看正态分布?
  all_diff <- topTable(fit1, 
                       adjust.method = 'fdr',
                       coef=1,
                       p.value = 1,
                       lfc = log(1,2),
                       number = Inf,
                       sort.by = 'logFC')
  head(all_diff)
  # 所有差异加一列ID改成gene名字
  all_diff$ID <- rownames(all_diff)
  # 加一列表示logP
  all_diff$logP<- -log10(all_diff$adj.P.Val)
  #保存为csv
  write.csv(all_diff,file="outdata/mitocluster_lnc_all_diff.csv")
}


######简单火山图#########

{
  
  FCvalue <- 0.5
  padj <- 0.05
  #将adj.P.Val<0.05,logFC>0.5的基因设置为显著上调基因
  #将adj.P.Val<0.05,logFC<0.5的基因设置为显著下调基因
  
  all_diff$group = "not-significant"
  all_diff$group[which((all_diff$adj.P.Val< padj) & (all_diff$logFC> FCvalue))] = "up-regulated"
  all_diff$group[which((all_diff$adj.P.Val< padj) & (all_diff$logFC< -FCvalue))] = "down-regulated"
  #查看上调和下调基因的数目
  table(all_diff$group)
  # down-regulated not-significant    up-regulated 
  # 45           15282              37 
  all_diff$group<-as.factor(all_diff$group)
  all_diff$logP<- -log10(all_diff$adj.P.Val)
  #火山图，加上颜色
 ggscatter(all_diff,x="logFC",y="logP",
                 color = "group",
                 palette = c("green","gray","red"),
                 size = 1.5)+theme_base() + 
    geom_hline(yintercept = 1.30,linetype = "dashed") + 
    geom_vline(xintercept = c(-FCvalue,FCvalue),linetype = "dashed")
}
#加上排名前十的基因
{
  all_diff$label = ""
  #对差异基因的p值由小到大排序，学会一个order算法！
  all_diff<-all_diff[order(all_diff$adj.P.Val),]
  
  #高表达的基因中，选择adj.P.val最小的10个
  up.genes<-head(all_diff$ID[which(all_diff$group == "up-regulated")],10)
  #低表达的基因中，选择adj.P.val最小的10个
  down.genes<-head(all_diff$ID[which(all_diff$group == "down-regulated")],10)
  #以上两步合并加入label中
  all_diff.top10.genes <- c(as.character(up.genes),as.character(down.genes))
  all_diff$label[match(all_diff.top10.genes,all_diff$ID)] <- all_diff.top10.genes
}

#### 火山图最终成图########
p <- ggplot(data = all_diff, 
            aes(x = logFC, y = logP)) +
  geom_point(alpha=0.7, size=2, aes(color=group) ) +
  scale_color_manual(values = c("#00468BCC","gray","#ED0000CC"))+
  geom_vline(xintercept=c(-FCvalue,FCvalue),lty="dashed",col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty="dashed",col="black",lwd=0.8) +
  geom_text_repel(aes(label = label), box.padding = unit(1, "lines"), 
                  point.padding = unit(1, "lines"), show.legend = F, 
                  segment.color = 'black', size = 3,max.overlaps = 40)+
  theme_minimal()+
  theme(axis.text.x = element_text(size = 15, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.position = "right"
  )
p
library(export)
graph2ppt(file="output/plots/DEanalysis.pptx", 
          width = 8.5, aspectr = 1.5, append = T)



########热图绘制


# 热图绘制 --------------------------------------------------------------------

#data presentation 数据准备好了
{
  library(tidyverse)
  DEG <- all_diff |> dplyr::select(ID,logFC,adj.P.Val)
  #要求DEG有4列，geneid、FC、regulation、pval
  DEG <- DEG[is.finite(DEG$logFC),]#去掉无限小数的行
  DEG <- DEG[abs(DEG$logFC)> FCvalue & DEG$adj.P.Val< padj,] # 筛选logFC大于1的
  names(DEG) <- c("genes","fold_change","p_value")
  DEG$regulation <- "up"
  DEG$regulation[DEG$fold_change<0] <- "down"
  table(DEG$regulation)
}

#pheatmap 做热图
{
  library(scales)
  library(pheatmap)
  # dd表示总的表达矩阵
  dd <- lnc
  DEG
  DEG1<- DEG[order(abs(DEG$fold_change),decreasing = T),]#按照FC排序
  DEG1 <- DEG1[1:40,]#只取前40个基因，要求DEG有4列，geneid、FC、regulation、pval
  dd1=dd[rownames(dd) %in% DEG1$genes, ]# 在表达矩阵中选取这些基因
  dd2=apply(dd1,1,rescale)         ##归一化
  dd2=t(dd2)                     ##转置回来
  pheatmap(dd2,cluster_rows = T,
           # scale = "column",
           cluster_cols = F,show_colnames = F)  # 做热图
}


load(file = "outdata/LUAD_cln_clean.rds")
intersect(colnames(lnc), sample_group$sample)
phe3 <- phe_f
phe3 <- phe3[match(x = colnames(dd), table = phe3$sample), ]
dim(phe3)
phe3$group <- sample_group$group[match(phe3$sample, sample_group$sample,nomatch = 0)]
intersect(phe3$sample, sample_group$sample)
phe3_anno <- phe3 %>% dplyr::select(group, age, gender, smoke_group,ajcc_T, ajcc_N, radiotherapy)
rownames(phe3_anno) <- colnames(dd2)

pheatmap(dd2,cluster_rows = T,cluster_cols = F,show_colnames = F,
         annotation_col = phe3_anno)

graph2ppt(file="output/plots/DEanalysis.pptx", width = 8.5, aspectr = 1,append = T)

