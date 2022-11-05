rm(list = ls())
load(file = "outdata/step4.2_mrna.rds")
load(file = "outdata/step4_group_data.rds")
library(pacman)
p_load(limma,ggpubr,tidyverse,ggthemes)
########### 设计limma比较分组矩阵以及PCA################
{
  p_load(FactoMineR,factoextra) # PCA用到的包
  allid497
  
  sample_group$group |> table()
  
  group_list <- c(rep(1,279), rep(0, 218))
  group <- factor(group_list, 
                  levels = c(0,1), 
                  labels = c("A","B"))
  group
  design <- model.matrix(~0+group)
  rownames(design) <- colnames(mrna)
  colnames(design) <- levels(group)
  design
  contrast.matrix <- makeContrasts(A - B,levels = design)
  contrast.matrix
  # voom(mrna, design, plot = T) 过滤掉低表达基因以后就会有一条很平滑的曲线
  ##### PCA图
  {
    dat <-  as.data.frame(t(mrna)) # 画PCA图时要求是行名是样本名，列名时探针名，因此此时需要转换。格式要求data.frame
    dat.pca <- PCA(dat, graph = F)
    # fviz_pca_ind按样本  fviz_pca_var按基因
    fviz_pca_ind(dat.pca,
                 geom.ind = "point", # c("point", "text)2选1
                 col.ind = group, # color by groups
                 # palette = c("#00AFBB", "#E7B800"),# 自定义颜色
                 addEllipses = T, # 加圆圈
                 legend.title = "Groups"# 图例名称
    )
    # export::graph2ppt(file="output/plots/PCA.pptx",width = 5.6,heig = 5, append = T)
    
    # plotMDS(mrna, col = as.numeric(group))
  }
}

#######进行limma 差异分析######
{
  fit <- lmFit(mrna,design)
  fit1 <- contrasts.fit(fit, contrast.matrix)
  fit1 <- eBayes(fit1)
  qqt(fit1$t, df = fit1$df.prior+fit1$df.residual, pch = 16, cex = 0.2)
  abline(0,1) #QQ图看正态分布?
  all_diff <- topTable(fit1, 
                       adjust.method = 'fdr',
                       coef=1,
                       p.value = 1,
                       lfc <- log(1,2),
                       number = Inf,
                       sort.by = 'logFC')
  head(all_diff)
  # 所有差异加一列ID改成gene名字
  all_diff$ID <- rownames(all_diff)
  # 加一列表示logP
  all_diff$logP<- -log10(all_diff$adj.P.Val)
  #保存为csv
  write.csv(all_diff,file="outdata/mitocluster_mrna_all_diff.csv")
}


######简单火山图#########

{
  
  FCvalue <- 1
  padj <- 0.05
  #将adj.P.Val<0.05,logFC>0.5的基因设置为显著上调基因
  #将adj.P.Val<0.05,logFC<0.5的基因设置为显著下调基因
  
  all_diff$group <-  "not-significant"
  all_diff$group[which((all_diff$adj.P.Val< padj) & (all_diff$logFC> FCvalue))] = "up-regulated"
  all_diff$group[which((all_diff$adj.P.Val< padj) & (all_diff$logFC< -FCvalue))] = "down-regulated"
  #查看上调和下调基因的数目
  table(all_diff$group)
  all_diff$group <- as.factor(all_diff$group)
  all_diff$logP <- -log10(all_diff$adj.P.Val)
  #火山图，加上颜色
  ggscatter(all_diff,x="logFC",y="logP",
            color = "group",
            palette = c("green","gray","red"),
            size = 1.5)+theme_base()
  #再加上辅助线
  p <- ggscatter(all_diff,x="logFC",y="logP",
                 color = "group",
                 palette = c("green","gray","red"),
                 size = 1.5)+theme_base() + 
    geom_hline(yintercept = 1.30,linetype = "dashed") + 
    geom_vline(xintercept = c(-FCvalue,FCvalue),linetype = "dashed")
  p
  # dev.size("px")
  # ggsave(p,filename = "output/KEAP1_diffgene_number_volcanoplot.pdf")      ##,width = ,height = )
}
#加上排名前十的基因
{
  all_diff$label <-  ""
  #对差异基因的p值由小到大排序，学会一个order算法！
  all_diff <- all_diff[order(all_diff$adj.P.Val),]
  #
  all_diff$X <- rownames(all_diff)
  all_diff$X <- NULL
  #高表达的基因中，选择adj.P.val最小的10个
  up.genes <- head(all_diff$ID[which(all_diff$group == "up-regulated")],10)
  #低表达的基因中，选择adj.P.val最小的10个
  down.genes <- head(all_diff$ID[which(all_diff$group == "down-regulated")],10)
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
  ggrepel::geom_text_repel(aes(label = label), box.padding = unit(1, "lines"), 
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
graph2ppt(file="output/plots/DEanalysis.pptx", width = 7, aspectr = 1.5, append = T)



########热图绘制


# 热图绘制 --------------------------------------------------------------------

#data presentation 数据准备好了
{
  DEG <- all_diff |> dplyr::select(ID,logFC,adj.P.Val)
  #要求DEG有4列，geneid、FC、regulation、pval
  DEG <- DEG[is.finite(DEG$logFC),]#去掉无限小数的行
  DEG <- DEG[abs(DEG$logFC)> FCvalue & DEG$adj.P.Val< padj,] # 筛选logFC大于1的
  names(DEG) <- c("genes","fold_change","p_value")
  DEG$regulation <- "up"
  DEG$regulation[DEG$fold_change<0] <- "down"
  table(DEG$regulation)
  save(DEG,all_diff,file = "outdata/step4.3_enrich_Analysis.rds")
}

#pheatmap 做热图
{
  library(scales)
  library(pheatmap)
  # dd表示总的表达矩阵
  dd <- mrna
  DEG1<- DEG[order(abs(DEG$fold_change),decreasing = T),]#按照FC排序
  DEG1 <- DEG1[1:50,]#只取前40个基因，要求DEG有4列，geneid、FC、regulation、pval
  dd1=dd[rownames(dd) %in% DEG1$genes, ]# 在表达矩阵中选取这些基因
  dd2=apply(dd1,1,rescale)         ##归一化
  dd2=t(dd2)                     ##转置回来
  pheatmap(dd2,cluster_rows = T,
           # scale = "column",
           cluster_cols = F,show_colnames = F)  # 做热图
}



#############给热图加上注释信息，也就是survival信息###############

{  
  load(file = "outdata/LUAD_cln_clean.rds")
  phe3 <- phe_f
  rownames(phe3) <-  NULL
  phe3 <- phe3[match(x = allid497, table = phe3$sample), ]
  rownames(phe3) <- phe3$sample
  phe3$group <- sample_group$group[match(x = phe3$sample, table = sample_group$sample)]
  phe3_anno <- phe3 %>% dplyr::select(group, age, gender, smoke_group,ajcc_T, ajcc_N, radiotherapy)
  pheatmap(dd2,cluster_rows = T,cluster_cols = F,show_colnames = F,annotation_col = phe3_anno)
  
  # annotation_col = data.frame( CellType = factor(rep(c("CT1", "CT2"), 5)), Time = 1:5 )
  # rownames(annotation_col) = paste("Test", 1:10, sep = "")
  # annotation_row = data.frame( GeneClass = factor(rep(c("Path1", "Path2", "Path3"), c(10, 4, 6))))
  # rownames(annotation_row) = paste("Gene", 1:20, sep = "")
  graph2ppt(file="output/plots/DEanalysis.pptx", width = 7, aspectr = 1,append = T)
  # pdf(file="output/Top_40genes.pdf",width = 3,height = 5.5)
  # pheatmap(dd2,cutree_rows = 2,cutree_cols = 2,fontsize_row = 8 )
  # dev.off()
}

