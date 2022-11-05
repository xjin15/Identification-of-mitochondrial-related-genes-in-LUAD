rm(list = ls())
# mirna -------------------------------------------------------------------
load(file = "outdata/step4_group_data.rds")


############### 导入mirna的数据并做整理构建表达矩阵 ##############
p_load(limma,edgeR, ggplot2,ggthemes, ggpubr, tidyverse, ggrepel)
############# mirna 数据加载##########
# mirna <- read.csv(file = "../maf分组/data/miRNA_HiSeq_gene", 
#                   header = T, sep = "\t", stringsAsFactors = F, 
#                   row.names = 1,check.names = F)
mirna <- data.table::fread(file = "../maf分组/data/miRNA_HiSeq_gene") |> 
  column_to_rownames("sample")
  

# TCGA的名称.改成-
# colnames(mirna_stemloop) <- gsub('\\.', '-', colnames(mirna_stemloop))
mirna <- mirna
# 确认没有重复的列名（样本名称）
stopifnot({length(colnames(mirna))==length(unique(colnames(mirna)))})
# 只保留fpkm样品中结尾01（初发癌）的样品 564减少到510
mirna <- mirna[,substring(text = colnames(mirna), 
                                        first = 14) == "01"]
colnames(mirna) <- paste0(colnames(mirna), "A")
dim(mirna) # 去之前有mirna 2228个,448个样本 
# 去掉在所有样品都不表达的基因（mirna=0）
# 先转换成array 矩阵而不 是dataframe
rt <- as.matrix(mirna) 
rt <- replace_na(rt,0) # 用数字0替换NA值
rt <- rt[rowSums(rt) > 0, ]  
dim(rt) # 去以后剩mirna 2221个
### 挑选出需要的sample
intersample <- intersect(allid,colnames(rt))
rt1 <-rt[ , intersample] # 429个样本
length(intersect(intersample,sample_group$sample[1:279])) #前253个是clusterA
length(intersect(intersample,allid[280:497]))# 后176是clusterB

############## id转换############
name <- rownames(rt1)
p_load(miRBaseVersions.db) # 一个miRNA数据库R包。可选择对应的版本和ID号以及名字
items <- select(miRBaseVersions.db,
                keys = name,
                keytype = "MIMAT",
                columns = c("ACCESSION","NAME","VERSION"))
# 选择版本为20.0 以及name
id_name <- items[items$VERSION == 20.0, c("ACCESSION","NAME")]
# match一下我们要的
id_name <- id_name[match(x=rownames(rt1), table = id_name$ACCESSION), ]
rt2 <- cbind(rt1,id_name) %>% dplyr::select(ACCESSION, NAME, everything())
stopifnot({rownames(rt2)==rt2[,1]})
rownames(rt1) <- rt2[,2]
write.csv(rt1, file = "outdata/miRNA_rpm.txt", row.names = T)
############   用limma做差异分析      ###########
dt <- read.csv(file = "outdata/miRNA_rpm.txt", header = T,
               check.names = F, 
               row.names = 1,
               stringsAsFactors = F)
rt <- as.matrix(dt)
######## 创建分组信息
{
  group_list <- c(rep(0,253), rep(1, 176))
  group <- factor(group_list, 
                  levels = c(0,1), 
                  labels = c("A","B"))
  group |> table()
  design <- model.matrix(~0+group)
  rownames(design) <- colnames(dt)
  colnames(design) <- levels(group)
  design
  contrast.matrix <- makeContrasts(B - A,levels = design)
  contrast.matrix
}

#######进行limma 差异分析######

{
  fit <- lmFit(rt,design)
  fit1 <- contrasts.fit(fit, contrast.matrix)
  fit1 <- eBayes(fit1)
  qqt(fit1$t, df = fit1$df.prior+fit1$df.residual, pch = 16, cex = 0.2)
  # Zero sample variances detected, have been offset away from zero 
  #这只是告诉你至少有一个gene它的表达值在所有的样品中没有任何变化，这个gene会在计算中被忽略.
  abline(0,1) #QQ图看正态分布?
  all_diff <- topTable(fit1, adjust.method = 'fdr',coef=1,p.value = 1,lfc = log(1,2),number = Inf,sort.by = 'logFC')
  head(all_diff)
  # 所有差异加一列ID改成gene名字
  all_diff$ID <- rownames(all_diff)
  # 加一列表示logP
  all_diff$logP<- -log10(all_diff$adj.P.Val)
  #保存为csv
  write.csv(all_diff,file="outdata/mitocluster_mirna_all_diff.csv")
  
}

##############火山图绘制####################
{ 
  ## 得出所有的差异基因
  all_diff$group = "not-significant"
  FCvalue <- 0.5
  padj <- 0.05
  #将adj.P.Val<0.05,logFC>0.5的基因设置为显著上调基因
  #将adj.P.Val<0.05,logFC<0.5的基因设置为显著下调基因
  all_diff$group[which((all_diff$adj.P.Val< padj) & (all_diff$logFC> FCvalue))] = "up-regulated"
  all_diff$group[which((all_diff$adj.P.Val< padj) & (all_diff$logFC< -FCvalue))] = "down-regulated"
  #查看上调和下调基因的数目
  table(all_diff$group)
  # down-regulated not-significant    up-regulated 
  #  34            2122              65 
  all_diff$group<-as.factor(all_diff$group)
  all_diff$logP<- -log10(all_diff$adj.P.Val)
  #火山图，加上颜色
   ggscatter(all_diff,x="logFC",y="logP",
                 color = "group",
                 palette = c("green","gray","red"),
                 size = 1.5)+theme_base() + 
    geom_hline(yintercept = 1.30,linetype = "dashed") + 
    geom_vline(xintercept = c(-FCvalue,FCvalue),linetype = "dashed")
  
  # dev.size("px")
  # ggsave(p,filename = "output/KEAP1_diffgene_number_volcanoplot.pdf")      ##,width = ,height = )
}
#加上排名前十的基因
{
  all_diff$label = ""
  #对差异基因的p值由小到大排序，学会一个order算法！
  all_diff<-all_diff[order(all_diff$adj.P.Val),]
  #
  #高表达的基因中，选择adj.P.val最小的10个
  up.mirs<-head(all_diff$ID[which(all_diff$group == "up-regulated")],10)
  #低表达的基因中，选择adj.P.val最小的10个
  down.mirs<-head(all_diff$ID[which(all_diff$group == "down-regulated")],10)
  #以上两步合并加入label中
  all_diff.top10.mirs <- c(as.character(up.mirs),as.character(down.mirs))
  all_diff$label[match(all_diff.top10.mirs,all_diff$ID)] <- all_diff.top10.mirs
}

# 画火山图,repel=T表示字体不会重叠#   一运行就直接Rsession破灭
{
    p <- ggplot(data = all_diff, 
              aes(x = logFC, y = logP)) +
    geom_point(alpha=0.7, size=2, aes(color=group) ) +
    scale_color_manual(values = c("#00468BCC","gray","#ED0000CC"))+
    geom_vline(xintercept=c(-0.5,0.5),lty="dashed",col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty="dashed",col="black",lwd=0.8) +
    geom_text_repel(aes(label = label), box.padding = unit(1, "lines"), 
                    point.padding = unit(1, "lines"), show.legend = F, 
                    segment.color = 'black', size = 3,max.overlaps = 30)+
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
}



###### 绘制热图#####
#data presentation 数据准备好了
{
  DEG <- read.csv("outdata/mitocluster_mirna_all_diff.csv")
  DEG <- DEG[,c("X","logFC","adj.P.Val")] #要求DEG有4列，geneid、FC、regulation、pval
  DEG <- DEG[is.finite(DEG$logFC),]#去掉无限小数的行
  DEG <- DEG[abs(DEG$logFC)>0.5 & DEG$adj.P.Val<0.05,] # 筛选logFC大于1.5的
  names(DEG) <- c("genes","fold_change","p_value")
  DEG$regulation <- "up"
  DEG$regulation[DEG$fold_change<0] <- "down"
}

#pheatmap 做热图
{
  library(scales)
  library(pheatmap)
  dd <- read.csv(file = "outdata/miRNA_rpm.txt", header = T,
                 check.names = F, ## 不检查列名！
                 row.names = 1,
                 stringsAsFactors = F) # dd表示总的表达矩阵
  DEG1<- DEG[order(abs(DEG$fold_change),decreasing = T),]#按照FC排序
  DEG1 <- DEG1[1:30,]#只取前40个基因，要求DEG有4列，geneid、FC、regulation、pval
  dd1=dd[rownames(dd) %in% DEG1$genes, ]# 在表达矩阵中选取这些基因
  dd2=apply(dd1,1,rescale)         ##归一化
  dd2=t(dd2)                     ##转置回来
  pheatmap(dd2,cluster_rows = T,cluster_cols = F,show_colnames = F,scale = "none")  # 做热图
  # pdf(file="output/Top_40genes.pdf",width = 3,height = 5.5)
  # pheatmap(dd2,cutree_rows = 2,cutree_cols = 2,fontsize_row = 8 )
  # dev.off()
}

load(file = "outdata/LUAD_cln_clean.rds")
intersect(colnames(mirna), sample_group$sample)
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
