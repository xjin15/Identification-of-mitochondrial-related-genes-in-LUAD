#################用R包预测miRNA的靶基因，然后对靶基因进行富集分析#################
# BiocManager::install("multiMiR",ask = F,update = F)
library(multiMiR)
# db.ver = multimir_dbInfoVersions()
# db.ver[,1:3]

# The default is to search validated interactions in human
example1 <- get_multimir(mirna = 'hsa-miR-18a-3p', summary = TRUE)
names(example1)
# Check which types of associations were returned
table(example1@data$type)

# Detailed information of the validated miRNA-target interaction
head(example1@data)
dim(example1@data)

# Which interactions are supported by Luciferase assay?
example1@data[grep("Luciferase", example1@data[, "experiment"]), ]
example1@summary[example1@summary[,"target_symbol"] == "KRAS",]
summary(example1)
#查询多个mirna的靶基因
multimir_results <- get_multimir(org = "hsa",
                                 mirna   = all_diff.top10.mirs,
                                 table   = 'validated',
                                 summary = T,
)
table(multimir_results@data$type)
dim(multimir_results@data)
head(multimir_results@data)
multimir_results@data[grep("Luciferase", multimir_results@data[, "experiment"]), ]

example1 <- get_multimir(org = "hsa",
                         mirna   = up.mirs,
                         target = down.genes,
                         table   = 'all',
                         summary = T ,
                         predicted.cutoff = 10,
                         predicted.cutoff.type = "p",
                         use.tibble = T)
table(example1@data$type)
result <- select(example1, keytype = "type", keys = "validated", columns = columns(example1))
unique_pairs <- 
  result[!duplicated(result[, c("mature_mirna_id", "target_entrez")]), ]

result


example2 <- get_multimir(org = "hsa",
                         mirna   = down.mirs,
                         target = up.genes,
                         table   = 'all',
                         summary = T ,
                         predicted.cutoff = 10,
                         predicted.cutoff.type = "p",
                         use.tibble = T)
table(example2@data$type)
result2 <- select(example2, keytype = "type", keys = "validated", columns = columns(example2))
unique_pairs2 <- 
  result2[!duplicated(result2[, c("mature_mirna_id", "target_entrez")]), ]

result2

# ############富集分析

#enrichment analysis 富集分析：可以用R语言的clusterprofiler或者用david网站

{
  #   if (!requireNamespace("BiocManager", quietly = TRUE))
  #     install.packages("BiocManager")
  #   
  #   BiocManager::install("clusterProfiler")
  #   if (!requireNamespace("BiocManager", quietly = TRUE))
  #     install.packages("BiocManager")
  #   
  #   BiocManager::install("org.Hs.eg.db")
  #   BiocManager::install("DO.db")
  ##安装clusterProfiler,org.Hs.eg.db,
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(cowplot)
  keytypes(org.Hs.eg.db)
  
  
  gene <- bitr(DEG$genes, fromType ="SYMBOL",##把差异表达基因名称转换成geneID和SYMBOL
               toType = "ENTREZID", 
               OrgDb = org.Hs.eg.db)
  geneList <- bitr(rownames(dd), fromType = "SYMBOL",##dd里面是所有的基因
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
  ego <- enrichGO(gene           = gene$ENTREZID,
                  universe       = names(geneList$ENTREZID),
                  OrgDb          = org.Hs.eg.db,
                  ont            = "BP", 
                  pAdjustMethod  = "BH",
                  pvalueCutoff   = 0.01,
                  qvalueCutoff   = 0.05,
                  readable       = T )              ##GO分析
  
  kk <- enrichKEGG(gene          = gene$ENTREZID,
                   organism      = 'hsa',
                   pvalueCutoff  = 0.05,
                   pAdjustMethod = "BH", 
                   qvalueCutoff  = 0.05)            ##KEGG分析  
  
  p1 <- dotplot(object = ego, showCategory = 10, orderBy = "x") + ggtitle("dotplot for GOBP")
  p2 <- dotplot(object = kk,  showCategory = 10, orderBy = "x") + ggtitle("dotplot for KEGG")
  p1
  graph2ppt(file="output/plots/mrna_go_plot.pptx")
  p2
  graph2ppt(file="output/plots/mrna_kegg_plot.pptx")
  
  # ggsave(pp,filename = "output/plots/DEG_enrichment.pdf",width = 12,height = 3.8)
  write.table(DEG,file = "output/DEG_keap1.xls",sep = "\t",quote = F,row.names = F)##输出差异表达基因
}