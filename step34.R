#################使用ceRNA做分析##############
 # BiocManager::install("GDCRNATools")
rm(list = ls())
library(GDCRNATools)
library(tidyverse)
# ceoutput <-  gdcCEAnalysis(lnc = deLNC,
#               pc = dePC,
#               deMIR = NULL,
#               lnc.targets = 'starBase',
#               pc.targets = 'starBase',
#               rna.expr = rnaExpr,
#               mir.expr = mirExpr)
# 这是示例

### 不是只要差异mrna和lncrna 的表达矩阵，表达矩阵的参数不用调节#############
# rnaExpr_up <- rnaExpr[c(deLNC_up1, dePC_up1), ]
# rnaExpr_down <- rnaExpr[c(deLNC_down1, dePC_down1), ]
# rnaExpr <- read_csv(file = "../maf分组/outdata/mrna_fpkm_finish.txt") %>% column_to_rownames(.,colnames(.)[1])

load("outdata/LUAD_TCGA_ALLfpkm.rds")
load("outdata/step4_group_data.rds")
fpkm[1:3,1:3]
rnaExpr <- fpkm[,sample_group$sample]
rownames(rnaExpr) <- base::gsub("\\..*", "", rownames(rnaExpr))
rnaExpr[1:3,1:3]
rnaExpr <- as.data.frame(rnaExpr)
###################最后一步#################
mirExpr <- read.csv(file = "outdata/miRNA_rpm.txt",
                    header = T,
                    check.names = F, 
                    row.names = 1,
                    stringsAsFactors = F) %>% as.data.frame()
mirExpr[1:3,1:3]
################取出dePC和deLNC#####################
DEG <- read.csv("outdata/mitocluster_lnc_all_diff.csv")
DEG <- DEG[,c(1,2,5)] #要求DEG有4列，geneid、FC、regulation、pval
DEG <- DEG[is.finite(DEG$logFC),]#去掉无限小数的行
DEG <- DEG[abs(DEG$logFC)> 1 & DEG$P.Value<0.05,] # 筛选logFC大于1.5的
names(DEG) <- c("genes","fold_change","p_value")
DEG$regulation <- "up"
DEG$regulation[DEG$fold_change<0] <- "down"
table(DEG$regulation)
deLNC_up <- DEG$genes[DEG$regulation == "up"]
deLNC_down <- DEG$genes[DEG$regulation == "down"]
deLNC_all <- DEG$genes

DEG <- read.csv("outdata/mitocluster_mrna_all_diff.csv")
DEG <- DEG[,c(1,2,5)] #要求DEG有4列，geneid、FC、regulation、pval
DEG <- DEG[is.finite(DEG$logFC),]#去掉无限小数的行
DEG <- DEG[abs(DEG$logFC) > 0.5 & DEG$P.Value<0.05,] # 筛选logFC大于1.5的
names(DEG) <- c("genes","fold_change","p_value")
DEG$regulation <- "up"
DEG$regulation[DEG$fold_change<0] <- "down"
dePC_up <- DEG$genes[DEG$regulation == "up"]
dePC_down <- DEG$genes[DEG$regulation == "down"]
dePC_all <- DEG$genes

DEG <- read.csv("outdata/mitocluster_mirna_all_diff.csv")
DEG <- DEG[,c(1,2,5)] #要求DEG有4列，geneid、FC、regulation、pval
DEG <- DEG[is.finite(DEG$logFC),]#去掉无限小数的行
DEG <- DEG[abs(DEG$logFC)> 0.25 & DEG$P.Value<0.05,] # 筛选logFC大于1.5的
names(DEG) <- c("genes","fold_change","p_value")
DEG$regulation <- "up"
DEG$regulation[DEG$fold_change<0] <- "down"
deMIR_all <- DEG$genes

#################### 一个教训又要来做ID转换了 ##############
## 要从symbol 转成 gene id
gtf <- read.csv(file = "outdata/gtfmRNA22.txt", header = T,sep = "\t")
deLNC_all1 <- gtf$gene_id[match(x = deLNC_all, table = gtf$gene_name)]
dePC_all1 <- gtf$gene_id[match(x = dePC_all, table = gtf$gene_name)]
# deLNC_up1 <- gtf$gene_id[match(x = deLNC_up, table = gtf$gene_name)]
# deLNC_down1 <- gtf$gene_id[match(x = deLNC_down, table = gtf$gene_name)]
# dePC_up1 <- gtf$gene_id[match(x = dePC_up, table = gtf$gene_name)]
# dePC_down1 <- gtf$gene_id[match(x = dePC_down, table = gtf$gene_name)]


# ceRNA最后一步 ---------------------------------------------------------------
ceOutput_star <- gdcCEAnalysis(lnc = deLNC_all1,
                               pc = dePC_all1,
                               deMIR = deMIR_all,
                               lnc.targets = 'starBase',
                               pc.targets = 'starBase',
                               rna.expr = rnaExpr,
                               mir.expr = mirExpr)
## miRcode 做不出来,
ceOutput_miRcode <- gdcCEAnalysis(lnc = deLNC_all1,
                                  pc = dePC_all1,
                                  deMIR = deMIR_all,
                                  lnc.targets = 'miRcode',
                                  pc.targets = 'miRcode',
                                  rna.expr = rnaExpr,
                                  mir.expr = mirExpr)
## spongescan 也做不出来1
ceOutput_spongescan <- gdcCEAnalysis(lnc = deLNC_all1,
                                     pc = dePC_all1,
                                     deMIR = deMIR_all,
                                     lnc.targets = 'spongeScan',
                                     pc.targets = 'spongeScan',
                                     rna.expr = rnaExpr,
                                     mir.expr = mirExpr)

ceOutput2 <- ceOutput_star[ceOutput_star$hyperPValue<0.01 &
                             ceOutput_star$corPValue<0.01 & ceOutput_star$regSim != 0,]
ceOutput_final <- ceOutput_star

edges <- gdcExportNetwork(ceNetwork = ceOutput_final, net = 'edges')
nodes <- gdcExportNetwork(ceNetwork = ceOutput_final, net = 'nodes')

write.table(edges, file='outdata/ceRNA_cyto_edges.txt', sep='\t', quote=F,row.names = F)
write.table(nodes, file='outdata/ceRNA_cyto_nodes.txt', sep='\t', quote=F,row.names = F)
write.csv(ceOutput_star, file = "outdata/ceRNA_output.csv")
cln_data <- read_csv("outdata/cln_xena.csv") %>% arrange(group,sample) %>% .[,-3]
shinyCorPlot(gene1    = deLNC_all1, 
             gene2    = dePC_all1, 
             rna.expr = rnaExpr, 
             metadata = cln_data)

genelist <- gtf$gene_name[match(x = ceOutput_star$Genes, gtf$gene_id)]

genelist2 <- all_diff.top10.genes

setdiff(genelist,genelist2)
intersect(genelist,genelist2)



{
  # ceOutput_up <-  gdcCEAnalysis(lnc = deLNC_up1,
  #                            pc = dePC_up1,
  #                            deMIR = NULL,
  #                            lnc.targets = 'starBase',
  #                            pc.targets = 'starBase',
  #                            rna.expr = rnaExpr_up,
  #                            mir.expr = mirExpr)
  
  
  # ceOutput_up_all <-  gdcCEAnalysis(lnc = deLNC_up1,
  #                               pc = dePC_up1,
  #                               deMIR = NULL,
  #                               lnc.targets = 'starBase',
  #                               pc.targets = 'starBase',
  #                               rna.expr = rnaExpr,
  #                               mir.expr = mirExpr)
  
  
  # ceOutput_down <- gdcCEAnalysis(lnc = deLNC_down1,
  #                                pc = dePC_down1,
  #                                deMIR = NULL,
  #                                lnc.targets = 'starBase',
  #                                pc.targets = 'starBase',
  #                                rna.expr = rnaExpr_down,
  #                                mir.expr = mirExpr)
  # ceOutput_down_all <- gdcCEAnalysis(lnc = deLNC_down1,
  #                                pc = dePC_down1,
  #                                deMIR = NULL,
  #                                lnc.targets = 'starBase',
  #                                pc.targets = 'starBase',
  #                                rna.expr = rnaExpr,
  #                                mir.expr = mirExpr)
  
  # ceOutput_all_all <- gdcCEAnalysis(lnc = deLNC_all1,
  #                                    pc = dePC_all1,
  #                                    deMIR = NULL,
  #                                    lnc.targets = 'starBase',
  #                                    pc.targets = 'starBase',
  #                                    rna.expr = rnaExpr,#########这个参数不用强行改，就原始的表达矩阵就好了
  #                                    mir.expr = mirExpr)
}