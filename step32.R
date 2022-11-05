rm(list = ls())
p_load(tidyverse,ggplot2, limma, export)

library(pheatmap)
# 准备数据
load("outdata/step4.2_mrna.rds")
load("outdata/step4_group_data.rds")
mrna <- t(mrna)
######0.免疫基因数据读取：包括固有免疫功能、抗原提呈能力##########
{
{
  immugeneMHC <- data.table::fread(
    "Antigen presentation
HLA-A
HLA-B
HLA-C
HLA-DPA1
HLA-DPB1
HLA-DQA1
HLA-DQA2
HLA-DQB1
HLA-DQB2
HLA-DRA
HLA-DRB1
HLA-DRB3
HLA-DRB4
HLA-DRB5
MICA
MICB
") %>% 
    pull()
  
} # 读取抗原呈递的MHC
{
  immugeneINHI <- data.table::fread(
    "inhibitory
  ADORA2A
ARG1
BTLA
CD274
CD276
CTLA4
EDNRB
HAVCR2
IDO1
IL10
IL13
IL4
KIR2DL1
KIR2DL2
KIR2DL3
LAG3
PDCD1
SLAMF7
TGFB1
TIGIT
VEGFA
VEGFB
C10orf54
VTCN1
"
  ) %>% pull()
} # 读取免疫检查点 抑制免疫作用的
{
immugeneSTIM <- data.table::fread(
  "Stimulatory
  BTN3A1
BTN3A2
CCL5
CD27
CD28
CD40
CD40LG
CD70
CD80
CX3CL1
CXCL10
CXCL9
ENTPD1
GZMA
HMGB1
ICAM1
ICOS
ICOSLG
IFNA1
IFNA2
IFNG
IL12A
IL1A
IL1B
IL2
IL2RA
ITGB2
PRF1
SELP
TLR4
TNF
TNFRSF14
TNFRSF18
TNFRSF4
TNFRSF9
TNFSF4
TNFSF9
"
) %>% pull()} # 读取免疫检查点 促进免疫作用的


immugene <- c('IRF3','MYD88','TICAM1','TLR3','TLR5','TLR7',
              'TLR8','DDX58','IFIH1','MAVS','CLEC7A','CLEC4E',
              'CD209','CLEC10A','NLRP3','AIM2','PYCARD', #innate免疫共17个分子initiation of innate immunity
              'HLA-A','HLA-B','HLA-C','TAP1','TAP2','B2M','HLA-DPA1',
              'HLA-DPB1','HLA-DQA1','HLA-DQB1','HLA-DQB2',
              'HLA-DRB4','HLA-DRB6','HLA-E','HLA-F','HLA-J', #抗原提呈系列MHC-I/II antigen-presenting process的HLA共16个分子
              "ADORA2A", "ARG1", "BTLA", "CD274", "CD276", "CTLA4", "EDNRB", 
              "HAVCR2", "IDO1", "IL10", "IL13", "IL4", "KIR2DL1", "KIR2DL2", 
              "KIR2DL3", "LAG3", "PDCD1", "SLAMF7", "TGFB1", "TIGIT", "VEGFA", 
              "VEGFB", "C10orf54", "VTCN1", # 免疫检查点抑制基因 24个
              "BTN3A1", "BTN3A2", "CCL5", "CD27", 
              "CD28", "CD40", "CD40LG", "CD70", "CD80", "CX3CL1", "CXCL10", 
              "CXCL9", "ENTPD1", "GZMA", "HMGB1", "ICAM1", "ICOS", "ICOSLG", 
              "IFNA1", "IFNA2", "IFNG", "IL12A", "IL1A", "IL1B", "IL2", "IL2RA", 
              "ITGB2", "PRF1", "SELP", "TLR4", "TNF", "TNFRSF14", "TNFRSF18", 
              "TNFRSF4", "TNFRSF9", "TNFSF4", "TNFSF9" ) #免疫检查点激活基因37个

immportgenes <- data.table::fread("../data/IMMUNE_data/GeneList_immport.txt") 
}
save(immugene,immugeneINHI,immportgenes,immugeneMHC,immugeneSTIM,file = "../data/IMMUNE_data/ALL_immu_genesets.Rdata")



# fpkm数据整理 ----------------------------------------------------------------
load('outdata/LUAD_TCGA_ALLfpkm.rds')
# 基因名改成symbol，相同基因名的

fpkm <- as.matrix(fpkm)
rownames(fpkm)
genesid <- rownames(fpkm)
genesid <- substring(genesid,1,15)
genesname <- gtf_gene$gene_name[match(x=genesid,table = gtf_gene$gene_id,nomatch = 0)]
genesname |> duplicated() |> table() # 2096个重复的基因名
rownames(fpkm) <- genesname


# 表达矩阵去重复（相同基因选择平均值最大的一行）
{
  ### 表达矩阵中取重复的基因名称取平均表达量高的那一行
  expr <- fpkm
  #计算行平均值，按降序排列
  index=order(rowMeans(expr[,]),decreasing = T)
  #调整表达谱的基因顺序
  expr_ordered=expr[index,]
  #对于有重复的基因，保留第一次出现的那个，即行平均值大的那个
  rownames(expr_ordered) |> duplicated() |> table()
  keep=!duplicated(rownames(expr_ordered))
  table(keep)
  #得到最后处理之后的表达谱矩阵
  expr_max=expr_ordered[keep,]
  expr_max[1:3,1:5]
  
}
fpkm <- expr_max
fpkm[1:3,1:3]
dim(fpkm)
save(fpkm,file = "outdata/LUAD_TCGA_fpkm_symbol_maxgene.Rdata")

rm(list = ls())
load('../data/IMMUNE_data/ALL_immu_genesets.Rdata')
load('outdata/LUAD_TCGA_fpkm_symbol_maxgene.Rdata')

