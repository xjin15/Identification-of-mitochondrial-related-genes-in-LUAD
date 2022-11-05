# GSVA分析24种细胞比例
################ GSVA 分析 ##########
##### 1数据准备+GSVA  注意：TCGAgs必须是的变量必须是因子######
# 因子很容易转换成任何其他数据类型
# 但是其他的数据类型，比如字符型无法直接转换成因子
rm(list = ls())
load(file = "outdata/step4.2_mrna.rds")
load(file = "outdata/LUAD_TCGA_ALLfpkm.rds")
# 直接使用mrna的表达矩阵就可以预测免疫细胞的比例了。如果也可以试试FPKM和mRNA预测结果的差异

library(tidyverse)
library(GSVA)
library(pheatmap)

TCGAgs <- read.csv("maf分组/outdata/TCGAgs.csv", stringsAsFactors = F, 
                   row.names = 1,na.strings = "")
###
TCGAgs<-as.list(TCGAgs)
TCGAgs<-lapply(TCGAgs, function(x) x[!is.na(x)])#去除list中的NA
# 所谓表达矩阵，矩阵矩阵要转换为matrix
fpkm <- as.matrix(fpkm)
rownames(fpkm)
genesid <- rownames(fpkm)
genesid <- substring(genesid,1,15)
genesname <- gtf_gene$gene_name[match(x=genesid,table = gtf_gene$gene_id,nomatch = 0)]
genesname |> duplicated() |> table() # 2096个重复的基因名

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

fpkm_unique <- expr_max
fpkm_unique_497 <- fpkm_unique[,allid497]
# 得到GSVA得分
res_es_gsva <- gsva(fpkm_unique_497, TCGAgs, 
                    min.sz = 1,
                    max.sz = Inf,
                    mx.diff=T, #mx.diff=FALSE es值是一个双峰的分布T代表近似正态分布
                    verbose=FALSE, 
                    parallel.sz=0,
                    method= "gsva")
# 先做个热图看看
pheatmap(res_es_gsva, cluster_cols = F,cluster_rows = T,show_colnames = F,
         scale = "none")
graph2ppt(file = "output/plots/GSVA.pptx",append = T)
res_es_gsva1 <- res_es_gsva # fpkm
allgsva <- res_es_gsva
allgsva <- as.matrix(allgsva)
# 保存为rds文件
# save(allgsva,allid497,file = "outdata/step10_GSVA24cell.rds")

  