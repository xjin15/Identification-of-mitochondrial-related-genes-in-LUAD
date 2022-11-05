# 预测药物疗效

# 1.认识gdsc数据库的表达量矩阵和药物的IC50 -----------------------------------------------
rm(list = ls()) 
## 魔幻操作，一键清空~
options(stringsAsFactors = F)
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
dir='../data/GDSC/Training Data/'
GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
dim(GDSC2_Expr)  # 17419   805    805中细胞系的17419个基因的表达矩阵
GDSC2_Expr[1:4, 1:4]
boxplot(GDSC2_Expr[,1:4])
df=melt(GDSC2_Expr[,1:4])
head(df)
p1=ggboxplot(df, "Var2", "value") +th
p1
# Read GDSC2 response data. rownames() are samples, colnames() are drugs. 
GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
dim(GDSC2_Res)#805 198 805种细胞系，198个药物
GDSC2_Res[1:4, 1:4]
p2=ggboxplot(melt(GDSC2_Res[ , 1:4]), "Var2", "value") +th ; p2
# IMPORTANT note: here I do e^IC50 since the IC50s are actual ln values/log transformed already,
# and the calcPhenotype function Paul #has will do a power transformation (I assumed it would be better to not have both transformations)
# RES值实际上是IC50值的log或者ln化以后的数值
GDSC2_Res <- exp(GDSC2_Res) 
p3=ggboxplot(melt(GDSC2_Res[ , 1:4]), "Var2", "value") +th ; p3
library(patchwork)
p1+p2+p3


# 就是看针对具体的细胞系来说，那些药物有奇效那些药物是打酱油。
# 代码如下所示：
ggboxplot(melt(GDSC2_Res[ 1:4 ,]), "Var1", "value") +th
# 箱线图不知道具体的药物
round(apply(GDSC2_Res[ 1:4 ,], 1, function(x){
  return(c(
    head(sort(x)),
    tail(sort(x))
  ))
}),2)

apply(GDSC2_Res[ 1:4 ,], 1, function(x){ 
  names(x)=gsub('_[0-9]*','',colnames(GDSC2_Res))
  return(c(
    names(head(sort(x))),
    names(tail(sort(x)))
  ))
})


# 2.使用oncoPredict ---------------------------------------------------------
rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
p_load(oncoPredict)
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
dir='../data/GDSC/Training Data/'
GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res) 
####### 读入自己的表达矩阵
load("outdata/LUAD_TCGA_fpkm_symbol_maxgene.Rdata")
load("outdata/step4_group_data.rds")
dim(fpkm)
fpkm[1:4,1:4]
rownames(GDSC2_Expr) %>% length() #17419
rownames(fpkm)%>% length() #58387
intersect(rownames(GDSC2_Expr), rownames(fpkm)) %>% length() #17133
setdiff(rownames(GDSC2_Expr), rownames(fpkm))
# fpkm矩阵中比GDSC少了很多200多个基因，主要是一些lncRNA
### 读取gtf文件按
# gtf_gene <- data.table::fread(file = "../maf分组/output/gtfmRNA22.txt", header = T,sep = "\t")
# intersect(rownames(GDSC2_Expr), gtf_gene$gene_name) %>% length()
#17133 
# 说明这个v22的注释文件里面就是少了一部分基因。
# 试试最新的v38注释版本
# AnnoData = rtracklayer::import('../data/gencode.v38.annotation.gtf.gz')
# intersect(rownames(GDSC2_Expr), AnnoData$gene_name) %>% length() 
#16302

## 没办法，只能用v22的注释版本，比训练集确实少了200多个基因。##
gnames <- intersect(rownames(GDSC1_Expr),rownames(fpkm))
fpkm_luad <- fpkm[gnames,sample_group$sample]
testExpr <- fpkm_luad
GDSC1_Expr = readRDS(file=file.path(dir,'GDSC1_Expr (RMA Normalized and Log Transformed).rds'))
GDSC1_Res = readRDS(file = file.path(dir,"GDSC1_Res.rds"))
GDSC1_Res <- exp(GDSC1_Res) 
save(GDSC2_Expr,GDSC2_Res,testExpr,file = "outdata/step15GDSC2predict.Rdata")
save(GDSC1_Expr,GDSC1_Res,testExpr,file = "outdata/step15GDSC1predict.Rdata")
load(file = "outdata/step15GDSC1predict.Rdata")
###### 一个函数预测药物,注意：运行时间需要半小时以上,本质上是岭回归

calcPhenotype(trainingExprData = GDSC1_Expr,
              trainingPtype = GDSC1_Res,
              testExprData = testExpr,
              batchCorrect = 'eb',  #   "eb" for ComBat
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10,
              printOutput = TRUE,
              removeLowVaringGenesFrom = 'rawData' )

# 前面的R包 oncoPredict的核心函数calcPhenotype运行完毕后，
# 会在当前工作目录下面输出 calcPhenotype_Output 文件夹，
# 里面有一个 DrugPredictions.csv的文件，
library(data.table)
testPtype <- fread('calcPhenotype_Output/DrugPredictions.csv', data.table = F)
testPtype[1:4, 1:4]

# 3.pRphonetic包 -----------------------------------------------------------


