
############################# 1 数据准备##########################
rm(list = ls())
# GenVisR
p_load(maftools, ggplot2, reshape2, data.table, ggpubr, export)
load(file = "outdata/step4_group_data.rds")
load(file = "outdata/step8_LUADmaf.rds")  # 加载我们已经分好组的maf文件
#读入TCGA-LUAD.maf
# 这里有个isTCGA选项，选择T，它会自动精简我们的样本barcode成前12位
# subsetmaf 给maf文件按照我们自己的样本取子集的maf
{# subsetmaf 把原来的TCGA-LUADmaf文件根据TSB取子集。
# TCGAdataMAF<-read.maf("../maf分组/data/LUAD.maf",isTCGA = T)

# tsb_A <- sample_group$sample[sample_group$group == 1] |> substring(1,12)
# tsb_B <- sample_group$sample[sample_group$group == 2] |> substring(1,12)
# 
# ClusterA_maf <- subsetMaf(maf = TCGAdataMAF,tsb = tsb_A)
# ClusterB_maf <- subsetMaf(maf = TCGAdataMAF,tsb = tsb_B)
# save(TCGAdataMAF,ClusterA_maf,ClusterB_maf,file = "outdata/step8_LUADmaf.rds")
# 
  }



# 显示整个maf文件中突变基因之间的相关性
somaticInteractions(maf = TCGAdataMAF, top = 25, 
                    pvalue = c(0.05, 0.01)  # pvalue表示星号和＃ 的阈值
                    
                    )

graph2ppt(file="output/plots/SNP.pptx",append = T)
# 肿瘤通路
OncogenicPathways(maf = TCGAdataMAF)
graph2ppt(file="output/plots/SNP.pptx",append = T)



# 2突变组oncoplot等 -----------------------------------------------------------

oncoplot(maf = ClusterA_maf, top = 20, )
graph2ppt(file="output/plots/SNP.pptx",append = T)

somaticInteractions(maf = ClusterA_maf,)
graph2ppt(file="output/plots/SNP.pptx",append = T)
graph2pdf(file="output/plots/oncoplot_Amut.pdf")

OncogenicPathways(maf = ClusterA_maf)
graph2ppt(file="output/plots/SNP.pptx",append = T)



# 3 B组瀑布图等 ----------------------------------------------------------------



oncoplot(maf = ClusterB_maf, top = 20, )
graph2ppt(file="output/plots/SNP.pptx",append = T)

somaticInteractions(maf = ClusterB_maf,)
graph2ppt(file="output/plots/SNP.pptx",append = T)
graph2pdf(file="output/plots/oncoplot_BWild.pdf",)

OncogenicPathways(maf = ClusterB_maf)
graph2ppt(file="output/plots/SNP.pptx",append = T)



# 4 比较AB两组的突变基因 -----------------------------------------------------------


compareMAF_AB <- mafCompare(m1 = ClusterA_maf, m2 = ClusterB_maf, 
                            m1Name = 'Cluster_A', m2Name = 'Cluster_B', 
                            minMut = 10,
                            pathways = T, # 如果设置了pathway=T，就会对比两组基因突变在一些信号通路上的差异
                            useCNV = F, # 设置了PATHWAY=T ,CNV就不能设置了
)#最小mut数默认为5的基因
# 展示两组差异突变基因所在的通路
forestPlot(mafCompareRes = compareMAF_AB, 
           # color = c('royalblue', 'maroon'), 
           geneFontSize = 0.8)
graph2ppt(file = "output/plots/SNP.pptx",append = T,) # 

## 设置了fdr以后，设置P值就没用了 
# 比较差异的基因
compareMAF_AB_moren <- mafCompare(m1 = ClusterA_maf, m2 = ClusterB_maf, 
                                  m1Name = 'Cluster_A', m2Name = 'Cluster_B', 
                                  minMut = 10,
                                  pathways = F, # 如果设置了pathway=T，就会对比两组基因突变在一些信号通路上的差异
                                  useCNV = T, # 设置了PATHWAY=T ,CNV就不能设置了
)#最小mut数默认为5的基因
forestPlot(mafCompareRes = compareMAF_AB_moren, 
           # color = c('royalblue', 'maroon'), 
           fdr = 0.01,geneFontSize = 0.8)
## 设置了fdr以后，设置P值就没用了 
graph2ppt(file = "output/plots/SNP.pptx",append = T,) 


# graph2pdf(file = "output/plots/野生vs突变8低频基因的cobarplot.pdf")
compareMAF_AB_moren$results

deSNP <- compareMAF_AB_moren$results %>%
  dplyr::arrange(., order(adjPval,decreasing = F)) %>%
  dplyr::filter(is.finite(or)) %>%
  dplyr::select(Hugo_Symbol) %>% head( ., n=25) %>% as.data.frame() %>% .[,1]
deSNP

# 比较两组差异的突变图
coBarplot(genes = deSNP[1:8], #设置要展示的基因，默认是突变前五的，但这样就不是突变差异明显的那些的，所以还是要自己设置、
          m1 = ClusterA_maf, 
          m2 = ClusterB_maf, 
          m1Name = "Cluster_A", 
          m2Name = "Cluster_B",)
graph2ppt(file = "output/plots/SNP.pptx",append = T,) 


# print(compareMAF_AB)
# summary_compareMAF<-as.data.frame(compareMAF_AB$results)

# 5 加上cnv的瀑布图 -------------------------------------------------------------
# TCGA肺腺癌的CNV整合数据已经做好。直接读取。GISTIC 分数。1代表扩增-1代表缺失0不变
TCGAcnvmelt<-read.csv("outdata/TCGAcnvmelt.csv",row.names = 1)
#整合入data
class(TCGAcnvmelt$Sample_name)
TCGAcnvmelt$Sample_name<-as.character(TCGAcnvmelt$Sample_name)#注意这里要改为文字格式，要不然下面会分成两类！
# TCGAdataMAF<-read.maf("outdata/TCGAdata_mut_MAF_maftools.maf")
TCGAcnvmelt$Sample_name <- substring(TCGAcnvmelt$Sample_name,1,12)

TCGAcnvmeltA <-  TCGAcnvmelt[TCGAcnvmelt$Sample_name %in% ClusterA_maf@data$Tumor_Sample_Barcode,]
TCGAcnvmeltB <-  TCGAcnvmelt[TCGAcnvmelt$Sample_name %in% ClusterB_maf@data$Tumor_Sample_Barcode,]

ClusterA_maf_plusCNV <- read.maf(ClusterA_maf@data, cnTable = TCGAcnvmeltA, isTCGA = T)
ClusterB_maf_plusCNV <- read.maf(ClusterB_maf@data, cnTable = TCGAcnvmeltB, isTCGA = T)

# 注意这里read.maf, 第一个必须是maf文件里面的data，然后cnTabe就用我们前面整理好的CNV文件
# 生成一个新的maf以后，这个新的maf需要在后面write.mafsummary

# 突变总览
oncoplot(maf = ClusterA_maf_plusCNV)
graph2ppt(file = "output/plots/SNP.pptx",append = T,) 
graph2pdf(file = "output/plots/oncoplot_Amut_CNV.pdf")

plotmafSummary(maf = ClusterA_maf_plusCNV, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
graph2ppt(file = "output/plots/SNP.pptx",append = T,) 

oncoplot(maf = ClusterB_maf_plusCNV)
graph2ppt(file = "output/plots/SNP.pptx",append = T,) 
graph2pdf(file = "output/plots/oncoplot_BWild_CNV.pdf")

plotmafSummary(maf = ClusterB_maf_plusCNV, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
graph2ppt(file = "output/plots/SNP.pptx",append = T,) 

# 一个测试:1 和 2 是加了CNV write出去再重新读进来的。和之前的一样，而且文件还小了。
# 因此建议加上CNV以后先write。mafsummary再重新read，这样文件小一些
# 新的测试表明：不能先write再重新读入，这样会让CNV数据没有了。

# plotmafSummary(maf = ClusterA_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
# plotmafSummary(maf = ClusterA_maf_plusCNV, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
# plotmafSummary(maf = ClusterA_maf_plusCNV_1, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)



# genes = c("TP53","TTN","MUC16","CSMD3","RYR2","LRP1B","USH2A","ZFHX4","KRAS","FLG","SPTA1","XIRP2","NAV3","ZNF536","ANK2","CSMD1","FAT3","COL11A1","MUC17","PCDH15"))#显示特定基因 
#最终确定的基因(根据total的数目选择，或者根据与研究相关的基因进行选择)
# oncoplot(maf = TCGAdataMAFA.plus.cn, 
#          keepGeneOrder = T,
#          genes = c("TTN","TP53","CSMD3","RYR2","USH2A","MUC16","ZFHX4","SPTA1",
#                    "LRP1B","NAV3","KRAS","FLG","XIRP2","STK11","COL11A1",
#                    "APOB","TNR","RYR3","PCLO","FAT3","ADGRG4"))#显示特定基因 
# oncoplot(maf = TCGAdataMAFA.plus.cn, 
#          genes = c( "AKR1C1", "AKR1C2", "ZNF724P", "ZNF730", "TMPRSS7", 
#   "ZNF98", "SCN1A", "AKR1C3", "ZNF676", "SMARCA4", "SNX19", 
#   "RPSAP58", "SPICE1", "ZNF681", "CD200R1L", "ABHD10", "ATG3", 
#   "TSPAN16", "FAM181B", "PRSS55", "STK11", "DDIAS", "ZNF728"))




# 6 比较加上CNV的瀑布图 -----------------------------------------------------------
compareMAF_AB_CNV <- mafCompare(m1 = ClusterA_maf_plusCNV, m2 = ClusterB_maf_plusCNV, 
                                m1Name = 'clusterA', m2Name = 'clusterB', minMut = 10, useCNV = T)#最小mut数默认为5的基因
forestPlot(mafCompareRes = compareMAF_AB_CNV, fdr = 0.0001, geneFontSize = 0.8)
# 加上CNV数据以后，比较出来太多差异的基因了。几百上千个，没意义了呀
## 设置了fdr以后，设置P值就没用了
graph2ppt(file = "output/plots/SNP.pptx",append=T) # 记得把KEAP1和 NFE2L2删除 

print(compareMAF_AB)
summary_compareMAF<-as.data.frame(compareMAF_AB$results)
write.csv(summary_compareMAF,"outdata/summary_compareMAF.csv")

# deSNP <- compareMAF_AB$results %>%
#   dplyr::arrange(., order(or,decreasing = T)) %>% 
#   dplyr::filter(is.finite(or)) %>% 
#   dplyr::select(Hugo_Symbol) %>% head( , n=25) %>% as.data.frame() %>% .[,1] 
# deSNP
##得到差异明显的前8个基因，但是这些基因是突变频率比较低的
deSNP <- c("STK11","EGFR","GRIN2B","SPEF2", "SNTG2","BRWD3","OR6N1","ADGRB1")
###这里我确定了10个基因，是前10突变的，看起来有差距的10个基因。
genes <- c("TTN","RYR2","CSMD3","USH2A","SPTA1","ZFHX4","NAV3","EGFR","SPEF2","SNTG2")
deSNP <- genes
## 画出瀑布图
oncoplot(maf = TCGAdataMAFA.plus.cn, 
         keepGeneOrder = T,
         genes = deSNP)
graph2ppt(file = "output/plots/突变组8低频基因的oncoplot.pptx")
graph2pdf(file = "output/plots/突变组8低频基因的oncoplot")

oncoplot(maf = TCGAdataMAFB.plus.cn, 
         keepGeneOrder = T,
         genes = deSNP)#显示特定基因
graph2ppt(file = "output/plots/野生组8低频基因的oncoplot")
graph2pdf(file = "output/plots/野生组8低频基因的oncoplot.pdf")

coOncoplot(genes = deSNP, m1 = TCGAdataMAFA.plus.cn, m2 = TCGAdataMAFB.plus.cn, m1Name = 'MUT', m2Name = 'WILD', removeNonMutated = TRUE)
graph2pdf(file = "output/plots/野生vs突变8低频基因的cooncoplot.pdf")

coBarplot(genes = deSNP, m1 = TCGAdataMAFA.plus.cn, m2 = TCGAdataMAFB.plus.cn, 
          m1Name = 'MUT', m2Name = 'WILD', yLims = c(30,30) )
graph2pdf(file = "output/plots/野生vs突变8低频基因的cobarplot.pdf")

#将比较结果可视化，绘制森林图 keap1肯定明显，因为这是我的分组条件
forestPlot(mafCompareRes = compareMAF_AB, pVal = 0.1, 
           color = c('royalblue', 'maroon'),
           fdr = 0.1,
           geneFontSize = 0.8)














# MAF<-as.data.frame(TCGAdataMAF@data)
# waterfall(MAF, mainRecurCutoff=.17,#[0,1]范围内的数值，基于样本中特定基因的突变频率选择plot的基因。
#           mainXlabel=F,#布尔型，是否显示X轴的sample name；默认为FALSE
#           mainDropMut=F,#布尔型，是否从mutation type legend中去除不使用的”mutation type”，默认为FALSE
#           plotMutBurden = T,#布尔型，是否在顶端展示mutBurden的sub-plot，默认为TRUE
#           #plotSamples = ,#指定plot的样本
#           #maxGenes=100,#指定最多可以plot的基因数
#           rmvSilent = F,#是否去除沉默突变，默认为FALSE
#           fileType = "MAF",#读入X的数据格式，分为：”MGI”, “MAF”, “Custom”
#           out = "plot") #输出的类型,包含 “data”, “grob”, 和”plot”三种,默认为 “plot”
# # graph2ppt(file="output/plots/野生组瀑布图.pptx")
# graph2pdf(file="output/plots/野生组瀑布图2", width = 12, aspectr = 0.7)

# MAF <- as.data.frame(TCGAdataMAF@data)
# waterfall(MAF, mainRecurCutoff=.23,#[0,1]范围内的数值，基于样本中特定基因的突变频率选择plot的基因。
#           mainXlabel=F,#布尔型，是否显示X轴的sample name；默认为FALSE
#           mainDropMut=F,#布尔型，是否从mutation type legend中去除不使用的”mutation type”，默认为FALSE
#           plotMutBurden = T,#布尔型，是否在顶端展示mutBurden的sub-plot，默认为TRUE
#           #plotSamples = ,#指定plot的样本
#           #maxGenes=100,#指定最多可以plot的基因数
#           rmvSilent = F,#是否去除沉默突变，默认为FALSE
#           fileType = "MAF",#读入X的数据格式，分为：”MGI”, “MAF”, “Custom”
#           out = "plot") #输出的类型,包含 “data”, “grob”, 和”plot”三种,默认为 “plot”
# ## 也可以选择自己想展示的基因来画图
# waterfall(MAF, #mainRecurCutoff=.17,#[0,1]范围内的数值，基于样本中特定基因的突变频率选择plot的基因。
#           mainXlabel=F,#布尔型，是否显示X轴的sample name；默认为FALSE
#           mainDropMut=F,#布尔型，是否从mutation type legend中去除不使用的”mutation type”，默认为FALSE
#           plotMutBurden = T,#布尔型，是否在顶端展示mutBurden的sub-plot，默认为TRUE
#           #plotSamples = ,#指定plot的样本
#           #maxGenes=100,#指定最多可以plot的基因数
#           plotGenes = c("TTN","TP53","CSMD3", "USH2A","MUC16","ZFHX4","SPTA1",
#                         "LRP1B","NAV3","KRAS","FLG","XIRP2","STK11","COL11A1",
#                         "APOB","TNR","RYR3","PCLO","FAT3","ADGRG4"),
#           # geneOrder = c("TTN","TP53","CSMD3", "USH2A","MUC16","ZFHX4","SPTA1",
#           #               "LRP1B","NAV3","KRAS","FLG","XIRP2","STK11","COL11A1",
#           #               "APOB","TNR","RYR3","PCLO","FAT3","ADGRG4"),
#           rmvSilent = F,#是否去除沉默突变，默认为FALSE
#           fileType = "MAF",#读入X的数据格式，分为：”MGI”, “MAF”, “Custom”
#           out = "plot") #输出的类型,包含 “data”, “grob”, 和”plot”三种,默认为 “plot” 
# # graph2ppt(file="output/plots/突变组瀑布图.pptx")
# graph2pdf(file="output/plots/突变组瀑布图2", width = 12, aspectr = 0.7)


