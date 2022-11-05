########## cnv分析，最后得到差不多是cnv的箱线图或者小提琴图，和snv分析流程相似
####### 1.下载CNV的score数据 ######
rm(list = ls())
p_load(dplyr,maftools,export)
load("outdata/step4_group_data.rds")
load("outdata/step8_LUADmaf.rds")

####################    4.写入mut组的maf并做cnv分析   #############
TCGAcnvmelt<-read.csv("../mitoGENE/outdata/TCGAcnvmelt.csv",row.names = 1)
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


#展示重点变量的总结信息
# #Shows sample summry.
# getSampleSummary(ClusterA_maf_plusCNV)
# SampleSummaryA<-getSampleSummary(ClusterA_maf_plusCNV)
# #Shows gene summary.
# getGeneSummary(ClusterA_maf_plusCNV)
# GeneSummary<-getGeneSummary(ClusterA_maf_plusCNV)
# #Shows all fields in MAF
# getFields(ClusterA_maf_plusCNV)
# Fields<-getFields(ClusterA_maf_plusCNV)
# #shows clinical data associated with samples
# getClinicalData(ClusterA_maf_plusCNV)
# ClinicalData<-getClinicalData(ClusterA_maf_plusCNV)
# #Writes maf summary to an output file with basename laml.
write.mafSummary(maf = ClusterA_maf_plusCNV, basename = 'outdata/ClusterA_maf_plusCNV')
# plotmafSummary


# ClusterA_maf_plusCNV_1<-read.maf(maf = 'outdata/ClusterA_maf_plusCNV_maftools.maf')
# plotmafSummary(maf = ClusterA_maf_plusCNV, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

# 一个测试:1 和 2 是加了CNV write出去再重新读进来的。和之前的一样，而且文件还小了。
# 因此建议加上CNV以后先write。mafsummary再重新read，这样文件小一些
# 新的测试表明：不能先write再重新读入，这样会让CNV数据没有了。
 
 # plotmafSummary(maf = ClusterA_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
 # plotmafSummary(maf = ClusterA_maf_plusCNV, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
 # plotmafSummary(maf = ClusterA_maf_plusCNV_1, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)


#瀑布图

oncoplot(maf = ClusterA_maf_plusCNV_1, top = 20,
         # genesToIgnore = "KEAP1",
         )#,
oncoplot(maf = ClusterA_maf_plusCNV, top = 20,
         # genesToIgnore = "KEAP1",
         )#,
graph2ppt(file = "output/plots/SNP.pptx",append = T)
graph2pdf(file = "output/plots/oncoplot_Amut_CNV.pdf")


oncoplot(maf = ClusterB_maf_plusCNV, top = 20,
         # genesToIgnore = "KEAP1",
         )#,
graph2ppt(file = "output/plots/SNP.pptx",append = T)
graph2pdf(file = "output/plots/oncoplot_BWild_CNV.pdf")

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
graph2ppt(file = "output/plots/突变组加上cnv的oncoplot.pptx")
graph2pdf(file = "output/plots/突变组加上cnv的oncoplot.pdf")

################ 两个MAF相互比较 ############
compareMAF_AB_CNV <- mafCompare(m1 = ClusterA_maf_plusCNV, m2 = ClusterB_maf_plusCNV, 
                            m1Name = 'clusterA', m2Name = 'clusterB', minMut = 10, useCNV = T)#最小mut数默认为5的基因
forestPlot(mafCompareRes = compareMAF_AB_CNV, fdr = 0.0001, geneFontSize = 0.8)
# 加上CNV数据以后，比较出来太多差异的基因了。几百上千个，没意义了呀
## 设置了fdr以后，设置P值就没用了
graph2ppt(file = "output/plots/SNP.pptx") # 记得把KEAP1和 NFE2L2删除 

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

