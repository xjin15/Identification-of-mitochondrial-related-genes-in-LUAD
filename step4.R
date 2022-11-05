rt1_TP <- rt1[,1:510]
# 肿瘤组织的fpkm值

DEMTgene <- c("MT.ND1", "MT.ND5", "MT.TV", "MT.CYB", "MT.TL2", "MT.ND3", 
  "MT.TS2", "MT.TM", "MT.TT", "MT.TF", "MT.TN", "MT.TE", "MT.TR", 
  "MT.ND6")
options(scipen = 10)
# 一致性聚类 -------------------------------------------------------------------
rt1_TP
mads=apply(rt1_TP,1,mad) 
means=rowMeans(rt1_TP)
means
library(ALL)
data(ALL)
d=exprs(ALL)
dim(d)
d[1:5,1:5]  #共128个样品，12625个探针数据
##对这个芯片表达数据进行简单的normalization，
# 取在各个样品差异很大的那些gene或者探针的数据来进行聚类分析
mads=apply(d,1,mad)   #计算每个基因的标准差
# mad函数是用来估计标准差的。
# 标准差 = R语言中的mad函数，得到的是估计标准差。
# 绝对中位差实际求法是用原数据减去中位数后得到的新数据的绝对值的中位数。
# 但绝对中位差常用来估计标准差，估计标准差=1.4826*绝对中位差。
# R语言中返回的是估计的标准差。
mads <- mads[mads>0]
names(mads)
mads <- mads[intersect(names(mads),DEMTgene)]
names(mads)

#d=d[rev(order(mads)),] # 取差异最大的前5000个基因
rt1_TP <- rt1_TP[rev(order(mads)),] %>% as.matrix()
#sweep函数减去中位数进行标准化
# d = sweep(d,1, apply(d,1,median,na.rm=T))
dd = sweep(rt1_TP,1, apply(rt1_TP,1,median,na.rm=T))

library(scales)
dd <- t(apply(rt1_TP,1,rescale))
#也可以对这个d矩阵用DESeq的normalization 进行归一化，取决于具体情况
### 一致性聚类
library(ConsensusClusterPlus)
title="output/MITO"  #设置图片输出路径

#结果将会输出k从2-6各个情况下的分型情况，聚类的方法用的是 hc ，抽样比例为0.8，最后输出png图
results = ConsensusClusterPlus(dd,maxK=6,
                               reps=1000,pItem=0.8, # 抽样比例为0.8
                               pFeature=1,title=title,clusterAlg="hc", # 聚类的方法用的是 hc 
                               distance="pearson",
                               plot="png")
#这里设置的maxK=6、reps=50，但是实际上需要更高的reps（如1000）和更高的maxK（如20）

###查看结果
#ConsensusClusterPlus输出的是一个列表，其中列表对应于来自Kth集群的结果，例如，results[[2]]是k=2的结果。
View(results)

#输出K=2时的一致性矩阵，consensusMatrix输出一致矩阵。
results[[2]][["consensusMatrix"]][1:5,1:5]

#           [,1]      [,2]      [,3]      [,4]      [,5]
# [1,] 1.0000000 1.0000000 0.9655172 1.0000000 1.0000000
# [2,] 1.0000000 1.0000000 0.8857143 1.0000000 1.0000000
# [3,] 0.9655172 0.8857143 1.0000000 0.9166667 0.8823529
# [4,] 1.0000000 1.0000000 0.9166667 1.0000000 1.0000000
# [5,] 1.0000000 1.0000000 0.8823529 1.0000000 1.0000000

# hclust选项
results[[2]][["consensusTree"]]

# Call:
# hclust(d = as.dist(1 - fm), method = finalLinkage)
# Cluster method   : average
# Number of objects: 128

#consensusClass可以看到各个样品被分到了哪个类别里面去
results[[2]][["consensusClass"]][1:5]

# 01005 01010 03002 04006 04007
#    1     1     1     1     1
##### 计算聚类一致性 (cluster-consensus) 和样品一致性 (item-consensus)#####
icl <- calcICL(results, title = title,plot = "png")


## 返回了具有两个元素的list，分别查看一下
dim(icl[["clusterConsensus"]])
# [1] 20  3

icl[["clusterConsensus"]]

#       k cluster clusterConsensus
# [1,] 2       1        0.9079483
# [2,] 2       2        0.7584326
# [3,] 3       1        0.6246200
# [4,] 3       2        0.9111359
# [5,] 3       3        0.9864123
# [6,] 4       1        0.8908356
# [7,] 4       2        0.8869606
# [8,] 4       3        0.6663949
# [9,] 4       4        0.9829523
# [10,] 5       1        0.8612347
# [11,] 5       2        0.8848722
# [12,] 5       3        0.5568284
# [13,] 5       4        0.8390983
# [14,] 5       5        1.0000000
# [15,] 6       1        0.8256498
# [16,] 6       2        0.9377737
# [17,] 6       3        0.6496445
# [18,] 6       4        0.7267928
# [19,] 6       5        0.6982017
# [20,] 6       6        1.0000000

dim(icl[["itemConsensus"]])   #128*（2+3+4+5+6）=2560
# [1] 2560    4

icl[["itemConsensus"]][1:5,]

#   k cluster  item itemConsensus
# 1 2       1 01007    0.05261929
# 2 2       1 01003    0.05551604
# 3 2       1 02020    0.04554248
# 4 2       1 04018    0.06059102
# 5 2       1 09002    0.06779347


# 分成2组 --------------------------------------------------------------------


groupmito <- results[[2]][["consensusClass"]]
cln <- read.csv(file = "data/LUAD_survival.csv")
cln <- cln[cln$sample %in% names(groupmito),]
cln$group <- groupmito[base::match(cln$sample, names(groupmito),nomatch = 0)]
cln$group %>% table()

# 检查group分组信息情况
table(cln$group)
# group必须是因子，代表有无突变
str(cln)
cln$group <- as.factor(cln$group)
# write.csv(cln, file = "outdata/cln_xena.csv", row.names = F)
################## 读入整理好的生存数据，做KM图###############
# rm(list = ls())
# cln <- read.csv("outdata/cln_xena.csv")
# 以下代码直接运行即可
cln$OS.time <- cln$OS.time / 30

# data1 <- cln %>% filter(OS.time<=62)
data1 <- cln
p_load(survival,survminer)
fit<-survfit(Surv(OS.time,OS)~group, data = data1) 

surv_summary(fit) #查看生存率及其标准误
surv_median(fit = fit,combine = F) # 查看中位生存时间
pp <- ggsurvplot(fit, data=data1, linetype = 1, palette = c("jco"),size=1,
                 surv.scale = c("percent"),pval = TRUE,
                 legend.title = '', legend.labs=c("Mut","Wild"),
                 break.time.by =12,
                 xlim = c(0,60),risk.table = TRUE,
                 risk.table.title = "Patients at risk",
                 ylab = "Overall Survival, %",xlab = "Months",
                 font.x = c(20,"plain","black"),font.y = c(20,"plain","black"),
                 font.tickslab = c(20,"plain","black"),
                 risk.table.fontsize = 6,font.legend =c(20,"plain","black"),
                 font.main = c(20,"plain","black"),pval.size = 8)

######risk table 的修改
pp$table <- ggpar(
  pp$table,
  font.title    = c(13, "bold.italic", "green"),  ###risk table 标题的修改
  font.subtitle = c(15, "bold", "pink"),  ###risk table小标题的修改
  font.caption  = c(11, "plain", "darkgreen"), ####插入字的修改
  font.x        = c(18, "plain", "black"), ### risk table x的修改
  font.y        = c(18, "plain", "black"),### risk table y的修改
  font.xtickslab = c(18, "plain", "black"),### risk table x 坐标轴的修改
  font.ytickslab = c(18),  ### risk table y 坐标轴的修改
  legend=c(0.8,0.88),
  censor.size=3
)
print(pp)
graph2ppt(file="output/plots/kmplot_os_by_xena_survival.pptx")
