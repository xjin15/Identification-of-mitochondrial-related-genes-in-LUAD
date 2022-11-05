rm(list = ls())
# 读入差异基因和临床数据
deg_mito <- read.table(file = "outdata/deg_mito.tsv") |> as.vector() |> unlist()
load(file = "outdata/step4.2_mrna.rds")
load(file = "outdata/LUAD_cln_clean.rds")

# 选择有mrna矩阵的差异基因
meta <- phe_f |> select(sample,OS,OS.time) |> 
  rename(event = OS,
         time = OS.time) 
meta <- meta[sample %in% colnames(mrna),]
mrna_all <- as.data.frame(mrna)
mrna_degmito <- mrna_all[deg_mito,]
exprSet_hub <- mrna_degmito
library(tinyarray) 
# 花花老师的R包，可以批量算截断值和生存分析。好用。
# 如果有问题就看批量分析，那边有解决报错的方法
# 批量计算截断值

point_cut(exprSet_hub,meta)

p_KMsig <- surv_KM(exprSet_hub,
             meta,
             cut.point = T, # 批量计算截断值
             pvalue_cutoff = 0.9,
             min_gn = 0.1) %>% as.data.frame()


p_COXsig_deg1747 <- surv_cox(exprSet = exprSet_hub,
                     meta = meta,
                     cut.point = T,
                     pvalue_cutoff = 1,
                     HRkeep = "all",
                     continuous = F,
                     min_gn = 0.1)

exprSet_hub <- mrna_all[rownames(DEG_more),]
p_COXsig_deg1747 <- surv_cox(exprSet = exprSet_hub,
                     meta = meta,
                     cut.point = T,
                     pvalue_cutoff = 1,
                     HRkeep = "all",
                     continuous = F,
                     min_gn = 0.1)

write.table(p_COXsig,file = "outdata/deg_more_cox_analysis.tsv",
            sep = "\t",row.names = T)


