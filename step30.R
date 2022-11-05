# 干性指数数据准备
rm(list = ls())
# 赵队数据
stem_ZMN <- read.csv("..//maf分组/data/StemnessScore_ZMN.csv") %>% select(-1)
stem_ZMN$Sample <- str_sub(string = stem_ZMN$Sample,start = 1,end = 16  )

# CELL文章mRNAsi数据
stem_mRNAsi <- read.csv("..//maf分组/data/CELL_for_mRNAsi.csv") 
colnames(stem_mRNAsi)[1] <- "Sample"
stem_mRNAsi$Sample <- str_sub(string = stem_mRNAsi$Sample,start = 1,end = 16  )
# CELL 文章mDNAsi数据
stem_mDNAsi <- read.csv("..//maf分组/data/CELL_for_mDNAsi.csv") 
colnames(stem_mDNAsi)[1] <- "Sample"
stem_mDNAsi$Sample <- str_sub(string = stem_mDNAsi$Sample,start = 1,end = 16  )

save(stem_ZMN,stem_mDNAsi,stem_mRNAsi, file = "outdata/step11_stemness.rds")
