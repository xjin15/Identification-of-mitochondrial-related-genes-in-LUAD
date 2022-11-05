# 干性指数
rm(list = ls())
# 导入文件 干性指数ZMN和cell--------------------------------------------------------------------
p_load(tidyverse,export,ggpubr)
stem <- read.csv("..//maf分组/data/StemnessScore_ZMN.csv") %>% select(-1)
load("outdata/step4_group_data.rds")
stem$Sample <- str_sub(string = stem$Sample,start = 1,end = 16  )
stem <- stem[stem$Annotation == "TP",]
ids <- intersect(allid, stem$Sample)

stem_ZMN_497 <- stem[match(ids, stem$Sample),]
stem_ZMN_497$group <- "ClusterA"
stem_ZMN_497$group[280:497] <- "ClusterB"

p_load(tidyverse)
stem_TCGA <- read.csv("..//maf分组/data/CELL_for_mRNAsi.csv") 
colnames(stem_TCGA)[1] <- "Sample"
stem_TCGA$Sample <- str_sub(string = stem_TCGA$Sample,start = 1,end = 16  )
stem_LUAD <- stem_TCGA[stem_TCGA$cancer.type == "LUAD" & stem_TCGA$sample.type == 1,]
ids <- intersect(allid, stem_LUAD$Sample)


stem_TCGA_491 <- stem_LUAD[match(ids, stem_LUAD$Sample),]
stem_TCGA_491$group <- "ClusterA"
stem_TCGA_491$group[match(allid[280:497],stem_TCGA_491$Sample)] <- "ClusterB"
stem_TCGA_491$group |> table()
# 画图 ----------------------------------------------------------------------
######赵队干性指数的图##########
library(ggpubr)
ggviolin(x="group",y="stemnessScore", fill = "group",data = stem_ZMN_497, size = 0.1,
         palette = c("#00CCFF","#FF3333"),add = "boxplot", 
         add.params = list(group = "group"))+
  xlab("")+
  ylab("")+
  # scale_y_continuous(limits =c(0,1.5),breaks = seq(0,1.5,0.3))+ 
  stat_compare_means(aes(group=group),label = "p.format",hide.ns = T,
                     bracket.size = 20,label.x = 1.35,label.y = 1.2,size = 8)+
  theme_gray()+
  theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5))+
  labs(title = "RNAsi Plot")+
  theme(plot.title = element_text(size = 20, color = "black", 
                                  hjust = 0.5, 
                                  angle = 0))+
  theme(legend.text = element_text(size = 14, face = "bold"), 
        legend.title = element_text(size = 14, face = "bold")) +
  labs(x = NULL, y = NULL)

library(export)
graph2ppt(file = "output/plots/Stemness.pptx",append = T,)

######## mRNAsi干性指数图##########
library(ggpubr)
ggviolin(x="group",y="mRNAsi", fill = "group",data = stem_TCGA_491, size = 0.1,
         palette = c("#00CCFF","#FF3333"),add = "boxplot", 
         add.params = list(group = "group"))+
  xlab("")+
  ylab("")+
  stat_compare_means(aes(group=group),label = "p.format",hide.ns = T,
                     bracket.size = 20,label.x = 1.35,label.y = 0.6,size = 8)+
  theme_gray()+
  theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5))+
  labs(title = "RNAsi Plot")+
  theme(plot.title = element_text(size = 20, color = "black", 
                                  hjust = 0.5, 
                                  angle = 0))+
  theme(legend.text = element_text(size = 14, face = "bold"), 
        legend.title = element_text(size = 14, face = "bold")) +
  labs(x = NULL, y = NULL)

library(export)
graph2ppt(file = "output/plots/Stemness.pptx",append = T,)

########## TCGA干性指数图EREG.mRNAsi##########
library(ggpubr)
ggviolin(x="group",y="EREG.mRNAsi", fill = "group",data = stem_TCGA_491, size = 0.1,
         palette = c("#00CCFF","#FF3333"),add = "boxplot", 
         add.params = list(group = "group"))+
  xlab("")+
  ylab("")+
  stat_compare_means(aes(group=group),label = "p.format",hide.ns = T,
                     bracket.size = 20,label.x = 1.35,label.y = 1,size = 8)+
  theme_gray()+
  theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5))+
  labs(title = "EREG.mRNAsi Plot")+
  theme(plot.title = element_text(size = 20, color = "black", 
                                  hjust = 0.5, 
                                  angle = 0))+
  theme(legend.text = element_text(size = 14, face = "bold"), 
        legend.title = element_text(size = 14, face = "bold")) +
  labs(x = NULL, y = NULL)

library(export)
graph2ppt(file = "output/plots/Stemness.pptx",append = T,)


# 2.mDNAsi ------------------------------------------------------------------

stem_DNAsi <- read.csv("..//maf分组/data/CELL_for_mDNAsi.csv") 
colnames(stem_DNAsi)[1] <- "Sample"
stem_DNAsi$Sample <- str_sub(string = stem_DNAsi$Sample,start = 1,end = 16  )
stem_DNAsi_LUAD <- stem_DNAsi[stem_DNAsi$cancer.type == "LUAD" & stem_DNAsi$sample.type == 1,]
ids <- intersect(allid, stem_DNAsi_LUAD$Sample)

stem_437d <- stem_DNAsi_LUAD[match(ids, stem_DNAsi_LUAD$Sample),]
stem_437d$group <- "ClusterA"
stem_437d$group[match(allid[280:497],stem_437d$Sample)] <- "ClusterB"
### ###画DNAsi图#########
library(ggpubr)
ggviolin(x="group",y="mDNAsi", fill = "group",data = stem_437d, size = 0.1,
         palette = c("#00CCFF","#FF3333"),add = "boxplot", 
         add.params = list(group = "group"))+
  xlab("")+
  ylab("")+
  # scale_y_continuous(limits =c(0,1.5),breaks = seq(0,1.5,0.3))+ 
  stat_compare_means(aes(group=group),label = "p.format",hide.ns = T,
                     bracket.size = 20,label.x = 1.35,label.y = 0.4,size = 8)+
  theme_gray()+
  theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5))+
  labs(title = "DNAsi Plot")+
  theme(plot.title = element_text(size = 20, color = "black", 
                                  hjust = 0.5, 
                                  angle = 0))+
  theme(legend.text = element_text(size = 14, face = "bold"), 
        legend.title = element_text(size = 14, face = "bold")) +
  labs(x = NULL, y = NULL)

library(export)
graph2ppt(file = "output/plots/Stemness.pptx",append = T,)


