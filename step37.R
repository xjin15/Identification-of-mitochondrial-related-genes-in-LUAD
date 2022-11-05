# 1.读取预测好的耐药性数据1 -----------------------------------------------------------
rm(list = ls())

p_load(tidyverse)
library(data.table)
testPtype <- fread('../mitoGENE/calcPhenotype_Output/DrugPredictions_by_GDSC1.csv', data.table = F) %>% column_to_rownames("V1")
testPtype[1:4, 1:4]
#读取分组数据
load("outdata/step4_group_data.rds")

ptype_longdt <- testPtype %>% 
  rownames_to_column(var = "sample") %>% 
  pivot_longer(cols = 2:368,names_to = "Drug") %>% 
  mutate(Cluster = ifelse(sample %in% sample_group$sample[1:279],"ClusterA","ClusterB")) %>% 
  mutate(Cluster = as.factor(Cluster))


# 2.画图 --------------------------------------------------------------------
data=ptype_longdt
library(export)
# 分组直方图？
library(ggpubr)
# 根据不同药物，比较两组IC50差异
x <- compare_means(value ~ Cluster, data, group.by = "Drug")
xsig <- x[x$p.adj<0.05,]
datasig <- data[data$Drug %in% xsig$Drug,] 


# 肺癌的靶向药物有哪些

LUADdrugs <- c("Gefitinib_1010","Afatinib_1032","Crizotinib_1083","Erlotinib_1168","Trametinib_1372","Osimertinib_1919",
"Savolitinib_1936","Cisplatin_1005","Docetaxel_1007","Docetaxel_1819","Vinorelbine_2048","Paclitaxel_1080",
"Gemcitabine_1190")
data1 <- data[data$Drug %in% LUADdrugs,] %>% 
  rename(name = Drug,
         group = Cluster)

p <- ggboxplot(data1,x = "name",y = "value",color = "group")+ 
  stat_compare_means(aes(group = group),label =  "p.signif", label.x = 1.5) 
p
graph2ppt(file = "output/plots/GDSC.pptx", append = T)

p + facet_wrap(~name,scales = "free")
graph2ppt(file = "output/plots/GDSC.pptx", append = T)

ggviolin(data1,x="name",y="value",fill = "group",size = 0.1,
          palette = c("#00CCFF","#FF3333"))+
  xlab("LUAD Drugs")+ylab("IC50")+
  scale_y_log10()+
  stat_compare_means(aes(group=group),label = "p.signif",hide.ns = T,
                     bracket.size = 10)+
  theme_classic()+
  theme(axis.text = element_text(size = 8), 
        axis.text.x = element_text(size = 8, angle = 45, vjust = 0.6,hjust = 0.5))
  
graph2ppt(file = "output/plots/GDSC.pptx", append = T,width = 12, asp=1.5)
library(ggstatsplot)



ggboxplot(data1,x="name",y="value",fill = "group",size = 0.1,
         palette = c("#00CCFF","#FF3333"))+
  xlab("LUAD Drugs")+ylab("IC50")+
  scale_y_log10()+coord_flip()+
  stat_compare_means(aes(group=group),label = "p.signif",hide.ns = T,
                     bracket.size = 10)+
  theme_classic()+
  theme(axis.text = element_text(size = 8),)


graph2ppt(file = "output/plots/GDSC.pptx",append = T,width = 15, asp = 1.6)

ggviolin(data1[data1$name == "Cisplatin_1005",],x="group",y="value",fill = "group",size = 0.1,
         palette = c("#00CCFF","#FF3333"))+
  xlab("LUAD Drugs")+ylab("IC50")+
  scale_y_log10()+
  stat_compare_means(aes(group=group),label = "p.signif",hide.ns = T,
                     bracket.size = 10)+
  theme_classic()+
  theme(axis.text = element_text(size = 8), 
        axis.text.x = element_text(size = 8, angle = 45, vjust = 0.6,hjust = 0.5))

# 3.读取GDSC2预测数据 -----------------------------------------------------------

p_load(tidyverse)
library(data.table)
testPtype <- fread('calcPhenotype_Output/DrugPredictions_by_GDSC2.csv', data.table = F) %>% column_to_rownames("V1")
testPtype[1:4, 1:4]
#读取分组数据
load("outdata/step4_group_data.rds")

ptype_longdt <- testPtype %>% 
  rownames_to_column(var = "sample") %>% 
  pivot_longer(cols = 2:199,names_to = "Drug") %>% 
  mutate(Cluster = ifelse(sample %in% sample_group$sample[1:279],"ClusterA","ClusterB")) %>% 
  mutate(Cluster = as.factor(Cluster))

data=ptype_longdt
library(export)
# 分组直方图？
library(ggpubr)
# 根据不同药物，比较两组IC50差异
x <- compare_means(value ~ Cluster, data, group.by = "Drug")
xsig <- x[x$p<0.05,]
datasig <- data[data$Drug %in% xsig$Drug,] 
LUADdrugs <- c("Gefitinib_1010","Afatinib_1032","Crizotinib_1083","Erlotinib_1168","Trametinib_1372","Osimertinib_1919",
               "Savolitinib_1936","Cisplatin_1005","Docetaxel_1007","Docetaxel_1819","Vinorelbine_2048","Paclitaxel_1080",
               "Gemcitabine_1190")
data1 <- data[data$Drug %in% LUADdrugs,] %>% 
  rename(name = Drug,
         group = Cluster)
# 箱线图

ggboxplot(data1,x="name",y="value",fill = "group",size = 0.1,
          palette = c("#00CCFF","#FF3333"))+
  xlab("LUAD Drugs")+ylab("IC50")+
  scale_y_log10()+coord_flip()+
  stat_compare_means(aes(group=group),label = "p.signif",hide.ns = T,
                     bracket.size = 10)+
  theme_classic()+
  theme(axis.text = element_text(size = 8),)

graph2ppt(file = "output/plots/GDSC.pptx",append = T,width = 15, asp = 1.6)

# 对顺铂的IC50 预后好IC50反而更高，更耐药？
ggviolin(data1[data1$name == "Cisplatin_1005",],x="group",y="value",fill = "group",size = 0.1,
         palette = c("#00CCFF","#FF3333"),add = "boxplot",add.params = list(fill = "white"))+
  xlab("LUAD Drugs")+ylab("IC50")+
  scale_y_log10()+
  stat_compare_means(aes(group=group),label = "p.signif",hide.ns = T,
                     bracket.size = 10)+
  theme_classic()+
  theme(axis.text = element_text(size = 8), 
        axis.text.x = element_text(size = 8, angle = 45, vjust = 0.6,hjust = 0.5))


# Osimertinib_1919 奥西替尼
ggviolin(data1[data1$name == "Osimertinib_1919",],x="group",y="value",fill = "group",size = 0.1,
 palette = c("#00CCFF","#FF3333"),add = "boxplot",add.params = list(fill = "white"))+
  # xlab("LUAD Drugs")+
  ylab("IC50")+
  # scale_y_log10()+
  stat_compare_means(aes(group=group),hide.ns = T,
                     label = "p.signif",label.x.npc = "centre",
                     bracket.size = 10)+
  theme_classic()+
  theme(axis.text = element_text(size = 8), 
        axis.text.x = element_text(size = 8, angle = 45, vjust = 0.6,hjust = 0.5))


# 4.挑选有差异的，并且A的IC50低的药物 ---------------------------------------------------

dataA <- data[data$Cluster == "ClusterA",]
dataA_wide <- pivot_wider(dataA,names_from = Drug,values_from = value)
medianA <- apply(dataA_wide[,3:ncol(dataA_wide)],MARGIN = 2,median)

dataB <- data[data$Cluster == "ClusterB",]
dataB_wide <- pivot_wider(dataB,names_from = Drug,values_from = value)
medianB <- apply(dataB_wide[,3:ncol(dataB_wide)],MARGIN = 2,median)

setequal(names(medianA),names(medianB))

DrugAless <- names(medianA)[medianA < medianB]
drugmedian <- data.frame(clusterA = medianA,
                         clusterB = medianB) %>% 
  rownames_to_column("drugs")

drugmedian_less <- drugmedian[drugmedian$drugs %in% DrugAless,]
intersect(drugmedian_less$drugs,LUADdrugs)
data_plot <- datasig[datasig$Drug %in% DrugAless,] %>% 
  rename(name = Drug,
         group = Cluster)

ggboxplot(data_plot[data_plot$name %in% LUADdrugs,],x="name",y="value",fill = "group",size = 0.1,
          palette = c("#00CCFF","#FF3333"))+
  xlab("Drugs")+ylab("IC50")+
  scale_y_log10()+
  coord_flip()+
  stat_compare_means(aes(group=group),label = "p.signif",hide.ns = T,
                     bracket.size = 10)+
  theme_classic()+
  theme(axis.text = element_text(size = 8),)


ggviolin(data_plot[data_plot$name %in% LUADdrugs,],x="group",y="value",fill = "group",size = 0.1,
         palette = c("#00CCFF","#FF3333"),add = "boxplot",add.params = list(fill = "white"),
         facet.by = "name")+
  # xlab("LUAD Drugs")+
  ylab("IC50")+
  # scale_y_log10()+
  stat_compare_means(aes(group=group),hide.ns = T,
                     label = "p.signif",label.x.npc = "centre",
                     bracket.size = 10)+
  theme_classic()+
  theme(axis.text = element_text(size = 8), 
        axis.text.x = element_text(size = 8, angle = 45, vjust = 0.6,hjust = 0.5))

graph2ppt(file = "output/plots/GDSC.pptx",append = T,width = 15, asp = 1.6)
