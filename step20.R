rm(list=ls())
pacman::p_load(tidyverse)
load(file = "outdata/step4_group_data.rds")
# BiocManager::install("PoisonAlien/TCGAmutations")
# TCGAmutations包整合了TCGA中全部样本的maf文件
# devtools::install_github(repo = "PoisonAlien/TCGAmutations")
p_load(TCGAmutations)
tmp=as.data.frame(tcga_available())
dt <- TCGAmutations::tcga_load(study = "LUAD")
dt <- dt@data
dt1 <- as.data.frame(table(dt$Tumor_Sample_Barcode))
names(dt1) <- c('Barcode', 'Freq')
dt1$tmb <- dt1$Freq/38
names(dt1)
write.csv(dt1, file = 'outdata/LUAD_TMB.csv')

rt <- dt1 |> mutate(sample = as.character(Barcode))
rt$Barcode <- NULL
rt$sample <- substring(rt$sample,first = 1,last = 16)

rt1 <- rt[rt$sample %in% sample_group$sample, ] 
rt1$group <- sample_group$group[match(rt1$sample, sample_group$sample,nomatch = 0)]

TCGAtmb <- rt1 |> select(group, tmb, sample) 
rownames(TCGAtmb) <- TCGAtmb$sample
TCGAtmb$sample <- NULL

library(tidyverse)
#变成长格式
TCGAtmbbox <- reshape2::melt(data = TCGAtmb,id.vars=c("group"))
TCGAtmbbox$group <- as.factor(TCGAtmbbox$group)
TCGAtmbbox %>% filter(group == 1) %>% summarise(mean = mean(value))
mutval <- TCGAtmbbox %>% filter(group == 1) %>% select(value)  
wildval <- TCGAtmbbox %>% filter(group == 2) %>% select(value)  
apply(mutval, 2, median)
apply(wildval, 2, median)
#箱线图绘制
ggplot(TCGAtmbbox)+  
  geom_boxplot(aes(x= variable,y= value,fill=group),
               width = 0.6,        #宽度
               position = position_dodge(0.8),
               outlier.size = 1, outlier.color = "black"#箱外点大小颜色
  )+
  scale_fill_manual(values = c("#FF3333","#00CCFF"),
                    breaks = c("mut","wild"),#分组的东西
                    labels = c("Mut","Wild"))+
  xlab("Mut  VS.  Wild")+
  ylab("Total Mutation Burden")+
  theme(axis.text.x = element_text(size = 15, hjust = 0.5, vjust = 0.5))# +
  #scale_y_continuous(limits =c(0, 150), breaks = seq(0,150,50))


####小提琴图绘制
#  scale_y_continuous(limits =c(0,1600),breaks = seq(0,1600,500))+ 
library(ggpubr)
p <- ggpubr::ggviolin(x="variable",y="value",fill = "group",data = TCGAtmbbox,size = 0.1,
              palette = c("#00CCFF","#FF3333"),add = "boxplot",position = position_dodge(1), 
              add.params = list(group = "group"))+
  xlab("")+
  ylab("")+
  #scale_y_continuous(limits =c(-200, 1300),breaks = seq(-200,1300,500))+ 
  stat_compare_means(aes(group=group),label = "p.format", hide.ns = T, size = 6,
                     #bracket.size = 20,label.x = 1.35,label.y = 1200,size = 9
                     )+
  theme_gray()+
  theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5))+
  labs(title = "Tumor Mutation Burden")+ #Mutation Load Plot
  theme(plot.title = element_text(size = 20, color = "black", 
                                  hjust = 0.5, 
                                  angle = 0)) 
p+  theme(legend.text = element_text(size = 14, face = "bold"), 
          legend.title = element_text(size = 14,face = "bold")) +
    labs(x = NULL, y = NULL)

library(export)
graph2ppt(file="output/plots/TMB.pptx", width = 8, aspectr = 1.5)
############ 更好看的箱线图
# TCGAtmbbox$value <- as.numeric(TCGAtmbbox$value)
# compare_means(value ~ group, data = TCGAtmbbox)
# ggboxplot(TCGAtmbbox,x="variable",y="value",fill = "group",size = 0.1,
#           palette = c("#00CCFF","#FF3333"))+
#   xlab("")+ylab("Log(Copy number)")+
#   theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 0.5, vjust = 0.5))+
#   scale_y_continuous(limits = c(-50, 200),breaks = seq(-50,200,50))+
#   stat_compare_means(aes(group=group),label = "p.signif",hide.ns = T,
#                      bracket.size = 10)+
#   theme_gray()

ggviolin(x="variable",y="value",fill = "group",data = TCGAtmbbox,size = 0.1,
         palette = c("#00CCFF","#FF3333"),add = "boxplot",position = position_dodge(1), 
         add.params = list(group = "group"))+
  xlab("")+
  ylab("")+
  #scale_y_continuous(limits =c(-200, 1300),breaks = seq(-200,1300,500))+ 
  stat_compare_means(aes(group=group),label = "p.format", hide.ns = T,size = 9)+
                     # bracket.size = 20,label.x = 1.35,label.y = 1200,
  theme_cleveland()+
  theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5))+
  labs(title = "Tumor Mutation Burden")+ #Mutation Load Plot
  theme(plot.title = element_text(size = 20, color = "black", 
                                  hjust = 0.5, 
                                  angle = 0)) +
  theme(legend.text = element_text(size = 14, 
                                   face = "bold"), 
         legend.title = element_text(size = 14,face = "bold")) +
   labs(x = NULL, y = NULL)
graph2ppt(file="output/plots/TMB.pptx", append=T, width = 8, aspectr = 1.5)

