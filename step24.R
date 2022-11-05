rm(list = ls())
# 导入
load(file = "outdata/LUAD_TCGA_ALLfpkm.rds")
load(file = "outdata/step4.3_enrich_Analysis.rds")

p_load(tidyverse)

fpkm1 <- fpkm[,substring(text = colnames(fpkm), 
                         first = 14, last = 16) == "01A"]
fpkm[1:4,1:4]
fpkm_NT <- fpkm |>  dplyr::select(ends_with("11A"))
fpkm_TP <- fpkm |>   dplyr::select(ends_with("01A"))
fpkm_0111 <- fpkm %>% 
  dplyr::select(ends_with("01A") | ends_with("11A"))
dim(fpkm_0111)


# 转ID
rt <- fpkm_0111 %>% 
  rownames_to_column("rn")
rt <- rt %>% 
  mutate(rn = substr(rn,start = 1,stop = 15))
gn <- gtf_gene$gene_name[match(rt$rn,table = gtf_gene$gene_id,nomatch = 0)]
rt$rn <- gn
rt <- rt %>% dplyr::select(-1)
rt <- as.matrix(rt)
rownames(rt) <- gn
rt[1:4,1:4]

deg_mito <- read.table(file = "outdata/deg_mito.tsv") |> as.vector() |> unlist()

DEGgene <- rownames(rt)

# 求行平均值


rowMeans(rt[DEGgene,]) %>% format(.,scientific=F) 
rowSums(rt[DEGgene,]) %>% format(.,scientific=F)
rt <- rt[rowSums(rt) > 5,]
dim(rt)  #57364 568
rt1 <- rt %>% as.data.frame() %>% 
  filter(rownames(rt) %in% DEGgene)


rt1_NT <- rt1[,511:568]
rt1_TP <- rt1[,1:510]
rt2 <- t(rt1) %>% as.data.frame()
rt2$group <- "NT"
rt2$group[1:510] <- "TP"
rt2 <- rt2 |> dplyr::select(group,everything()) |> 
  mutate(group = as.factor(group))
rt_plot <- rt2 %>% pivot_longer(cols = 2:ncol(rt2),
                                names_to = "DEGgene",
                                values_to = "value")

library(ggpubr)
ggboxplot(rt_plot,x="DEGgene",y="value",fill = "group",size = 0.1,)
ggboxplot(rt_plot,x="DEGgene",y="value",fill = "group",size = 0.1,
          palette = c("#00CCFF","#FF3333"))+
  xlab("")+ylab("")+
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 0, vjust = 0))+
  #scale_y_continuous(limits = c(0, 10),breaks = seq(0,10,2))+
  stat_compare_means(aes(group=group),label = "p.signif",hide.ns = T,
                     bracket.size = 10)+
  theme_gray()+
  theme(axis.text = element_text(size = 8), 
        axis.text.x = element_text(size = 8, angle = 45, vjust = 0.6)) 
library(export)
####只选择有意义的
x <- compare_means(value ~ group, group.by = "DEGgene",data = rt_plot)
write.table(x,file = "outdata/所有基因在LUAD肿瘤和癌旁的差异.tsv",
            sep = "\t",row.names = F,col.names = T)


{
# #####后面的代码可以不用
# xnosig <- x[which(x$p > 0.05),]
# rt_plotsig <- rt_plot %>% 
#   filter(!DEGgene %in% xnosig$MTgene)
# ggboxplot(rt_plotsig,x="group",y="value",fill = "group",size = 0.1,
#           facet.by = "DEGgene",
#           palette = c("#00CCFF","#FF3333"))+
#   xlab("")+ylab("InfiltrationScore")+
#   theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 0, vjust = 0))+
#   # scale_y_continuous(limits = c(-1, 1),breaks = seq(-1,1,0.2))+
#   stat_compare_means(aes(group=group),label = "p.format",hide.ns = T,
#                      bracket.size = 10)+
#   theme_gray()+
#   theme(axis.text = element_text(size = 8), 
#         axis.text.x = element_text(size = 8, angle = 45, vjust = 0.6))
# 
# 
# ggboxplot(rt_plotsig,x="MTgene",y="value",fill = "group",size = 0.1,
#           palette = c("#00CCFF","#FF3333"))+
#   xlab("")+ylab("InfiltrationScore")+
#   theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 0, vjust = 0))+
#   scale_y_continuous(limits = c(-5, 15),breaks = seq(-5,15,5))+
#   stat_compare_means(aes(group=group),label = "p.format",hide.ns = T,
#                      bracket.size = 10)+
#   theme_gray()+
#   theme(axis.text = element_text(size = 8), 
#         axis.text.x = element_text(size = 8, angle = 45, vjust = 0.6))
# 
# ggviolin(rt_plotsig,x="group",y="value",fill = "group",size = 0.1,
#          facet.by = "MTgene",
#          palette = c("#00CCFF","#FF3333"))+
#   xlab("")+ylab("InfiltrationScore")+
#   theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 0, vjust = 0))+
#   scale_y_continuous(limits = c(-5, 15),breaks = seq(-5,15,5))+
#   stat_compare_means(aes(group=group),label = "p.format",hide.ns = T,
#                      bracket.size = 10)+
#   theme_gray()+
#   theme(axis.text = element_text(size = 8), 
#         axis.text.x = element_text(size = 8, angle = 45, vjust = 0.6))
# 
# rt_plotsig$MTgene %>% unique() %>% dput()
}


