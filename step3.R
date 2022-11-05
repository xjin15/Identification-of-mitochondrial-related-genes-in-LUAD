p_load(tidyverse)

mtgene1626 <- data.table::fread(file = "data/human mitochondria genes.tsv") %>% 
  rename(esid = V1,
         Symbol = V3,
         cat  = V4,
         chr  = V5,
        )
mtgene1157 <- data.table::fread(file = "data/mitocarta.csv")

fpkm <-read.csv(file = "data/TCGA-LUAD.htseq_fpkm.tsv", 
                 header = T, sep = "\t", stringsAsFactors = F, 
                 row.names = 1, check.names = F)#第一列做行名
save(fpkm,mtgene1626,mtgene1157,file = "outdata/fpkm.rds")
load(file = "outdata/fpkm.rds")
fpkm[1:4,1:4]
str(fpkm)
gtf_gene <- read.csv(file = "../maf分组/output/gtfmRNA22.txt", header = T,sep = "\t")
gtf_gene$gene_type %>% unique()
MTgene <- gtf_gene$gene_name %>% grep(pattern = "MT-",value = T)
MTgene <- MTgene[-1]
ncol(fpkm)
fpkm1 <- fpkm[,substring(text = colnames(fpkm), 
             first = 14, last = 16) == "01A"]
fpkm_NT <- fpkm %>% 
  select(ends_with("11A"))
fpkm_TP <- fpkm %>% 
  select(ends_with("01A"))
fpkm_0111 <- fpkm %>% 
  select(ends_with("01A") | ends_with("11A"))
dim(fpkm_0111)
# 转ID
rt <- fpkm_0111 %>% 
  rownames_to_column("rn")
rt <- rt %>% 
  mutate(rn = substr(rn,start = 1,stop = 15))
gn <- gtf_gene$gene_name[match(rt$rn,table = gtf_gene$gene_id,nomatch = 0)]
rt$rn <- gn
rt <- rt %>%select(-1)
rt <- as.matrix(rt)
rownames(rt) <- gn
rt[1:4,1:4]
rowMeans(rt[MTgene,]) %>% format(.,scientific=F) 
rowSums(rt[MTgene,]) %>% format(.,scientific=F)
rt <- rt[rowSums(rt)>0.05,]
dim(rt)
rt1 <- rt %>% as.data.frame() %>% 
  filter(rownames(rt) %in% MTgene)


rt1_NT <- rt1[,511:568]
rt1_TP <- rt1[,1:510]
rt2 <- t(rt1) %>% as.data.frame()
rt2$group <- "NT"
rt2$group[1:510] <- "TP"
rt_plot <- rt2 %>% pivot_longer(cols = starts_with("MT"),
                                names_to = "MTgene",
                                values_to = "value")

library(ggpubr)
ggboxplot(rt_plot,x="MTgene",y="value",fill = "group",size = 0.1,)
ggboxplot(rt_plot,x="MTgene",y="value",fill = "group",size = 0.1,
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
x <- compare_means(value ~ group, group.by = "MTgene",data = rt_plot)
xnosig <- x[which(x$p > 0.05),]
rt_plotsig <- rt_plot %>% 
  filter(!MTgene %in% xnosig$MTgene)
ggboxplot(rt_plotsig,x="group",y="value",fill = "group",size = 0.1,
          facet.by = "MTgene",
          palette = c("#00CCFF","#FF3333"))+
  xlab("")+ylab("InfiltrationScore")+
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 0, vjust = 0))+
  # scale_y_continuous(limits = c(-1, 1),breaks = seq(-1,1,0.2))+
  stat_compare_means(aes(group=group),label = "p.format",hide.ns = T,
                     bracket.size = 10)+
  theme_gray()+
  theme(axis.text = element_text(size = 8), 
        axis.text.x = element_text(size = 8, angle = 45, vjust = 0.6))


ggboxplot(rt_plotsig,x="MTgene",y="value",fill = "group",size = 0.1,
          palette = c("#00CCFF","#FF3333"))+
  xlab("")+ylab("InfiltrationScore")+
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 0, vjust = 0))+
   scale_y_continuous(limits = c(-5, 15),breaks = seq(-5,15,5))+
  stat_compare_means(aes(group=group),label = "p.format",hide.ns = T,
                     bracket.size = 10)+
  theme_gray()+
  theme(axis.text = element_text(size = 8), 
        axis.text.x = element_text(size = 8, angle = 45, vjust = 0.6))

ggviolin(rt_plotsig,x="group",y="value",fill = "group",size = 0.1,
         facet.by = "MTgene",
          palette = c("#00CCFF","#FF3333"))+
  xlab("")+ylab("InfiltrationScore")+
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 0, vjust = 0))+
  scale_y_continuous(limits = c(-5, 15),breaks = seq(-5,15,5))+
  stat_compare_means(aes(group=group),label = "p.format",hide.ns = T,
                     bracket.size = 10)+
  theme_gray()+
  theme(axis.text = element_text(size = 8), 
        axis.text.x = element_text(size = 8, angle = 45, vjust = 0.6))

rt_plotsig$MTgene %>% unique() %>% dput()


