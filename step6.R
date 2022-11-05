rm(list = ls())
p_load(tidyverse,ggplot2,ggpubr)

load(file = "outdata/LUAD_TCGA_ALLfpkm.rds")
fpkm[1:4,1:4]
str(fpkm)
mtgene1626$Symbol2 <- gtf_gene$gene_name[base::match(mtgene1626$esid,gtf_gene$gene_id,nomatch = 0)]
identical(mtgene1626$Symbol,mtgene1626$Symbol2)
setequal(mtgene1626$Symbol,mtgene1626$Symbol2)
setdiff(mtgene1626$Symbol,mtgene1626$Symbol2)

test1 <- mtgene1626 %>%
  select(starts_with("Symbol"))

test1 <- test1 %>% 
  mutate(de = ifelse(Symbol == Symbol2,0,1))
####上面的步骤我是看了一下GTF和mitocarta用的注释版本，明显不一样，但是不影响
# 我用的是esid来取

# fpkm提取线粒体基因1626用来聚类 -----------------------------------------------------

rownames(fpkm) <- substring(rownames(fpkm),1,15)
fpkm[1:4,1:4]
rt <- fpkm[mtgene1626$esid,]
rownames(rt) <- mtgene1626$Symbol
rt_NT <- rt %>% 
  select(ends_with("11A"))
rt_TP <- rt %>% 
  select(ends_with("01A"))
rt2 <- rt %>% 
  select(ends_with("01A") | ends_with("11A"))
rt2 <- rt2[rowMeans(rt) > 0.5,]

# 准备作图数据

rt3 <- t(rt2) %>% as.data.frame()
rt3$group <- "NT"
rt3$group[1:510] <- "TP"
rt_plot <- rt3 %>% pivot_longer(cols = -1508,
                                names_to = "MTgene",
                                values_to = "value")

# 选择比较有意义的来聚类 -------------------------------------------------------------
x <- compare_means(value ~ group, group.by = "MTgene",data = rt_plot)
xnosig <- x[which(x$p > 0.05),]
xsig <- x[which(x$p.adj < 0.05), ]
rt_plotsig <- rt_plot %>% 
  filter(!MTgene %in% xnosig$MTgene)
ggboxplot(rt_plotsig,x="MTgene",y="value",fill = "group")
str(xsig)
xsig1 <- xsig[,c(1,5,6,7,8,9)]
xsig1$func_cat <- mtgene1626$cat[match(xsig1$MTgene,mtgene1626$Symbol,no = 0)]
write.csv(xsig1,file = "outdata/mtgene_NTvsTP.csv")
siggene <- setdiff(rownames(rt2),xnosig$MTgene)

save(siggene,rt_TP,fpkm,file = "outdata/step3_julei.rds")
