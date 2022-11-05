### 加载数据和R包########
rm(list = ls())
load(file = "outdata/step4_group_data.rds")
load(file = "outdata/step3_julei.rds")
fpkm <- fpkm[,allid]
sample_group$group |> table() ## 第一组279，第二种类218
p_load(limma,edgeR,ggplot2,ggthemes,ggpubr,ggrepel,export,data.table,tidyverse)
# 设定阈值
rt <- as.matrix(fpkm)   
rt[1:3,1:3]
rt <- rt[rowSums(rt) > 0, ]    # 57492个基因  
dim(rt) # 57492,497样本
####### 基因ID转换及选择mrna #####
{
  
  gtf_gene <- read.csv(file = "../maf分组/output/gtfmRNA22.txt", header = T,sep = "\t")
  ## 选择proteincoding的蛋白质
  gtf_gene_protein <- gtf_gene %>% 
    dplyr::rename( Ensembl_ID = gene_id,
                   Symbol  = gene_name,
                   Biotype = gene_type  ) %>% 
    dplyr::filter(Biotype %in% c("IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_LV_gene", 
                                 "IG_M_gene", "IG_V_gene", "IG_Z_gene", 
                                 "nonsense_mediated_decay", "nontranslating_CDS", 
                                 "non_stop_decay", "polymorphic_pseudogene", 
                                 "protein_coding", "TR_C_gene", "TR_D_gene", "TR_gene", 
                                 "TR_J_gene", "TR_V_gene"))
  # 基因ID去掉点号后面的版本标识，准备转换ID。
  rownames(rt) <- substring(text = rownames(rt), first = 1, last = 15)
  # 确认没重复的基因ENSEMBL ID
  stopifnot(
    {
      duplicated(rownames(rt))== F
    }
  )
  # 查看行和列数
  dim(rt) #57492 497(mut111,wild382)
  matchgene <- intersect(rownames(rt), gtf_gene_protein$Ensembl_ID)
  mrna <- rt[matchgene, ]### 20039
  # 可以转换ID
  rownames(mrna) <- gtf_gene_protein$Symbol[match(x = rownames(mrna),table = gtf_gene_protein$Ensembl_ID)]
  # 重复的基因名字进行合并且取均值
  dim(mrna)
  mrna <- avereps(mrna) #19971 
}
allid497 <- allid

save(mrna,allid497,file = "outdata/step4.2_mrna.rds")
