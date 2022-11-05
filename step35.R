# 需要先开始做9.4批量生存分析和9.1 lASSO回归分析
# -------------------------------------------------------------------------


aaa <- p_COXsig_deg1747 %>% as.data.frame() 
aaa <- aaa %>% filter(p < 0.05)
rownames(aaa) %>% grep(pattern ="DA" ,value = T)
intersect(rownames(aaa),choose_gene_1se)
intersect(rownames(aaa),choose_gene_1se2)
intersect(rownames(aaa),choose_gene_min)
intersect(rownames(aaa),choose_gene_min2)
intersect(rownames(aaa),choose_gene_min2)

length(choose_gene_1se)
length(choose_gene_1se2)
length(choose_gene_min)
length(choose_gene_min2)


x=list(DEmRNAs =deg_mito,
       OS__genes = rownames(aaa),
       LASSO_reducted_Genes = c(choose_gene_min2,'DARS2','COX5B'))
draw_venn(x,name = 'The Venn Plot of hub genes')
graph2pdf('output/plots/Vennplot.pdf')
