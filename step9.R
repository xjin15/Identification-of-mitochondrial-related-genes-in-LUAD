rm(list = ls())
load(file = "outdata/step4.3_enrich_Analysis.rds")
library(pacman)
p_load(tidyverse,org.Hs.eg.db, clusterProfiler, DOSE, enrichplot)
p_load(enrichplot, ggupset, ggnewscale)


# 1.准备数据 ------------------------------------------------------------------
# 从差异基因中选择前100上调的基因
all_diff$group |> table()
###画更好看的GO和KEGG图
allDiff <- all_diff |> 
  rename(log2FoldChange = logFC,
         padj = adj.P.Val,
         g = group,
         gene_id = ID)

rownames(allDiff) <- NULL
## 定义差异基因标准
deg <- allDiff %>% 
  filter(abs(log2FoldChange) > 1) %>% 
  filter(padj < 0.05)

deg <- deg %>% 
  column_to_rownames("gene_id")

## 不同的阈值，筛选到的差异基因数量就不一样，后面的超几何分布检验结果就大相径庭。

logFC_t= 0
deg$g=ifelse(deg$padj > 0.05,'stable',
             ifelse( deg$log2FoldChange > logFC_t,'UP',
                     ifelse( deg$log2FoldChange < -logFC_t,'DOWN','stable') ))
table(deg$g)
head(deg)
deg$symbol=rownames(deg)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)

df %>% anyDuplicated()
head(df)
DEG=deg
head(DEG)

DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
head(DEG)
#save(DEG,file = 'anno_DEG.Rdata')

gene_up= DEG[DEG$g == 'UP','ENTREZID'] 
gene_down=DEG[DEG$g == 'DOWN','ENTREZID'] 
gene_diff=c(gene_up,gene_down)
gene_all=as.character(DEG[ ,'ENTREZID'] )
data(geneList, package="DOSE")
head(geneList)
boxplot(geneList)
boxplot(DEG$log2FoldChange)


go.up <- enrichGO(gene = gene_up,
                  OrgDb = "org.Hs.eg.db", 
                  ont="all",
                  universe = gene_all,
                  pvalueCutoff = 0.05, ##根据具体情况调整
                  qvalueCutoff =0.05)  ####根据具体情况调整
head(go.up)[,1:6]
dotplot(go.up ) + ggtitle("upregulate gene go analysis")
graph2ppt(file = "output/plots/Enrich_Analysis.pptx", width= 10, append = T, aspectr = 1.5)
#ggsave('go.up.dotplot.pdf')
go.down <- enrichGO(gene = gene_down,
                    OrgDb = "org.Hs.eg.db", 
                    ont="all",
                    universe = gene_all,
                    pvalueCutoff = 0.05, ##根据具体情况调整
                    qvalueCutoff =0.05)  ##根据具体情况调整
head(go.down)[,1:6]
dotplot(go.down ) + ggtitle("downregulate gene go analysis")
graph2ppt(file = "output/plots/Enrich_Analysis.pptx", width= 10, append = T, aspectr = 1.5)

#ggsave('go.down.dotplot.pdf')
go.diff <- enrichGO(gene = gene_diff,
                    OrgDb = "org.Hs.eg.db", 
                    ont="all",
                    universe = gene_all,
                    pvalueCutoff = 1, ##根据具体情况调整
                    qvalueCutoff = 1)  ##根据具体情况调整
head(go.diff)[,1:6]
dotplot(go.diff )
graph2ppt(file = "output/plots/Enrich_Analysis.pptx", width= 10, append = T, aspectr = 1.5)

#ggsave('go.diff.dotplot.pdf')

go_diff_dt <- as.data.frame(go.diff)
go_down_dt <- as.data.frame(go.down)
go_up_dt <- as.data.frame(go.up)

## 上面用的是富集分析中的pvalue,条目太多了，这里p.adjust
## 注意函数New.functions.R也要将pvalue改为p.adjust
down_go<-go_down_dt[go_down_dt$p.adjust<0.00001,]##根据具体情况调整
down_go$group=-1
up_go<-go_up_dt[go_up_dt$p.adjust<0.00001,]##根据具体情况调整
up_go$group=1

source('scripts/news.functions.R')
## go.kegg_plot()一步画图，只需要输入画图数据up.data，down.data
go_figure=go.kegg_plot(up.data=up_go,down.data=down_go)
print(go_figure)
graph2ppt(file = "output/plots/Enrich_Analysis.pptx", width= 12, append = T, aspectr = 1.5)



# 2.批量的Go和KEGG的富集分析 -------------------------------------------------------------
enrich_result <- list()
# 批量富集分析，上调基因和下调基因的

for (i in c("gene_up","gene_down","gene_diff")) {
 
  enrich_result[paste0("enrichGO_",i)] <- enrichGO(gene = base::eval(parse(text = i)),
                                                    OrgDb = "org.Hs.eg.db", 
                                                    ont="all",
                                                    universe = gene_all,
                                                    pvalueCutoff = 0.05, ##根据具体情况调整
                                                    qvalueCutoff =0.05)  ####根据具体情况调整
   print(paste0("enrichGO_",i,"was done"))
   enrich_result[paste0("enrichKEGG_",i)] <- enrichKEGG(
                                                        gene  =  base::eval(parse(text = i)),
                                                        organism     = 'hsa',
                                                        #universe     = gene_all,
                                                        # pvalueCutoff = 0.05,
                                                        # qvalueCutoff =0.05,
                                                        #use_internal_data = F #使用在线最新数据
                                                      )
   print(paste0("enrichKEGG_",i,"was done"))
}
######### KEGG 的上下调画到一幅图中#####

enrich_kegg_diff_dt <- as.data.frame(enrich_result$enrichKEGG_gene_diff)
enrich_kegg_up_dt <- as.data.frame(enrich_result$enrichKEGG_gene_up)
enrich_kegg_down_dt <- as.data.frame(enrich_result$enrichKEGG_gene_down)

## 上面用的是富集分析中的pvalue,条目太多了，这里p.adjust
## 注意函数New.functions.R也要将pvalue改为p.adjust
down_enrich_kegg<-enrich_kegg_down_dt[enrich_kegg_down_dt$p.adjust < 0.01,]##根据具体情况调整
down_enrich_kegg$group=-1
up_enrich_kegg<-enrich_kegg_up_dt[enrich_kegg_up_dt$p.adjust<0.01,]##根据具体情况调整
up_enrich_kegg$group=1

source('scripts/news.functions.R')
## go.kegg_plot()一步画图，只需要输入画图数据up.data，down.data
go_figure=go.kegg_plot(up.data=up_enrich_kegg,down.data=down_enrich_kegg)
print(go_figure)
graph2ppt(file = "output/plots/上下调的差异基因富集KEGG不同.pptx", width= 12, append = T, aspectr = 1.5)
######GO的上下调富集分析画到一幅图中####
enrich_go_diff_dt <- as.data.frame(enrich_result$enrichGO_gene_diff)
enrich_go_up_dt <- as.data.frame(enrich_result$enrichGO_gene_up)
enrich_go_down_dt <- as.data.frame(enrich_result$enrichGO_gene_down)

## 上面用的是富集分析中的pvalue,条目太多了，这里p.adjust
## 注意函数New.functions.R也要将pvalue改为p.adjust
down_enrich_go<-enrich_go_down_dt[enrich_go_down_dt$p.adjust < 0.000001
                                  | str_detect(enrich_go_down_dt$Description,pattern = "immune") ,]##根据具体情况调整
down_enrich_go$group=-1
up_enrich_go<-enrich_go_up_dt[enrich_go_up_dt$p.adjust<0.000001,]##根据具体情况调整
up_enrich_go$group=1

source('scripts/news.functions.R')
## go.kegg_plot()一步画图，只需要输入画图数据up.data，down.data
go_figure=go.kegg_plot(up.data=up_enrich_go,down.data=down_enrich_go)
print(go_figure)
graph2ppt(file = "output/plots/上下调的差异基因富集KEGG不同.pptx", width= 12, append = T, aspectr = 1.5)

dotplot(enrich_result$enrichGO_gene_diff)
dotplot(enrich_result$enrichKEGG_gene_diff)
dotplot(enrich_result$enrichGO_gene_up)
dotplot(enrich_result$enrichGO_gene_down)
dotplot(enrich_result$enrichKEGG_gene_up)
dotplot(enrich_result$enrichKEGG_gene_down)

# 3.GSEA富集分析数据准备  ------------------------------------------------------------------
geneList <- all_diff[order(all_diff$logFC,decreasing = T),"logFC",drop = F] |> 
  as.matrix() %>% .[,1]
geneList
geneList <- sort(geneList, decreasing = T)
anyDuplicated(names(geneList))

# 转换ID,会有少部分没有
ids <- bitr(
  geneID = names(geneList),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
) 

anyDuplicated(ids$SYMBOL) # 转换ID 出现了重复
ids[duplicated(ids$SYMBOL),] #重复的ID 名字
which(ids$SYMBOL %in% ids$SYMBOL[duplicated(ids$SYMBOL)] )# 显示重复ID位置
ids[ids$SYMBOL %in% ids$SYMBOL[duplicated(ids$SYMBOL)] , ]# 显示重复的symbolid


uni_ids <- ids[!duplicated(ids$SYMBOL),] # 对重复值只保留第一个出现的
uni_ids1 <- ids %>% distinct(SYMBOL,.keep_all = T) # distint函数
## distinct函数相当于对表格的unique函数，
# 如果不加keep参数则只提取出来这一列
# 加上keep参数就把其他的列也提取出来。

stopifnot(setequal(uni_ids,uni_ids1)) 
# 显示distinct函数和duplicated函数都是保留第一个出现的重复值

geneList <- geneList[names(geneList) %in% uni_ids$SYMBOL]
names(geneList) <- uni_ids$ENTREZID

geneList <- sort(geneList,decreasing = T)
# 最后确认一次geneList排好序。
head(geneList)
####注意：这里是ClusterB/ClusterA 所以NES＞1的属于在B中富集


# 4.Hallmark的GSEA ---------------------------------------------------------
# hallmark数据集准备
HALLMARK_df = msigdbr(species = "Homo sapiens",
                      category = "H") 
HALLMARK_df <- HALLMARK_df %>% 
dplyr::select(gs_name,gene_symbol,)
HALLMARK_df$gs_name <- substring(HALLMARK_df$gs_name,first = 10)

geneList_s <- all_diff[order(all_diff$logFC,decreasing = T),"logFC",drop = F] |> 
  as.matrix() %>% .[,1]
geneList_s <- sort(geneList_s, decreasing = T)
geneList_s
#### 还没弄好代码
gse_h <- GSEA(
  geneList     = geneList_s,      #基因list
  TERM2GENE = HALLMARK_df,
  # TERM2NAME = ,
  # nPerm      = 1000,       现在已经不推荐使用了
  minGSSize    = 10,            #基因集最小
  maxGSSize    = 500,           #基因集最大
  pvalueCutoff = 0.05,          # p值阈值
  exponent     = 1,             # step权重默认即可
  by           = "fgsea",       # 使用的方法
  verbose      = T  ,            # 过程可视化，随意选T
  eps          = 0
)

gseaplot2(gse_h, geneSetID = 1, title = gse_h$Description[1])


down_gsea_hallmark<-gse_h@result[gse_h@result$p.adjust < 0.01
                                 & gse_h@result$NES < 0,]##根据具体情况调整
down_gsea_hallmark$group=-1
up_gsea_hallmark<-gse_h@result[gse_h@result$p.adjust < 0.01
                                 & gse_h@result$NES > 0,]##根据具体情况调整
up_gsea_hallmark$group=1

source('scripts/news.functions.R')
## go.kegg_plot()一步画图，只需要输入画图数据up.data，down.data
gse_hallmark_figure=go.kegg_plot(up.data=up_gsea_hallmark,down.data=down_gsea_hallmark)
print(gse_hallmark_figure)
export::graph2ppt(file = "output/plots/GSEA_Analysis_new.pptx",width = 15,aspectr =1.2)







