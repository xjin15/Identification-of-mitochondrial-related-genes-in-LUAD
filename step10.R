rm(list = ls())
# 缩短名称函数 有些名称实在太长了
# shorten_names <- function(x, n_word=4, n_char=30){
#   if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 30))
#   {
#     if (nchar(x) > 30) x <- substr(x, 1, 30)
#     x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
#                      collapse=" "), "...", sep="")
#     return(x)
#   }
#   else
#   {
#     return(x)
#   }
# }
# 数据准备 -------------------------------------------------------------

# 输入筛到的基因，一般选择logFC和adj.p符合条件的基因做GO和KEGG，
# 也有的是选上调和下调分别前100个基因
# 为了完善自己的代码，我选择做上调100，下调100和全部200个基因做GO 和 KEGG，
# 然后是用所有的去做GSEA分析
 
load(file = "outdata/step4.3_enrich_Analysis.rds")
pacman::p_load(tidyverse,org.Hs.eg.db, clusterProfiler, DOSE, enrichplot)
pacman::p_load(enrichplot, ggupset, export, ggnewscale)
geneList <- all_diff[order(all_diff$logFC,decreasing = T),"logFC",drop = F] |> 
  as.matrix() %>% .[,1]
geneList <- sort(geneList, decreasing = T)
geneList
anyDuplicated(names(geneList))
# 转换ID，对genelist进行
{
ids <- bitr(
  geneID = names(geneList),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
) 

anyDuplicated(ids$SYMBOL) # 转换ID 出现了重复
ids[duplicated(ids$SYMBOL),] #重复的ID 名字
which(ids$SYMBOL %in% ids$SYMBOL[duplicated(ids$SYMBOL)] )
ids[ids$SYMBOL %in% ids$SYMBOL[duplicated(ids$SYMBOL)] , ]

uni_ids <- ids[!duplicated(ids$SYMBOL),]
uni_ids1 <- ids %>% distinct(SYMBOL,.keep_all = T)
stopifnot(setequal(uni_ids,uni_ids1))

geneList <- geneList[names(geneList) %in% uni_ids$SYMBOL]
names(geneList) <- uni_ids$ENTREZID
geneList <- sort(geneList,decreasing = T)
}

geneList # 转换成了名字是entrezid的Named num向量

## 下调基因 和 上调基因，这是按照 logFC 挑选出来的

mydf <- data.frame(Entrez=names(geneList), FC=geneList)
mydf <- mydf[abs(mydf$FC) > 1,]
mydf$group <- "upregulated"
mydf$group[mydf$FC < 0] <- "downregulated"


#用compareCluster -------------------------------------------------------------------------
# 利用Y叔的compareCluster包 比较上下调基因的GO 富集分析
# compareGO分析
compare_ego <- compareCluster(Entrez~group, 
                              data=mydf, 
                              fun="enrichGO",
                              OrgDb='org.Hs.eg.db')
head(compare_ego)
dotplot(compare_ego) #点图，会显示两个

compare_ego <- setReadable(compare_ego,OrgDb='org.Hs.eg.db') # id转成symbol
cnetplot(compare_ego)
emapplot(pairwise_termsim(simplify(compare_ego)))  # emaplot


### compare kegg分析
compare_eke <- compareCluster(Entrez~group, 
                              data=mydf, 
                              fun="enrichKEGG",OrgDb='org.Hs.eg.db'
                              )

head(compare_eke)
dotplot(compare_eke) #点图，会显示两个
compare_eke <- setReadable(compare_eke,OrgDb='org.Hs.eg.db',keyType = "ENTREZID") # id转成symbol
cnetplot(compare_eke)



## comparedo 分析
compare_edo <- compareCluster(Entrez~group, 
                              data=mydf, 
                              fun="enrichDO"#,OrgDb='org.Hs.eg.db'
                              )

head(compare_edo)
dotplot(compare_edo) #点图，会显示两个


#  compare pathway 分析   enrichPathway 有问题做不了分析
# compare_epa <- compareCluster(Entrez~group, 
#                               data=mydf, 
#                               fun = "enrichPathway",
#                               OrgDb='org.Hs.eg.db'
# )
# 
# head(compare_epa)
# dotplot(compare_epa) #点图，会显示两个


## KEGG显示

gene_down <- DEG[DEG$regulation == "up",c("genes","fold_change")]
gene_up <- DEG[DEG$regulation == "down",c("genes","fold_change")]

colnames(gene_down) <- c("symbol", "logFC")
colnames(gene_up) <- c("symbol", "logFC")
gene_sig <- DEG[,c("genes","fold_change")]
colnames(gene_sig) <- c("symbol", "logFC")
# 所有基因作为背景基因集
backgene <- bitr(geneID = names(geneList), fromType = "SYMBOL", 
                 toType = "ENTREZID",drop = T,
                 OrgDb = org.Hs.eg.db)
names(geneList) <- geneList_name$ENTREZID
names(geneList) %>% duplicated() %>% table()

CaseGene <- gene_sig$symbol

# 开始GO和KEGG富集分析及其可视化 --------------------------------------------------------
  CaseGeneSet <- bitr(geneID = CaseGene,
               fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = "org.Hs.eg.db"
                )
  Case_KEGG <- enrichKEGG(
                         gene  =  CaseGeneSet$ENTREZID,
                         organism     = 'hsa',
                         #universe     = gene_all,
                         # pvalueCutoff = 0.05,
                         # qvalueCutoff =0.05,
                         #use_internal_data = F #使用在线最新数据
                         )
  Case_KEGG@result$Description <- shorten_names(Case_KEGG@result$Description)
## GO分析
  Case_go <- enrichGO(gene = CaseGeneSet$ENTREZID, 
                 OrgDb = "org.Hs.eg.db", 
                 # pool = T, 分析子类型时是否使用所有的背景基因。T 则能富集到更多的途径
                 ont="all", #GO分析三个都要MF分子功能，CC细胞成分，BP生物学过程
                 # pvalueCutoff = 0.05, 默认0.05
                 # qvalueCutoff = 0.2, ，默认 0.2
                 #universe = names(geneList$ENTREZID), 我试过加不加universe没有区别
                 readable = T
                 ) 
  head(Case_go)
  Case_go@result$Description <- shorten_names(Case_go@result$Description)
##### kegg 分析可视化#######
  # 可视化主要有几个图：
  # 1.最经典的dotplot和barplot
  # 2.emaplot：enrichment map。将terms组织成网络，
      # 有重复基因的terms就会连接在一起并且靠近。这样的好处是可以识别功能terms的聚集
  # 3. cnetpplot category netplot
      # 把基因和富集到的通路联系，描述成一个网络，有助于了解哪些基因参与了多个富集的通路。
  # 4.山脊图。ridge plot
      # 感觉类似于barplot之类的，就是按基因集分组，
      # 密度图由每组内每个基因的fold change值的频率生成。
      # 有助于解释上调/下调的通路。
  Case_KEGG <- DOSE::setReadable(Case_KEGG, OrgDb='org.Hs.eg.db',keyType='ENTREZID')#加上symbol
  upsetplot(Case_KEGG, n = 10,title=" upsetplot of KEGG Enrichment ")
  dotplot(Case_KEGG,title=" Dotplot of KEGG Enrichment ")#绘制气泡图
  graph2ppt(file = "output/plots/Enrich_Analysis.pptx", width= 10, append = T, aspectr = 1.5)
  barplot(Case_KEGG, showCategory=15,title="Barplot of KEGG Enrichment")#绘制条形图
 
  
  cnetplot(Case_KEGG, categorySize="pvalue", colorEdge = TRUE) + ggtitle("cnetplot of KEGG analysis")
  graph2ppt(file = "output/plots/Enrich_Analysis.pptx", width= 10, append = T, aspectr = 1.5)
  
  # cnetplot(Case_KEGG,  circular = TRUE, colorEdge = TRUE) 
  #网络图，每个节点是一个富集到的pathway，
  # 若pathways之间有重叠的感兴趣基因，则自动将这两个通路用线连接。
  # 默认画top30个富集到的pathways, 节点大小对应该pathway下富集到的感兴趣基因个数，
  # 节点的颜色对应p.adjust的值，从小到大，对应蓝色到红色。
  Case_KEGG_emaplot <- pairwise_termsim(Case_KEGG)
  emapplot(Case_KEGG_emaplot) + ggtitle("Emaplot of KEGG analysis")
  # graph2ppt(file = "output/plots/Enrich_Analysis.pptx", width= 10, append = T, aspectr = 1.5)
  heatplot(Case_KEGG) + ggtitle("Heatplot of KEGG analysis")
  
  
  
 ##### GO分析可视化#######
  heatplot(Case_go) + ggtitle("Heatplot of GO analysis")
  upsetplot(Case_go, n = 10,title=" upsetplot of GO Enrichment ")
  dotplot(Case_go, title=" Dotplot of GO Enrichment ")#绘制气泡图
  barplot(Case_go, showCategory=15,title="Barplot of GO Enrichment")#绘制条形图
  dotplot(Case_go, split="ONTOLOGY")+ facet_grid(ONTOLOGY ~ ., scale="free") + 
    ggtitle("Dotplot of GO analysis")
  graph2ppt(file = "output/plots/Enrich_Analysis.pptx", width= 10, append = T, aspectr = 1.5)
  # cnetplot(Case_go,  circular = TRUE, colorEdge = TRUE) 
  cnetplot(Case_go, categorySize="pvalue",  colorEdge = TRUE)+ ggtitle("ctnet of GO analysis")
  graph2ppt(file = "output/plots/Enrich_Analysis.pptx", width= 10, append = T, aspectr = 1.5)
  Case_go_emaplot <- pairwise_termsim(Case_go)
  emapplot(Case_go_emaplot) + ggtitle("emaplot of GO analysis")
  # graph2ppt(file = "output/plots/Enrich_Analysis.pptx", width= 10, append = T, aspectr = 1.5)
#### 作好看的barplot of GO 
  {
    # sapply(data$Description, shorten_names)
    # test <- sapply(data$Description, shorten_names)
    # 
    # data <- Case_go@result   
    # data$ONTOLOGY<-gsub("BP", "Biological_Process", data$ONTOLOGY)
    # data$ONTOLOGY<-gsub("CC", "Cellular_Component", data$ONTOLOGY)
    # data$ONTOLOGY<-gsub("MF", "Molecular_Function", data$ONTOLOGY)
    # data$Description <- sapply(data$Description, shorten_names)
    # 
    # data <- separate(data = data, col = Term, into = c("GO_id", "GO_term"), sep = "~")
    # data <- subset(data,Count>3)  #数目很多时才做
    # data$GO_term = (sapply(levels(data$GO_term)[as.numeric(data$GO_term)],shorten_names))
    # data <- data[order(data[,1]),] #排序
    # data$GO_term<- as.character(data$GO_term)  #先转换成字符串
    # data$GO_term<-factor(data$GO_term,levels = c(data$GO_term)) #再强制加入因子
    # COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")
    # ggplot(data=data, aes(x=Description,y=Count, fill=ONTOLOGY)) + geom_bar(stat="identity", width=0.8) + coord_flip() +  xlab("GO term") + ylab("Num of Genes") +scale_fill_manual(values = COLS)+ theme_bw()
    # 
    
    
  } 

  
  
  
  
  
  


