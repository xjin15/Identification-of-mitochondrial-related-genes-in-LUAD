# 1.GSEA富集分析数据准备 ------------------------------------------------------------------
library(pacman)
p_load(fgsea,DOSE,tidyverse)
library(GOSemSim)

# 准备 geneList 文件 是降序的，有名字的数字型vector
## 这里用了管道符，先转成matrix，才能保留行名(R语言的基础通病)
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

# 2.GSEA的GO 和 KEGG --------------------------------------------------------

gse_go <- gseGO(
  geneList     = geneList,    #基因list
  OrgDb        = org.Hs.eg.db,  # 注释库
  ont          = "ALL",         #'ALL','BP','MF','CC'
  # keyType      = "SYMBOL",      # 默认为ENTREZID
  # nPerm      = 1000,       老版本的参数，现在已经不推荐使用了
  minGSSize    = 10,            #基因集最小
  maxGSSize    = 500,           #基因集最大
  pvalueCutoff = 0.05,          # p值阈值
  exponent     = 1,             # step权重默认即可
  by           = "fgsea",       # 使用的方法
  verbose      = T ,             # 过程可视化，随意选T
  eps          = 0,
)
# KEGG不像GO，KEGG必须把symbol转成ENTEREZID 
head(geneList)
gse_kegg <- gseKEGG(
  geneList     = geneList,      #基因list
  organism     = "hsa",         # 物种人
  keyType      = "kegg",        # 默认kegg就是ENTREZID
  # nPerm      = 1000,       现在已经不推荐使用了
  minGSSize    = 10,            #基因集最小
  maxGSSize    = 500,           #基因集最大
  pvalueCutoff = 0.05,          # p值阈值
  exponent     = 1,             # step权重默认即可
  by           = "fgsea",       # 使用的方法
  verbose      = T  ,            # 过程可视化，随意选T
  eps          = 0
)


# 3.msigdbr包准备gsea其他通路数据####

## 3.1熟悉msigdbr包 -----------------------------------------------------------
library(msigdbr)
msigdbr_species()

human <- msigdbr(species = "Homo sapiens")
human[1:4,1:4]
table(human[,1])
# H: hallmark gene sets 
# C1: positional gene sets 
# C2: curated gene sets 
# C3: motif gene sets 
# C4: computational gene sets 
# C5: GO gene sets 
# C6: oncogenic signatures 
# C7: immunologic signatures

table(human$gs_subcat)
#                     CGN             CGP              CM              CP     CP:BIOCARTA         CP:KEGG          CP:PID     CP:REACTOME 
# 214578           49342          426486           56738            4197            5515           16283            8722           99851 
# CP:WIKIPATHWAYS           GO:BP           GO:CC           GO:MF             HPO     IMMUNESIGDB  MIR:MIR_Legacy       MIR:MIRDB        TFT:GTRD 
# 33849          721379          115769          122674          464649         1066889           36362          397724          268672 
# TFT:TFT_Legacy             VAX 
# 170053           52075 


## 3.2选择GO和KEGG下载好 ----------------------------------------------------------------
KEGG_df = msigdbr(species = "Homo sapiens",category = "C2",subcategory = "CP:KEGG") %>% 
  dplyr::select(gs_exact_source,gene_symbol)
head(KEGG_df)
length(unique(KEGG_df$gene_symbol)) #基因数量
length(unique(KEGG_df$gs_exact_source)) #通路数量

GO_df = msigdbr(species = "Homo sapiens",category = "C5") %>% 
  dplyr::select(gene_symbol,gs_exact_source,gs_subcat)
dim(GO_df)

GO_df = GO_df[GO_df$gs_subcat!="HPO",]
table(GO_df$gs_subcat)

GO_df = GO_df[,c(2,1)]
head(GO_df)
# 基因数量
length(unique(GO_df$gene_symbol))
length(unique(GO_df$gs_exact_source))# terms数量

## 3.3做go和kegg -------------------------------------------------------------
ge = all_diff_sh1$logFC
names(ge) = rownames(all_diff_sh1)
ge = sort(ge,decreasing = T)
head(ge)
# head(geneList)
# 自己 获取的df是symbol，不是entrezid了
head(GO_df)
em <- GSEA(ge, TERM2GENE = GO_df)
em
#画个图来看看
gseaplot2(em, geneSetID = 1, title = em$Description[1])

# 4选择合适的category --------------------------------------------------------
# 肿瘤研究最常用的hallmarks
h <- msigdbr::msigdbr(species = "Homo sapiens", # 物种拉丁名
             category = "H") #此处以hallmark为例，你也可以选择MSigDB的其他注释


h <- dplyr::select(h, gs_name, gene_symbol) %>% #或entrez_gene
  as.data.frame %>% 
  split(., .$gs_name) %>%    # split 函数
  lapply(., function(x)(x$gene_symbol)) #或entrez_gene

# 在每个geneset里面去掉重复的基因
gs <- lapply(h, unique)

# 接下来去掉那些在两个或更多个pathways里出现过的genes
count <- table(unlist(gs))
keep <- names(which(table(unlist(gs)) < 2))
gs <- lapply(gs, function(x) intersect(keep, x))

# 过滤之后，很多pathway一个gene都不剩了，去掉这些
gs <- gs[lapply(gs, length) > 0]

# 预览过滤后的结果
head(gs)
save(gs, file = "outdata/hallmark.gs.RData")


#### 还没弄好代码
gse_h <- GSEA(
  geneList     = geneList,      #基因list
  TERM2GENE = NA,
  # TERM2NAME = gs,
  # nPerm      = 1000,       现在已经不推荐使用了
  minGSSize    = 10,            #基因集最小
  maxGSSize    = 500,           #基因集最大
  pvalueCutoff = 0.05,          # p值阈值
  exponent     = 1,             # step权重默认即可
  by           = "fgsea",       # 使用的方法
  verbose      = T  ,            # 过程可视化，随意选T
  eps          = 0
)


# 4.直接使用gsea官网下载的gmt文件 ----------------------------------------------------
MSigDB <- read.gmt(gmtfile = "../0000DATA/msigdb.v7.5.1.symbols.gmt")


GSEA.all <- GSEA(geneList = ge,
                  TERM2GENE=MSigDB,
                  # nPerm = 1000,
                  # minGSSize = minGSSize.gsea,
                  # maxGSSize = maxGSSize.gsea,
                  seed = T,
                  verbose = F,
                  pvalueCutoff = 0.05) 









# 5.GSEA结果的可视化######
write.csv(x=geneList, file = "outdata/genelistforgsea.txt")
dotplot(gse_go,
        showCategory = 10,
        split = "ONTOLOGY",
        label_format = 25,
        title = "Dotplot of GO GSEA analysis") + facet_grid(ONTOLOGY ~ ., scale="free") 
graph2ppt(file = "output/plots/GSEA_Analysis.pptx", width= 12, append = T, aspectr = 1.5)
# 山脊图
ridgeplot(gse_go,
          showCategory = 30,   # 展示多少个pathway
          fill = "p.adjust",   # 可以选用"pvalue", "p.adjust", "qvalue"
          core_enrichment = T, # 是否只展示核心的富集基因
          label_format = 20,   # 通路的名字超过20个字符就换行
          orderBy = "p.adjust",     # Y轴的排序
          decreasing = F       # 升序
)   + ggtitle(  "Ridge plot of GO GSEA analysis")
graph2ppt(file = "output/plots/GSEA_Analysis.pptx", width= 10, append = T, aspectr = 1.5)



##Y树自己画ggplot2
x <- gseGO(geneList,ont = "BP",OrgDb = org.Hs.eg.db,keyType = "SYMBOL")
x1 <- pairwise_termsim(x2)
# simplify 只能用于GOterm，而且一次只能用一个"BP""CC""MF".
# 因为这是其他某些搞语义的人精简的。
x2 <- simplify(x = x, 
               cutoff=0.7, 
               by="p.adjust", 
               select_fun=min,
               measure = "Wang",
               semData = NULL) 
# 对这个例子，x本来有2100多terms，
# simplify以后只剩下600多个terms 
# 但是感觉这儿simplify还是比较鸡肋的。
# 首先，只能用于GO的terms
# 其次 只能适用于一个一个的terms，比如ont =“BP”
emapplot(x1)
x3 <- pairwise_termsim(x)
emapplot(x)
upsetplot(x2)
dotplot(x2);dotplot(x)
cnetplot(x2)



data(gcSample)
x = compareCluster(gcSample, fun='enrichGO', ont='MF', OrgDb='org.Hs.eg.db')
y <- simplify(x, measure='Wang', semData=NULL) 


y <- arrange(x, abs(NES)) %>% # 先是排序，对x以NES的绝对值排序
  group_by(sign(NES)) %>%    # sign()函数根据数字正负返回值，-1，0，1
  slice(1:5)                 # 选取前5行
library(forcats)
library(ggplot2)
library(ggstance)
library(enrichplot)

ggplot(y, aes(NES, fct_reorder(Description, NES), fill=qvalues), showCategory=10) + 
  geom_col(orientation='y') + 
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  theme_minimal() + ylab(NULL)

##

gseaplot2(
  x = gse_go,
  geneSetID = 1:5,                 # 选择要展示的genesetid
  title = " GSEA plot of GO analysis",                     # 图片的标题
  color = "green",                 # 颜色选择
  base_size = 11,                  # 基本字体大小选择
  pvalue_table = F,                # 是否显示p值的表格
  ES_geom = "line",                # ES值的表现形式，dot或者line
  subplots = c(1,3),                  # 显示3张小图中的哪几张
  rel_heights = c(1.5, 0.5, 1)     # 3张图的相对高度
)
graph2ppt(file = "output/plots/GSEA_Analysis.pptx", width= 10, append = T, aspectr = 1.5)
tmp <- gse_go@result


dotplot( gse_kegg,
         showCategory = 10,
         label_format = 20,
         title = "Dotplot of  KEGG GSEA analysis") 
graph2ppt(file = "output/plots/GSEA_Analysis.pptx", width= 10, append = T, aspectr = 1.5)
# 山脊图
ridgeplot(gse_kegg,
          showCategory = 30,   # 展示多少个pathway
          fill = "p.adjust",   # 可以选用"pvalue", "p.adjust", "qvalue"
          core_enrichment = T, # 是否只展示核心的富集基因
          label_format = 20,   # 通路的名字超过20个字符就换行
          orderBy = "NES",     # Y轴的排序
          decreasing = F       # 不降序,升序，p值小的在前
) + ggtitle(  "Ridge plot of KEGG GSEA analysis")
graph2ppt(file = "output/plots/GSEA_Analysis.pptx", width= 10, append = T, aspectr = 1.5)


gseaplot2(
  x = gse_kegg,
  geneSetID = 1:5,                 # 选择要展示的genesetid
  title = " GSEA plot of KEGG analysis",                     # 图片的标题
  color = "green",                 # 颜色选择
  base_size = 11,                  # 基本字体大小选择
  pvalue_table = F,                # 是否显示p值的表格
  ES_geom = "line",                # ES值的表现形式，dot或者line
  subplots = 1:3,                  # 显示3张小图中的哪几张
  rel_heights = c(1.5, 0.5, 1)     # 3张图的相对高度
)
graph2ppt(file = "output/plots/GSEA_Analysis.pptx", width= 10, append = T, aspectr = 1.5)
tmp <- gse_kegg@result

