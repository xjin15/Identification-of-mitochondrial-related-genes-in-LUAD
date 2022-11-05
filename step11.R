
# 1.圆圈弦图 ------------------------------------------------------------------
####弦图的数据准备
go <- data.frame(Category = Case_go$ONTOLOGY,
                 ID = Case_go$ID,
                 Term = Case_go$Description,
                 Genes = gsub("/", ", ", Case_go$geneID),
                 adj_pval = Case_go$p.adjust)
# 显示准备一个新的dataframe，其实就是改改列名
library(GOplot) # 一个专门搞GO分析美化的R包
gene_sig1 <- gene_sig |> rename(ID = symbol)## 还用到genesig，一列是symbol一列是logFC。记得改名
circ <- circle_dat(go, gene_sig1) # 创建弦图数据
head(circ)
termNum = 5         # 展示5个GO分析得到的terms                      
geneNum = nrow(gene_sig1)                      
chord <- chord_dat(circ, gene_sig1[1:geneNum,], go$Term[1:termNum]) # 

### 出图了
GOChord(chord,
        space = 0.02,        
        gene.order = 'logFC',    
        gene.space = 0.25,      
        gene.size = 2)
# 还可以，配色一般
graph2ppt(file = "output/plots/Enrich_Analysis.pptx", width= 10, append = T, aspectr = 1)

# 2.泡泡图 -------------------------------------------------------------------

# 不好看
GOBubble(circ, title = 'Bubble_plot',
         colour = c('skyblue', 'pink', 'red'),
         display = 'multiple', labels = 3)

reduced_circ <- reduce_overlap(circ, overlap = 0.60) # 去掉太多重叠部分
GOBubble(reduced_circ, labels = 2) # 泡泡图在左边，表格在右边


# 3.聚类树图和环形图 --------------------------------------------------------------
# 聚类树图，先选择你要画的terms
process<-c("mitotic nuclear division","chromatid segretion","organelle fission")
 
GOCluster(circ,
          process,
          clust.by = 'term')

# 环形图
GOCircle(circ,nsub=5)
# 热图，丑不拉几的
GOHeat(chord, nlfc = 1, fill.col = c('gold', 'white', 'purple'))

