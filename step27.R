rm(list = ls())

# 1.准备数据 ------------------------------------------------------------------


# 做生存的lasso

# 需要TCGA表达矩阵expr、临床信息meta、分组group_list
load('outdata/LUAD_TCGA_fpkm_symbol_maxgene.Rdata')
load('outdata/step4_group_data.rds')
load('outdata/LUAD_cln_clean.rds')
exprSet = fpkm[,sample_group$sample]
meta <- phe_f %>% 
  select(sample,OS,OS.time) %>% 
  rename(ID = sample,
         event = OS) %>% as.data.frame()
rownames(meta) <- meta$ID
meta <- meta[sample_group$sample,]
## 必须保证生存资料和表达矩阵，两者一致
all(colnames(exprSet)==meta$ID)

# 2.构建lasso回归模型 -----------------------------------------------------------
# lasso回归的输入数据是表达矩阵(仅含tumor样本)和对应的生死
x=t(exprSet)
y=meta$event
library(glmnet)
set.seed(123)
model_lasso <- glmnet(x, y,nlambda=10, alpha=1)
print(model_lasso)
#> 
#> Call:  glmnet(x = x, y = y, alpha = 1, nlambda = 10) 
#> 
#>        Df    %Dev   Lambda
#>  [1,]   0 0.00000 0.127800
#>  [2,]   7 0.08689 0.076610
#>  [3,]  22 0.19490 0.045930
#>  [4,]  49 0.28660 0.027530
#>  [5,]  97 0.42080 0.016510
#>  [6,] 164 0.55460 0.009895
#>  [7,] 217 0.66240 0.005932
#>  [8,] 275 0.74850 0.003556
#>  [9,] 345 0.81660 0.002132
#> [10,] 403 0.86840 0.001278
# 这里是举例子，所以只计算了10个λ值，解释一下输出结果三列的意思。
# Df 是自由度
# 列%Dev代表了由模型解释的残差的比例，对于线性模型来说就是模型拟合的R^2(R-squred)。它在0和1之间，越接近1说明模型的表现越好，如果是0，说明模型的预测结果还不如直接把因变量的均值作为预测值来的有效。
# Lambda 是构建模型的重要参数。
# 解释的残差百分比越高越好，但是构建模型使用的基因的数量也不能太多，需要取一个折中值。


# 2.1 构建交叉验证挑合适的λ值 --------------------------------------------------------
cv_fit <- cv.glmnet(x=x, y=y, nlambda = 1000,alpha = 1)
plot(cv_fit)
# 两条虚线分别指示了两个特殊的λ值,
# 一个是lambda.min,一个是lambda.1se,
# 这两个值之间的lambda都认为是合适的。
# lambda.1se构建的模型最简单，即使用的基因数量少，
# 而lambda.min则准确率更高一点，使用的基因数量更多一点。
library(export)
library(eoffice)
eoffice::topptx(filename = 'output/plots/LASSO.pptx',append = T,)
# 2.2 用这两个λ值重新建模 ----------------------------------------------------------
model_lasso_min <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.min)
model_lasso_1se <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.1se)








