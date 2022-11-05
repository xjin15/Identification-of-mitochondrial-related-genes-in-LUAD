rm(list = ls())
# LASSO挑选合适的基因作为signature -------------------------------------------------
options(digits = 4)
p_load(ggplot2,ggthemes,ggpubr,ggrepel,export,tidyverse,glmnet,corrplot)
# 数据准备 --------------------------------------------------------------------
load(file = "outdata/step4.2_mrna.rds")
load(file = "outdata/step4_group_data.rds")
load(file = "outdata/LUAD_TCGA_ALLfpkm.rds")
rm(allid497)
allid <- as.data.frame(allid)
allid$X <- "clusterB"
allid$X[1:279] <- "clusterA" 
library(corrplot) #相关系数分析
library(glmnet) # 岭回归、 LASSO 、弹性网络模型
library(MASS) 

# fpkm和mrna转成行为样本，列为基因及其表达量的matrix。

mrna <- na.omit(mrna)
fpkm <- na.omit(fpkm)
mrna <- t(mrna)
fpkm <- t(fpkm)
mrna[1:4,1:4]
fpkm[1:4,1:4]

set.seed(123) #random number generator 
x1 <- as.matrix(mrna)  # 预测变量为多个因素，放在x中
x2 <- as.matrix(fpkm)
x2 <- x2[rownames(x1),]
y <- as.factor(allid$X)

load(file = "outdata/step4.3_enrich_Analysis.rds")
DEG_more <- read.table(file = "outdata/deg_more_logfc.tsv",
                       row.names = 1)
colnames(DEG_more) <- DEG_more[1,]
DEG_more <- DEG_more[-1,]
x <- x1[ ,rownames(DEG_more)]


######## 10折交叉验证挑选合适的λ值 #######
# 交叉验证

# cvlasso1 <- cv.glmnet(x = x1,y = y,
#                       family = "binomial", # 因变量是二分类
#                       alpha = 1, # lasso回归
#                       # type.measure = "auc", # 也可以不选
#                       nfolds = 10,)
# plot(cvlasso1,xvar = "lambda", label = T) # 回归系数随着logλ的变化曲线
# plot(cvlasso1)
# # 用整个mrna来预测，发现不行,选出来太多变量了
# cvlasso2 <- cv.glmnet(x = x2,y = y,
#                      family = "binomial", # 因变量是二分类
#                      alpha = 1, # lasso回归
#                      # type.measure = "auc", # 也可以不选
#                      nfolds = 10,)
# plot(cvlasso2,xvar = "lambda", label = T) # 回归系数随着logλ的变化曲线
# plot(cvlasso2)
# 用整个fpkm来预测，发现也不行,选出来太多变量了
x <- x1[, DEG$genes]
cvlasso <- cv.glmnet(x = x,y = y,
                     family = "binomial", # 因变量是二分类
                     alpha = 1, # lasso回归
                     type.measure = "class", # 也可以不选
                     nfolds = 5,)
plot(cvlasso,
     # xvar = "lambda", label = T
     ) # 回归系数随着logλ的变化曲线
lasso_deg <- glmnet(x = x,y = y,family = "binomial",alpha = 1,standardize = T,)
plot(lasso_deg,xvar = "lambda")

model_lasso_min <- glmnet(x = x, y= y,
                          alpha = 1,
                          lambda=cvlasso$lambda.min,
                          family = "binomial")
model_lasso_1se <- glmnet(x = x, y= y,
                          alpha = 1,
                          lambda=cvlasso$lambda.1se,
                          family = "binomial")
head(model_lasso_min$beta)
head(model_lasso_1se$beta)
choose_gene_min <- rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]
choose_gene_1se <- rownames(model_lasso_1se$beta)[as.numeric(model_lasso_1se$beta)!=0]

betamin <- model_lasso_min$beta[as.numeric(model_lasso_min$beta)!=0]
beta1se <- model_lasso_1se$beta[as.numeric(model_lasso_1se$beta)!=0]

length(choose_gene_1se);choose_gene_1se;beta1se
length(choose_gene_min);choose_gene_min;betamin
beta1se
graph2ppt(file = "output/plots/LASSO.pptx", append = T)
eoffice::topptx(filename = "output/plots/LASSO.pptx",append = T)
# 
# # 用差异基因more1747个预测，选出来74多个基因，不行 和线粒体相关基因做交叉
# DEG_mitogene <- intersect(rownames(DEG_more),mtgene1626$Symbol)
# x <- x[,DEG_mitogene]
# cvlasso <- cv.glmnet(x = x,y = y,
#                      family = "binomial", # 因变量是二分类
#                      alpha = 1, # lasso回归
#                      type.measure = "class", # 也可以不选
#                      nfolds = 5,)
# plot(cvlasso)
# 
# # 还是有很多，继续筛选DEG_MITOGENE_SURVIVAL
# load(file = "outdata/step3_julei.rds")
# # siggene 这里的siggene是线粒体基因1626个中，在TCGA肿瘤和癌旁组织有差异的基因
# DEG_mitogene <- intersect(DEG_mitogene,siggene)
# # 再筛选有生存差异的
# DEsurv <- data.table::fread(file = "outdata/所有基因在LUAD肿瘤和癌旁的差异.tsv")
# DEsurv_sig <- DEsurv |> filter(p.adj < 0.05)
# DEG_mitogene_desurv <- intersect(DEG_mitogene,DEsurv_sig$DEGgene)
# x_final <- x[,DEG_mitogene_desurv]
# cvlasso <- cv.glmnet(x = x_final,y = y,
#                      family = "binomial", # 因变量是二分类
#                      alpha = 1, # lasso回归
#                      type.measure = "class", # 也可以不选
#                      nfolds = 5,)
# plot(cvlasso)
# 
# # 39个特征 
# model_lasso_min <- glmnet(x = x, y= y, 
#                           alpha = 1, 
#                           lambda=cvlasso$lambda.min,
#                           family = "binomial")
# model_lasso_1se <- glmnet(x = x, y= y, 
#                           alpha = 1, 
#                           lambda=cvlasso$lambda.1se,
#                           family = "binomial")
# head(model_lasso_min$beta)
# head(model_lasso_1se$beta)
# choose_gene_min <- rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]
# choose_gene_1se <- rownames(model_lasso_1se$beta)[as.numeric(model_lasso_1se$beta)!=0]
# 
# betamin <- model_lasso_min$beta[as.numeric(model_lasso_min$beta)!=0]
# beta1se <- model_lasso_1se$beta[as.numeric(model_lasso_1se$beta)!=0]
# 
# length(choose_gene_1se);choose_gene_1se;beta1se
# length(choose_gene_min);choose_gene_min;betamin
# 
# choose_gene_min |> grep(pattern = "PNPT1",value = T,)
# choose_gene_min |> grep(pattern = "DA",value = T,)
# 
# cvlasso_2 <- cv.glmnet(x = x[,choose_gene_min],y = y,
#                        family = "binomial", # 因变量是二分类
#                        alpha = 1, # lasso回归
#                        type.measure = "class", # 也可以不选
#                        nfolds = 5,)
# plot(cvlasso_2)
# #
# cvlasso_2$lambda.1se  # 显示λ的值
# model_lasso_min_2 <- glmnet(x = x[,choose_gene_min], y= y, 
#                           alpha = 1, 
#                           lambda=cvlasso$lambda.min,
#                           family = "binomial")
# choose_gene_min_2 <- rownames(model_lasso_min_2$beta)[as.numeric(model_lasso_min_2$beta)!=0]
# choose_gene_min_2
# choose_gene_min_2 |> grep(pattern = "PNPT1",value = T,)
# choose_gene_min_2 |> grep(pattern = "DA",value = T,)

write.table(choose_gene_min,file = "outdata/lasso分组相关基因.txt",quote = F,row.names = F)
##### 自己预测自己 ######

lasso.prob <- predict(cvlasso, newx=x , s=c(cvlasso$lambda.min,cvlasso$lambda.1se) )
re=cbind(y,lasso.prob)
head(re)
lasso.prob
library(ROCR)
library(caret)
# 自己预测自己
#min
pred_min <- prediction(re[,2], re[,1])
auc_min = performance(pred_min,"auc")@y.values[[1]]
perf_min <- performance(pred_min,"tpr","fpr")
plot(perf_min,colorize=FALSE, col="blue") 
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
text(0.8,0.2, labels = paste0("AUC = ",round(auc_min,3)))
#1se
pred_1se <- prediction(re[,3], re[,1])
auc_1se = performance(pred_1se,"auc")@y.values[[1]]
perf_1se <- performance(pred_1se,"tpr","fpr")
plot(perf_1se,colorize=FALSE, col="red") 
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
text(0.8,0.2, labels = paste0("AUC = ",round(auc_min,3)))

plot(perf_min,colorize=FALSE, col="blue") 
plot(perf_1se,colorize=FALSE, col="red",add = T) 
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
text(0.8,0.3, labels = paste0("AUC of min = ",round(auc_min,3)),col = "blue")
text(0.8,0.2, labels = paste0("AUC of 1se = ",round(auc_1se,3)),col = "red")
library(export)
graph2ppt(file = "output/plots/LASSO.pptx", append = T)
########## 分成两组，训练集和验证集 ##########
library(caret)
set.seed(110)

sam <- createDataPartition(allid$X, p = .5,list = FALSE)
head(sam)

train <- x[sam,]
train_fapkm <- fpkm[sam,]
test <- x[-sam,]
train_meta <- allid[sam,]
test_meta <- allid[-sam,]
x_train = train
y_train = as.factor(train_meta$X)

cvlasso <- cv.glmnet(x = x,y = y,
                     family = "binomial", # 因变量是二分类
                     alpha = 1, # lasso回归
                     # type.measure = "auc", # 也可以不选
                     nfolds = 10,)

cv_fit <- cv.glmnet(x = x_train,y = y_train,
                    family = "binomial", # 因变量是二分类
                    alpha = 1, # lasso回归
                    # type.measure = "auc", # 也可以不选
                    nfolds = 10,)

plot(cv_fit)
graph2ppt(file = "output/plots/LASSO.pptx", append = T)

model_lasso_min <- glmnet(x=x, y=y, alpha = 1, family = "binomial",lambda=cv_fit$lambda.min)
model_lasso_1se <-glmnet(x=x, y=y, alpha = 1, family = "binomial",lambda=cv_fit$lambda.1se)
head(model_lasso_min$beta)
choose_gene_min2=rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0] %>% 
  substring(1,15)
choose_gene_1se2=rownames(model_lasso_1se$beta)[as.numeric(model_lasso_1se$beta)!=0] %>% 
  substring(1,15)
length(choose_gene_min2)
length(choose_gene_1se2)

lasso.prob2 <- predict(cv_fit, newx = test, type="response", 
                       s=c(cv_fit$lambda.min,cv_fit$lambda.1se)
                       )
re=cbind(as.factor(test_meta$X) ,lasso.prob2)
head(re)
# AUC图
# 训练集模型预测测试集
#min
pred_min <- prediction(re[,2], re[,1])
auc_min = performance(pred_min,"auc")@y.values[[1]]
perf_min <- performance(pred_min,"tpr","fpr")
plot(perf_min,colorize=FALSE, col="blue") 
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
text(0.8,0.2, labels = paste0("AUC = ",round(auc_min,3)))
#1se
pred_1se <- prediction(re[,3], re[,1])
auc_1se = performance(pred_1se,"auc")@y.values[[1]]
perf_1se <- performance(pred_1se,"tpr","fpr")
plot(perf_1se,colorize=FALSE, col="red") 
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
text(0.8,0.2, labels = paste0("AUC = ",round(auc_min,3)))

plot(perf_min,colorize=FALSE, col="blue") 
plot(perf_1se,colorize=FALSE, col="red",add = T) 
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
text(0.8,0.3, labels = paste0("AUC of min = ",round(auc_min,3)),col = "blue")
text(0.8,0.2, labels = paste0("AUC of 1se = ",round(auc_1se,3)),col = "red")
library(export)
graph2ppt(file = "output/plots/LASSO.pptx", append = T)
########### 用筛选到的特征基因做logistic回归分析#########
gtf_gene <- read.csv(file = "outdata/gtfmRNA22.txt", header = T,sep = "\t")
choose_gene_1se2 <- gtf_gene$gene_name[match(choose_gene_1se2,gtf_gene$gene_id,nomatch = 0)]
choose_gene_min2 <- gtf_gene$gene_name[match(choose_gene_min2,gtf_gene$gene_id,nomatch = 0)]
length(choose_gene_1se);choose_gene_1se;beta1se
length(choose_gene_min);choose_gene_min;betamin
length(choose_gene_1se2);choose_gene_1se2
length(choose_gene_min2);choose_gene_min2

# 最终选择 min2 的19个基因表达量为突变分组条件
choose_gene_min2
choose_gene_min2_noid=rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]
choose_gene_min2_noid
options(digits = 3)
betamin2 <- model_lasso_min$beta[as.numeric(model_lasso_min$beta)!=0] 
betamin2
fpkm[1:4,1:4]
## 写公式，算分数 #####
formular_call = paste(round(betamin2,digits = 4), choose_gene_min2_noid,sep = " * ",collapse= ") + (");formular_call
# 便于复制到word中
paste(round(betamin2,digits = 4), choose_gene_min2,sep = "×",collapse = "）+（")
paste(choose_gene_min2,collapse = "、")

## 计算fpkm中的score
fpkm[1:4,1:4]
fpkmlogi <- fpkm %>% 
  transmute(score = (-0.1058 * ENSG00000246526.2) + (-0.1257 * ENSG00000151632.15) + (-48.6069 * ENSG00000223421.1) + (-0.357 * ENSG00000108448.19) + (0.0587 * ENSG00000131759.16) + (0.0244 * ENSG00000130766.4) + (-19.895 * ENSG00000229755.1) + (-0.3798 * ENSG00000249199.1) + (-0.4665 * ENSG00000275788.1) + (-0.2482 * ENSG00000239035.1) + (-0.1185 * ENSG00000251264.1) + (-0.32 * ENSG00000224928.2) + (-0.1907 * ENSG00000177156.9) + (-0.3089 * ENSG00000235447.2) + (-1.0802 * ENSG00000243072.1) + (-0.4102 * ENSG00000254339.4) + (-1.3136 * ENSG00000256474.1) + (-0.5642 * ENSG00000236929.2) + (-0.3232 * ENSG00000281606.1) )

fpkmlogi <- fpkm %>% 
  as.data.frame() %>% 
  # rownames_to_column("ID") %>% 
  transmute(score = (-0.1058 * ENSG00000246526.2) + (-0.1257 * ENSG00000151632.15) 
            + (-48.6069 * ENSG00000223421.1) + (-0.357 * ENSG00000108448.19) +
              (0.0587 * ENSG00000131759.16) + (0.0244 * ENSG00000130766.4) + 
              (-19.895 * ENSG00000229755.1) + (-0.3798 * ENSG00000249199.1) + 
              (-0.4665 * ENSG00000275788.1) + (-0.2482 * ENSG00000239035.1) + 
              (-0.1185 * ENSG00000251264.1) + (-0.32 * ENSG00000224928.2) + 
              (-0.1907 * ENSG00000177156.9) + (-0.3089 * ENSG00000235447.2) + 
              (-1.0802 * ENSG00000243072.1) + (-0.4102 * ENSG00000254339.4) + 
              (-1.3136 * ENSG00000256474.1) + (-0.5642 * ENSG00000236929.2) + 
              (-0.3232 * ENSG00000281606.1)) 
fpkmlogi$group <- as.factor(allid$X)
# 建logistic模型
library(pROC)
logiroc <- roc(fpkmlogi$group,fpkmlogi$score) 

plot.roc(logiroc,print.thres=T,print.auc =T,main="ROC曲线",col="#008600")
plot(logiroc,
     print.thres=TRUE,
     main="ROC曲线",
     col="#008600")
graph2ppt(file = "output/plots/LASSO专利", append = T)


