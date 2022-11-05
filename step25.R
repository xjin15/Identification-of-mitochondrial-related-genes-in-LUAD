###### 需要做批量生存分析的基因挑出来
load(file = "outdata/step4.3_enrich_Analysis.rds")
DEG$genes
DEG_more <- all_diff |> dplyr::select(ID,logFC,adj.P.Val) |> 
  filter( abs(logFC) > 0.5 &
                      adj.P.Val < 0.05) 
names(DEG_more) <- c("genes","fold_change","p_value")
DEG_more$regulation <- "up"
DEG_more$regulation[DEG$fold_change<0] <- "down"
write.table(DEG_more,
            file = "outdata/deg_more_logfc.tsv",
            sep = "\t",
            )
table(DEG_more$regulation)

load(file = "outdata/step3_julei.rds")

intersect(DEG_more$genes,siggene)
deg_mito <- intersect(DEG_more$genes,siggene)
write.table(deg_mito,file = "outdata/deg_mito.tsv",quote = F,sep = "\t")

library(tinyarray)
data(meta1)

# 批量生存分析--by jimmy zeng -----------------------------------------------------------
rm(list=ls())
## 50 patients and 200 genes 
dat=matrix(rnorm(10000,mean=8,sd=4),nrow = 200)
rownames(dat)=paste0('gene_',1:nrow(dat))
colnames(dat)=paste0('sample_',1:ncol(dat))
os_months=abs(floor(rnorm(ncol(dat),mean = 50,sd=20)))
os_status=sample(rep(c('LIVING','DECEASED'),25))

library(survival)
my.surv <- Surv( os_months,os_status=='DECEASED')
## The status indicator, normally 0=alive, 1=dead. Other choices are TRUE/FALSE (TRUE = death) or 1/2 (2=death). 
## And most of the time we just care about the time od DECEASED . 

fit.KM=survfit(my.surv~1)
fit.KM
plot(fit.KM)

log_rank_p <- apply(dat, 1, function(values1){
  group=ifelse(values1>median(values1),'high','low')
  kmfit2 <- survfit(my.surv~group)
  #plot(kmfit2)
  data.survdiff=survdiff(my.surv~group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
})
names(log_rank_p[log_rank_p<0.05])

gender = ifelse(rnorm(ncol(dat))>1,'male','female')
age = abs(floor(rnorm(ncol(dat),mean = 50,sd=20)))
## gender and age are similar with group(by gene expression)

cox_results <- apply(dat , 1, function(values1){
  group=ifelse(values1>median(values1),'high','low')
  survival_dat <- data.frame(group=group,gender=gender,age=age,stringsAsFactors = F)
  m=coxph(my.surv ~ age + gender + group, data =  survival_dat)
  
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)
  HRse <- HR * se
  
  #summary(m)
  tmp <- round(cbind(coef = beta, se = se, z = beta/se, p = 1 - pchisq((beta/se)^2, 1),
                     HR = HR, HRse = HRse,
                     HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                     HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
  return(tmp['grouplow',])
  
})
cox_results[,cox_results[4,]<0.05]




# tcga批量生存分析 by 花花老师R包 tinyarray --------------------------------------------------------------------

rm(list = ls())
library(tinyarray)
data(exp_hub1)
data("exprSet_hub1")
data(meta1)
load(file = "outdata/step4.2_mrna.rds")
load(file = "outdata/LUAD_cln_clean.rds")

meta <- phe_f |> select(sample,OS,OS.time) |> 
  rename(event = OS,
         time = OS.time) 
meta <- meta[sample %in% colnames(mrna),]
exprSet_hub <- as.data.frame(mrna)
point_cut(exprSet_hub = exprSet_hub_nolog,meta = meta)

# 测试


exprSet_hub_nolog <- apply(exprSet_hub, 2, FUN = function(x){2^x - 1}) # 去log化




  if (ncol(exprSet_hub) != nrow(meta)) 
    stop("your exprSet_hub is not corresponds to meta")
  dat = cbind(t(exprSet_hub), meta)
  cut.point = c()
  for (i in 1:nrow(exprSet_hub)) {
    tryCatch({
      cut = survminer::surv_cutpoint(dat, 
                                     time = "time", 
                                     event = "event", 
                                     variables = rownames(exprSet_hub)[i])
      cut.point[[i]] = cut[["cutpoint"]][1, 1]
      },error=function(e){cat("ERROR :",conditionMessage(e), i,"\n")})
  }
  names(cut.point) = rownames(exprSet_hub)
  return(cut.point)

point_cut()



test_pointcut <- survminer::surv_cutpoint(data = dat,
                         time = "time",
                         event = "event",
                         variables = rownames(exprSet_hub)[1:200],
                         minprop = 0.1,
                         )
summary(test_pointcut)


meta$X_PATIENT <- substring(meta$sample,1,12)
surv_KM(exprSet_hub, meta, cut.point = T )

point_cut(exprSet_hub, meta)



group_list = make_tcga_group(exp_hub1);table(group_list)
## Warning in make_tcga_group(exp_hub1): NAs introduced by coercion
## group_list
## normal  tumor 
##    171    179
dim(exprSet_hub1)
## [1]   8 177
point_cut(exprSet_hub1,meta1)
surv_KM(exprSet_hub1,meta1)

surv_KM(exprSet_hub1,meta1,cut.point = T)
surv_cox(exprSet_hub1,meta1,cut.point = T)
k = exp_boxplot(log2(exp_hub1+1));k[[1]];length(k)
patchwork::wrap_plots(k,nrow = 2)
tmp = exp_surv(exprSet_hub1,meta1);length(tmp )
tmp
