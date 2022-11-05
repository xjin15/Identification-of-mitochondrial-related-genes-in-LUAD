do::rm_all()
load(file = "outdata/step5_sur+vars.rds")
p_load(ggplot2, tidyverse, survival, survminer,export)


# 单因素COX分析 -------------------------------------------------------------------
#### 数据准 备 #####
library(mice)
md.pattern(dt_sur)# 查看缺失值的分布
# 一种具有ncol（x）+ 1列的矩阵，其中每一行对应一个缺失的数据模式（1=观察到的，0=缺失的）
# 行和列按丢失信息量的增加进行排序。最后一列和行分别包含行和列计数。

library(VIM)
aggr(dt_sur)
complete.cases(dt_sur) %>% table()
# 
# dt_cox <- dt_sur |> 
#   mutate(
#     Radiotherapy = ifelse(is.na(Radiotherapy),"No/Unknown","Yes"),
#     Age = ifelse(is.na(Age),65, Age),
#     Age_group = ifelse(is.na(Age),"60-70", Age_group),
#     
#     )

dt_cox <- within(dt_sur, {
  Radiotherapy[is.na(Radiotherapy)] <- "YES" #NO/Unknow
  Age[is.na(Age)] <- median(Age,na.rm = T)
  Smoke[is.na(Smoke)] <- "NO"
  NNstage <- NULL
  NNstage[Nstage=="N2"|Nstage=="N3"|Nstage=="N1"] <- "N123" 
  NNstage[Nstage=="N0"] <- "N0" 
    })
complete.cases(dt_cox) %>% table()

dt_cox <- na.omit(dt_cox)


dt_cox <- dt_cox %>% 
  mutate(Sex = factor(Gender, levels = c("female","male"),labels = c("Female","Male")),
         group = factor(group, levels = c(1,2),labels = c("A","B")),
                  Age_median = factor(Age_median, levels = c("younger","older"),labels = c("Younger","Older")),
         Age_group = factor(Age_group,levels = c("<60", "60-70",">70")),
                  )
str(dt_cox)




# 单因素cox分析 ----------------------------------------------------------------


phe_f <- dt_cox |> 
  relocate(Sex,.before = Gender) |> 
  relocate(NNstage,.after = Nstage) |> 
  select(!c(Gender,sample)) |> 
  select(OS:PFS.time, everything())



data1 <- phe_f %>% filter(PFS.time <= 84) |> as.data.frame()
data2 <- phe_f %>% filter(OS.time <= 84) |> as.data.frame()
#1.构建函数
pfs_surv<- Surv(time = data1$PFS.time, event = data1$PFS==1)
data1$surv <- with(data1,pfs_surv)  
# data1$surv1 <- pfs_surv #这两个是等价的
cox <- coxph(surv ~ group, data = data1)

summary(cox)

Unicox<- function(x){ 
  FML <- as.formula(paste0 ("surv~",x))
  cox <- coxph(FML,data=unidata)
  cox1 <- summary(cox)
  # beta <- round(cox1$coef[,1],3) # betaβ值就是coef偏回归系数
  HR <- round(cox1$coefficients[,2],2)  # HR 就是exp(coef)  
  PValue <- round(cox1$coefficients[,5],3)  # p值就是Pr(>|z|)
  xname <- rownames(cox1$conf.int)
  CI <- paste(round(cox1$conf.int[,3],2),round(cox1$conf.int[,4],2), sep = '-')  # 95%可信区间就是lower .95和 upper .95
  Unicox<- data.frame('Variables' = xname,
                      'Hazard Ratio' = HR,
                      'P Value' = PValue,
                      'CI95' = CI,
                      'HR(95% CI)'=paste(HR,"(",CI,")",sep = ""))
  return(Unicox)}     
###### PFS单因素#####
#2.提取变量，构建数据框，批量生成表格
variable.names<- colnames(data1)[which(names(data1)=="group"):(ncol(data1)-1)] 
variable.names

unidata <- data1
pfs_Univar <- lapply(variable.names, Unicox)

pfs_Univar <- plyr::ldply(pfs_Univar,data.frame)
pfs_Univar #查看表格结果
#3 一些修饰
pfs_Univar$P.Value[which(pfs_Univar$P.Value == 0 )] <- "<0.001"
pfs_Univar$p1 <- "DELETE"
pfs_Univar$p1[which(pfs_Univar$P.Value <= 0.1) ] <- "<0.1"

pfs_Univar
# 挑选出pfs需要的因素
write.csv(pfs_Univar,"outdata/PFS单因素回归分析表.csv")



##用data2 做一遍 OS的分析
######### OS 单因素#########
os_surv<- Surv(time = data2$OS.time, event = data2$OS==1)
data2$surv <- with(data2,os_surv)  
cox <- coxph(surv~group, data=data2)
summary(cox)

variable.names<- colnames(data2)[which(names(data2)=="group"):(ncol(data2)-1)] 
variable.names
unidata <- data2
os_Univar <- lapply(variable.names, Unicox)
os_Univar <- plyr::ldply(os_Univar,data.frame)
os_Univar$P.Value[which(os_Univar$P.Value == 0 )] <- "<0.001"
os_Univar$p1 <- "DELETE"
os_Univar$p1[which(os_Univar$P.Value <= 0.1) ] <- "<0.1"
os_Univar 
write.csv(os_Univar,"outdata/OS单因素回归分析表.csv")


# OS多因素cox ------------------------------------------------------------------
os_multicox_variable <- paste( variable.names,collapse = " + ")
os_multicox_formula <- paste0("coxph(Surv(OS.time, OS) ~ ",os_multicox_variable, ", data = data2)")
os_multicox <- eval(parse(text = os_multicox_formula))

os_multicox_step_formula <- step(os_multicox,direction = "ba")
os_multicox_step_formula <- summary(os_multicox_step_formula)
os_multicox_step_formula <- os_multicox_step_formula$call
os_multicox_step_formula
if (F) {
  os_multicox <- eval(os_multicox_step_formula)
}
os_multicox <- coxph(formula = Surv(OS.time, OS) ~ group + Age_group + Smoke + 
                       Tstage + NNstage + Radiotherapy, data = data2)
summary(os_multicox)  
ggforest(model = os_multicox, data = data2, main ="Forest of OS multicox Analysis")
graph2ppt(file="output/plots/OS.pptx", width = 20,aspectr = 1,append = T)

#2 提取HR.P.95%CI
os_multi1<-tableone::ShowRegTable(os_multicox, 
                        exp=TRUE, 
                        digits=2, 
                        pDigits =3,
                        printToggle = TRUE, 
                        quote=FALSE, 
                        ciFun=confint)
#3 提取回归系数、统计量等                     
os_multi2<-broom::tidy(os_multicox)
os_multi2
#4 将两次提取结果合并
os_multi<-cbind(os_multi1, os_multi2)
os_multi
write.csv(os_multi,file="outdata/多因素coxtable.csv")

##### PFS 多因素cox########
pfs_multicox_variable <- paste( variable.names,collapse = " + ")
pfs_multicox_formula <- paste0("coxph(Surv(PFS.time, PFS) ~ ",pfs_multicox_variable, ", data = data1)")
pfs_multicox <- eval(parse(text = pfs_multicox_formula))

pfs_multicox_step_formula <- step(pfs_multicox,direction = "ba")
pfs_multicox_step_formula <- summary(pfs_multicox_step_formula)
pfs_multicox_step_formula <- pfs_multicox_step_formula$call
pfs_multicox_step_formula

# pfs_multicox <- coxph(Surv(PFS.time, PFS) ~ group+Age+Age_median+Age_group+Sex+
#                         Race+Smoke+Tstage+Nstage+NNstage+Mstage+
#                         Radiotherapy+Tumor_site, data = data1)
# step(pfs_multicox)
if(F){
  pfs_multicox <- eval(pfs_multicox_step_formula)
}
pfs_multicox <- coxph(Surv(PFS.time, PFS) ~ group + Age_median + Tstage + 
                        NNstage + Radiotherapy, data = data1)
summary(pfs_multicox)  
ggforest(model = pfs_multicox, data = data1, main ="Forest of PFS multicox Analysis")
graph2ppt(file="output/plots/OS.pptx", width = 20,aspectr = 1,append = T)

#2 提取HR.P.95%CI
pfs_multi1 <- tableone::ShowRegTable(pfs_multicox, 
                         exp=TRUE, 
                         digits=2, 
                         pDigits =3,
                         printToggle = TRUE, 
                         quote=FALSE, 
                         ciFun=confint)
#3 提取回归系数、统计量等                     
pfs_multi2<-tidy(pfs_multicox) #broom包
pfs_multi2
#4 将两次提取结果合并
pfs_multi<-cbind(pfs_multi1, pfs_multi2)
pfs_multi
write.csv(pfs_multi,file="outdata/PFS多因素回归分析表.csv")

# nomogram制作 --------------------------------------------------------------

# 数据整理，不能有缺失值。
dt1 <- data1[-ncol(data1)]
dt1 <- na.omit(dt1)
uu <- dt1
y<-Surv(uu$PFS.time,uu$PFS==1)
kmfit<-survfit(y~1, data=y)
coxmodel <- coxph(y ~ group+Age+Age_median+Age_group+Sex+
                    Race+Smoke+Tstage+Nstage+NNstage+Mstage+
                    Radiotherapy+Tumor_site, data=uu)
step(coxmodel) ### 这一步是挑选合适的因素进行cox多因素分析，但其实我在这之前已经挑选过了
coxmodel <- coxph(y ~ group + Age_median + Tstage + NNstage + Radiotherapy, data=uu)
summary(coxmodel)
library(rms)
dd<-datadist(uu)
options(datadist="dd")

f <- cph( y ~ group + Age_median + Tstage + NNstage + Radiotherapy, x=T, y=T,surv=T, data=uu)
surv<-Survival(f)
nom1 <- nomogram(f, fun=list(function(x) surv(12, x), function(x) surv(36, x), function(x) surv(60, x)), lp=F,funlabel=c("1-year survival ", "3-year survival ", "5-year survival "), maxscale=10, fun.at=c(0.95, 0.9, 0.8, 0.7, 0.6, 0.5,0.4,0.3,0.2,0.1,0.05))

plot(nom1,cex.axis = 1.3,cex.var =1.5,xfrac=.4,lmgp=0.5,tcl=-0.5)
graph2ppt(file="output/plots/nomo.pptx", width=10,height=10,append=T)

# 内部验证 
validate(f, method="boot", B=1000, dxy=T)
rcorrcens(Surv(PFS.time, PFS) ~ predict(f), data = uu) #
# Somers' Rank Correlation for Censored Data    Response variable:Surv(PFS.time, PFS)
#                C    Dxy  aDxy    SD    Z P   n
# predict(f) 0.356 -0.288 0.288 0.043 6.75 0 428
# C—index = 1-C==1-0.356=0.644


# 一致性检验
# 1年
f1 <- cph( y ~ group + Age_median + Tstage + NNstage + Radiotherapy, x=T, y=T,surv=T, data=uu,time.inc=12)
cal1 <- calibrate(f1, cmethod="KM", method="boot", u=12, m=100, B=1000)
# u=TIME.INC m = 样本量的三分之一到四分之一 b表示最大再抽样的样本量
{
  plot(cal1, lwd=2,cex.axis=1.5,tcl=-0.8, lty=1, errbar.col=c(rgb(0,118,192,maxColorValue=255)), xlab=" ",ylab = "",col=c(rgb(192,98,83,maxColorValue=255)))
  lines(cal1[,c("mean.predicted","KM")], type="b", lwd=3, col=c(rgb(192,98,83,maxColorValue=255)), pch=16)
  mtext("Nomogram predicted 1-year PFS ", side = 1, line = 2.5,cex=1.7)
  mtext("Actual 1-year PFS (proportion)  ", side = 2, line = 2.5,cex=1.7)
  box(lwd = 1)
  abline(0, 1, lty=3, lwd=2, col=c(rgb(0,118,192,maxColorValue=255)))
}
graph2ppt(file="output/plots/nomo.pptx", append = T)
# 3 年
f2<-cph( y ~ group + Age_median + Tstage + NNstage + Radiotherapy, x=T, y=T,surv=T, data=uu,time.inc=36)
cal2<-calibrate(f2, cmethod="KM", method="boot", u=36, m=100, B=1000)
{
  plot(cal2, lwd=2, cex.axis=1.5,tcl=-0.8,lty=1, errbar.col=c(rgb(0,118,192,maxColorValue=255)), xlab=" ",ylab = "",col=c(rgb(192,98,83,maxColorValue=255)))
  lines(cal2[,c("mean.predicted","KM")], type="b", lwd=3, col=c(rgb(192,98,83,maxColorValue=255)), pch=16)
  mtext("Nomogram predicted 3-year PFS ", side = 1, line = 2.5,cex=1.7)
  mtext("Actual 3-year PFS (proportion)  ", side = 2, line = 2.5,cex=1.7)
  box(lwd = 1)
  abline(0, 1, lty=3, lwd=2, col=c(rgb(0,118,192,maxColorValue=255)))
}
graph2ppt(file="output/plots/nomo.pptx", append = T)

# 5年  5年生存率做不了，没有这么多数据
f3 <- cph( y ~ group + Age_median + Tstage + NNstage + Radiotherapy, x=T, y=T,surv=T, data=uu,time.inc=60)
cal3 <- calibrate(f3, cmethod="KM", method="boot", u=60, m=100, B=1000)
{
plot(cal3, lwd=2,cex.axis=1.5,tcl=-0.8, lty=1, errbar.col=c(rgb(0,118,192,maxColorValue=255)), xlab=" ",ylab = "",col=c(rgb(192,98,83,maxColorValue=255)))
lines(cal3[,c("mean.predicted","KM")], type="b", lwd=3, col=c(rgb(192,98,83,maxColorValue=255)), pch=16)
mtext("Nomogram predicted 5-year PFS ", side = 1, line = 2.5,cex=1.7)
mtext("Actual 5-year PFS (proportion)  ", side = 2, line = 2.5,cex=1.7)
box(lwd = 1)
abline(0, 1, lty=3, lwd=2, col=c(rgb(0,118,192,maxColorValue=255)))
}
graph2ppt(file="output/plots/nomo.pptx", append = T)

################ 用data2做nomogram图###########
#需要加一个 R包
dt2 <- data2[-ncol(data2)]
dt2 <- na.omit(dt2)
uu <- dt2
y<-Surv(uu$OS.time,uu$OS==1)
kmfit<-survfit(y~1, data=y)
coxmodel <- coxph(y ~ group+Age+Age_median+Age_group+Sex+
                    Race+Smoke+Tstage+Nstage+NNstage+Mstage+
                    Radiotherapy+Tumor_site, data=uu)
step(coxmodel) ### 这一步是挑选合适的因素进行cox多因素分析，但其实我在这之前已经挑选过了
coxmodel <- coxph(y ~ group + Age_group + Smoke + Tstage + NNstage + 
                    Radiotherapy, data=uu)
library(rms)
dd<-datadist(uu)
options(datadist="dd")

f <- cph( y ~ group + Age_group + Smoke + Tstage + NNstage + 
            Radiotherapy, x=T, y=T,surv=T, data=uu)
surv<-Survival(f)
nom2 <- nomogram(f, fun=list(function(x) surv(12, x), function(x) surv(36, x), function(x) surv(60, x)), lp=F,funlabel=c("1-year survival ", "3-year survival ", "5-year survival "), maxscale=10, fun.at=c(0.95, 0.9, 0.8, 0.7, 0.6, 0.5,0.4,0.3,0.2,0.1,0.05))
plot(nom2,cex.axis = 1.3,cex.var =1.5,xfrac=.4,lmgp=0.5,tcl=-0.5)
graph2ppt(file="output/plots/nomo.pptx",append = T, width=10,height=10)

# 内部验证 C—index = 1-C==1-0.354=0.646
validate(f, method="boot", B=1000, dxy=T)
rcorrcens(Surv(OS.time, OS) ~ predict(f), data = uu) #C—index = 1-C==1-0.327=0.673
Cindex_os_nomo <- summary(coxmodel) 
Cindex_os_nomo <- Cindex_os_nomo$concordance[1] ### C-index 也可以这样计算
Cindex_os_nomo #        C 0.6728721 
# calibration_plot --------------------------------------------------------
# 一致性检验
for (i in c(12,36,60)) {
  print(i)
  f1 <- cph( y ~ T.Stage + N.Stage + M.Stage + RiskType, x=T, y=T,surv=T, data=uu,time.inc=i)
  cal1 <- calibrate(f1, cmethod="KM", method="boot", u=i, m=100, B=1000)
  {
    plot(cal1, lwd=2,cex.axis=1,tcl=-0.8, lty=1, errbar.col=c(rgb(0,118,192,maxColorValue=255)), xlab=" ",ylab = "",col=c(rgb(192,98,83,maxColorValue=255)))
    lines(cal1[,c("mean.predicted","KM")], type="b", lwd=3, col=c(rgb(192,98,83,maxColorValue=255)), pch=16)
    mtext(paste0("Nomogram predicted ",i/12,"-year OS "), side = 1, line = 2,cex=1.2)
    mtext(paste0("Actual ",i/12,"-year OS (proportion)"), side = 2, line = 2,cex=1.2)
    box(lwd = 1)
    abline(0, 1, lty=3, lwd=2, col=c(rgb(0,118,192,maxColorValue=255)))
  } 
  topptx(file="output/11森林图.pptx",append = T, width=6,height=6)
}




# 
# beta <- coef(ppcox)
# se <- sqrt(diag(vcov(ppcox)))
# HR <- exp(beta)
# HRse <- HR * se
# summary(ppcox)

################## 做nomogram图#################
library(regplot)
library(survival)
library(survminer)

obs <- data2[5,-22]
obs
nomo <- regplot(os_multicox, observation=obs, plots = c("density","boxes"),failtime =c(12,36,60), prfail = TRUE,
                boxcol="#ADD8E6",cexvars=1.2,cexscales=1.2,cexcats=1.0,droplines=TRUE, points = T,
                title = "Nomogram to predict OS in LUAD patients")
#res.cox表示模型,可以是广义线性模型（glm）,线性（lm）,生存分析（cox比例风险模型）
#observation指定某个患者各协变量的取值映射到相应的得分，并计算总得分
#failtime = c(12,36,60)计算其在1\3\5年的的累计事件发生概率
#本案例中lung是个生存类数据，status=1代表构建的是生存模型，因此prfail = T

graph2ppt(file="output/plots/nomo.pptx",width=12, aspectr = 1.5,append = T)


nomo

