rm(list = ls())
###################  整理phenotype的数据 还原OS.time的计算过程  ###############
p_load(ggplot2, ggpubr, dplyr, survival, survminer, stringr,export)
load(file = "outdata/step4_data.rds")
######### 导入xena的phenotype数据

phenotype_xena <- data.table::fread("../maf分组/data/luad_phtnotype_gdc_xena.csv",check.names = F)
colnames(phenotype_xena)[1] <- "sample"

# 在phenotype中选择有fpkm数据的
phenotype_xena[1:4,1:4]
# 有没有缺失值
phenotype_xena1 <- phenotype_xena
table(is.na(phenotype_xena1))# 暂时先不去掉
#
# #####提取status生存状态 death为1（出现结局事件为1），生存（删失数据）为0

phenotype_xena1$OS[phenotype_xena1$vital_status.demographic=='Dead'] <- 1
phenotype_xena1$OS[phenotype_xena1$vital_status.demographic=='Alive'] <- 0
table(is.na(phenotype_xena1$OS)) # 检查有无缺失值

# ######提取time生存时间

phenotype_xena1$days_to_death.demographic[which(phenotype_xena1$OS == 0 )]=0 #将活人的daystodeath定义为0
phenotype_xena1$days_to_last_follow_up.diagnoses[which(phenotype_xena1$OS == 1 )]=0 #将死人的daystofollowup定义为0
phenotype_xena1$OS.time <- phenotype_xena1$days_to_death.demographic + phenotype_xena1$days_to_last_follow_up.diagnoses

# ###生存时间除以30得到月份数
phenotype_xena1$OS.time <- phenotype_xena1$OS.time / 30
table(is.na(phenotype_xena1$OS.time))

# #######加一列mut代表分组信息 

phenotype_xena2 <- phenotype_xena1[ , c('sample', 'OS','OS.time')] #取OS和OS.time以及分组做KM分析
phenotype_xena2$group <- 'A'
phenotype_xena2$group[match(x =  (allid[1:279]), table = phenotype_xena2$sample)] <- 'B'

# setequal(phenotype_xena2$OS.time, cln_xena1$OS.time) # 查看phenotype和survival数据得到的OS.time是否相同
# 注意这里其实就算false也无所谓的，因为clnxena1是重新读进来的，有些无限小数保留成15位有效数字的大小了
# setequal(cln_xena1$OS.time, signif(x=phenotype_xena2$OS.time, digits = 15))


# ### 检查group分组信息情况

# ###group必须是因子，代表有无突变
str(phenotype_xena2)
phenotype_xena2$group <- as.factor(phenotype_xena2$group)
phenotype_xena2$group |> table()
phenotype_xena2$OS <- as.integer(phenotype_xena2$OS)
boxplot(phenotype_xena2$OS.time)
library(dplyr)

data1 <- phenotype_xena2 %>% filter(OS.time<=100)

fit<-survfit(Surv(OS.time,OS)~group, data = data1) 
pp_os <- ggsurvplot(fit, data=data1, linetype = 1, palette = c("jco"),size=1,
                    surv.scale = c("percent"),pval = TRUE,
                    legend.title = '', legend.labs=c("Mut","Wild"),
                    break.time.by =12,
                    xlim = c(0,60),risk.table = TRUE,
                    risk.table.title = "Patients at risk",
                    ylab = "Overall Survival, %",xlab = "Months",
                    font.x = c(20,"plain","black"),font.y = c(20,"plain","black"),
                    font.tickslab = c(20,"plain","black"),
                    risk.table.fontsize = 6,font.legend =c(20,"plain","black"),
                    font.main = c(20,"plain","black"),pval.size = 8)

######risk table 的修改
pp_os$table <- ggpar(
  pp_os$table,
  font.title    = c(13, "bold.italic", "green"),  ###risk table 标题的修改
  font.subtitle = c(15, "bold", "pink"),  ###risk table小标题的修改
  font.caption  = c(11, "plain", "darkgreen"), ####插入字的修改
  font.x        = c(18, "plain", "black"), ### risk table x的修改
  font.y        = c(18, "plain", "black"),### risk table y的修改
  font.xtickslab = c(18, "plain", "black"),### risk table x 坐标轴的修改
  font.ytickslab = c(18),  ### risk table y 坐标轴的修改
  legend=c(0.8,0.88),
  censor.size=3
)
print(pp_os) ##终于做出和survival一样的生存分析图了
graph2ppt(file = "output/plots/kmplot_os_by_xena_phenotype")




################################## 无进展生存期PFS生存曲线###############
# PFS 无进展生存期 LUAD数据库中的new_tumor_event_after_initial_treatment
# 把新生肿瘤时间的NA值设置为9999
phenotype_xena2$PFS.time <- phenotype_xena1$days_to_new_tumor_event_after_initial_treatment / 30
phenotype_xena2$PFS.time[which(is.na(phenotype_xena2$PFS.time))] <- 9999
phenotype_xena2$PFS.time <- ifelse(phenotype_xena2$PFS.time < phenotype_xena2$OS.time,  phenotype_xena2$PFS.time,  phenotype_xena2$OS.time)
phenotype_xena2$PFS <- phenotype_xena2$OS
phenotype_xena2$PFS[!is.na(phenotype_xena1$days_to_new_tumor_event_after_initial_treatment)] <- 1
table(phenotype_xena2$PFS)
#####无进展生存期数据挖掘完成

##########pfs 的生存曲线
table(is.na(phenotype_xena2$PFS)) # 检查有无缺失值

pfsdata <- phenotype_xena2[,c("sample","group", "OS", "OS.time", "PFS", "PFS.time" )]
str(pfsdata)
pfsdata$group <- as.factor(pfsdata$group)
pfsdata$PFS <- as.integer(pfsdata$PFS)
boxplot(pfsdata$PFS.time)
library(dplyr)

data1 <- pfsdata %>% filter(PFS.time<=62)

fit<-survfit(Surv(PFS.time,PFS)~group, data = data1) 
pp_pfs <- ggsurvplot(fit, data=data1, linetype = 1, palette = c("jco"),size=1,
                     surv.scale = c("percent"),pval = TRUE,
                     legend.title = '', legend.labs=c("Mut","Wild"),
                     break.time.by =12,
                     xlim = c(0,60),risk.table = TRUE,
                     risk.table.title = "Patients at risk",
                     ylab = "Progression Free Survival, %",xlab = "Months",
                     font.x = c(20,"plain","black"),font.y = c(20,"plain","black"),
                     font.tickslab = c(20,"plain","black"),
                     risk.table.fontsize = 6,font.legend =c(20,"plain","black"),
                     font.main = c(20,"plain","black"),pval.size = 8)

######risk table 的修改
pp_pfs$table <- ggpar(
  pp_pfs$table,
  font.title    = c(13, "bold.italic", "green"),  ###risk table 标题的修改
  font.subtitle = c(15, "bold", "pink"),  ###risk table小标题的修改
  font.caption  = c(11, "plain", "darkgreen"), ####插入字的修改
  font.x        = c(18, "plain", "black"), ### risk table x的修改
  font.y        = c(18, "plain", "black"),### risk table y的修改
  font.xtickslab = c(18, "plain", "black"),### risk table x 坐标轴的修改
  font.ytickslab = c(18),  ### risk table y 坐标轴的修改
  legend=c(0.8,0.88),
  censor.size=3
)
print(pp_pfs) 
surv_summary(fit) #查看生存率及其标准误
surv_median(fit,combine = F) # 查看中位生存时间


phenotype_xena1$PFS <- phenotype_xena2$PFS
phenotype_xena1$PFS.time <- phenotype_xena2$PFS.time

#######################准备数据phenotype数据制作tableone生存资料基线表#############
#sex,age,race,smoke,stage,T,N,M,Anatomic location
#sex(gender)，age(选择病理确诊年龄age_at_initial_pathologic_diagnosis)，
# race,
# Lifelong Non-smoker (less than 100 cigarettes smoked in Lifetime) = 1:从未吸烟
# Current smoker (includes daily smokers and non-daily smokers or occasional smokers) = 2：当前吸烟（包括规律吸烟者和偶尔吸烟者）
# Current reformed smoker for > 15 years (greater than 15 years) = 3：戒烟>15年
# Current reformed smoker for ≤15 years (less than or equal to 15 years) = 4：戒烟<15年
# Current reformed smoker, duration not specified = 5：戒烟，戒烟时间不详
# Smoker at Diagnosis = 6：确诊时吸烟
# Smoking History not documented = 7：不详

#smoke(tobacco_smoking_history)，stage(tumor_stage.diagnoses),grade,
# T,N,M(pathologic_T),lung(site_of_resection_or_biopsy.diagnoses),treatment.

phenotype_xena1 <- phenotype_xena1 |> mutate(PFS=phenotype_xena2$PFS,
                                             PFS.time = phenotype_xena2$PFS.time)

phenotype_xena3 <- phenotype_xena1[,c('sample', 'OS','OS.time','PFS','PFS.time',  
                                       'gender.demographic','age_at_initial_pathologic_diagnosis',
                                       'race.demographic', 'tobacco_smoking_history', 'tumor_stage.diagnoses',
                                       'site_of_resection_or_biopsy.diagnoses','pathologic_T',
                                       'pathologic_N', 'pathologic_M', 'radiation_therapy')]

colnames(phenotype_xena3)[which(names(phenotype_xena3) == "gender.demographic"):ncol(phenotype_xena3)] <- c('gender', 'age', 'race', 'smoke', 'stage',
                                                                                                            'site', 'ajcc_T', 'ajcc_N', 'ajcc_M','radiotherapy' ) 
is.na(phenotype_xena3) |> table()
phenotype_xena3 <- phenotype_xena3 |> filter(!str_detect(sample,"TCGA-17-Z0",))
# sex
table(is.na(phenotype_xena3$gender))
table(phenotype_xena3$gender,useNA="ifany")


# age
phenotype_xena3$age
which(is.na(phenotype_xena3$age))
table(is.na(phenotype_xena3$age)) # 10个病人缺失确诊年龄信息 要体现在基线表中
#age_medi 根据中位年龄分组
median_age <- median(phenotype_xena3$age, na.rm = T) #中位年龄为66
phenotype_xena3$age_median <- ifelse( phenotype_xena3$age >= median_age, 'older', 'younger' )
table(phenotype_xena3$age_median,useNA = "ifany")
# race race/ethnicity
which(is.na(phenotype_xena3$race))
table(phenotype_xena3$race, useNA = "ifany")
phenotype_xena3$race[phenotype_xena3$race=="black or african american"] <- "black"
phenotype_xena3$race[phenotype_xena3$race!='white' &
                       phenotype_xena3$race != "black"] <- "other"
table(phenotype_xena3$race)


# stage
table(phenotype_xena3$stage,useNA = "ifany")
phenotype_xena3$stage[phenotype_xena3$stage=='not reported'] <- NA
#phenotype_xena3$stage <- str_split(string = phenotype_xena3$stage, pattern = ' ', n = 2, simplify = T)[,2]
#n代表返回几个字符串 simplify为T表示返回字符串
phenotype_xena3$stage1 <- substring(text = phenotype_xena3$stage, first = 7)
table(phenotype_xena3$stage1)

# 把stage分为早期(1,2)和晚期(3,4)
# library(stringr)
# phenotype_xena3$stage2 <- str_split(string = phenotype_xena3$stage1,
#                                     pattern = 'a', n = 2, simplify = T)[,1]
# phenotype_xena3$stage2 <- str_split(string = phenotype_xena3$stage2,
#                                     pattern = 'b', n = 2, simplify = T)[,1]

phenotype_xena3$stage2 <- unlist(strsplit(phenotype_xena3$stage1, "a"))
phenotype_xena3$stage2 <- unlist(strsplit(phenotype_xena3$stage2, "b"))

phenotype_xena3$stage3[phenotype_xena3$stage2=='i' |  
                         phenotype_xena3$stage2=='ii' ] <- "early stage"  #1,2为早期
phenotype_xena3$stage3[phenotype_xena3$stage2=='iii' |  
                         phenotype_xena3$stage2=='iv' ] <- "later stage"  #3,4为晚期

phenotype_xena3 <- dplyr::rename(phenotype_xena3, stage_group= stage3) #改名

table(phenotype_xena3$stage2, useNA = "if")
table(phenotype_xena3$stage_group, useNA = "if")

# fix(phenotype_xena3) 
# 还有很多方法都可以。这种方法应该是比较简单的，一行代码搞定。


# smoke smokesta 分成吸烟和不吸烟
phenotype_xena3$smoke
table(phenotype_xena3$smoke, useNA = "ifany")
phenotype_xena3$smoke_group <- ifelse(phenotype_xena3$smoke == 1|
                                        phenotype_xena3$smoke == 3, "NO", "YES")
table(phenotype_xena3$smoke_group, useNA = "ifany")


# T,N,M T分期不进行改动
phenotype_xena3$ajcc_T <- substring(text = phenotype_xena3$ajcc_T, first = 1,last = 2 )
phenotype_xena3$ajcc_T[which(phenotype_xena3$ajcc_T=="" | 
                               phenotype_xena3$ajcc_T=="TX")] <- NA
table(phenotype_xena3$ajcc_T, useNA = "ifany") 

phenotype_xena3$ajcc_N <- substring(text = phenotype_xena3$ajcc_N, first = 1,last = 2 )
phenotype_xena3$ajcc_N[which(phenotype_xena3$ajcc_N==""| 
                               phenotype_xena3$ajcc_N=="NX")] <- NA
table(phenotype_xena3$ajcc_N, useNA = "ifany")

phenotype_xena3$ajcc_M <- substring(text = phenotype_xena3$ajcc_M, first = 1,last = 2 )
table(phenotype_xena3$ajcc_M, useNA = "ifany")
phenotype_xena3$ajcc_M[which(phenotype_xena3$ajcc_M=="")] <- NA
table(phenotype_xena3$ajcc_M, useNA = "ifany")


#######resection site
table(phenotype_xena3$site, useNA = "ifany")
phenotype_xena3$site2 <- "Other site"
phenotype_xena3$site2[which(phenotype_xena3$site == "Lower lobe, lung")] <- "Lower lobe"
phenotype_xena3$site2[which(phenotype_xena3$site == "Upper lobe, lung")] <- "Upper lobe"
phenotype_xena3$site2[which(phenotype_xena3$site == "Middle lobe, lung")] <- "Middle lobe"
phenotype_xena3 <- dplyr::rename(phenotype_xena3, resection_site = site2) #改名

###### radiation_therapy 
phenotype_xena3$radiotherapy 
table(phenotype_xena3$radiotherapy, useNA = "ifany")
phenotype_xena3$radiotherapy[which(phenotype_xena3$radiotherapy=="")] <- NA 
table(phenotype_xena3$radiotherapy, useNA = "ifany")


phe_f <- phenotype_xena3 %>% 
  select(c("sample", "OS", "OS.time","PFS", "PFS.time",  "gender", 
           "age", "age_median","race",  "smoke","smoke_group",
           "ajcc_T", "ajcc_N", "ajcc_M", "stage_group", "stage2",
           "resection_site","radiotherapy"))


phe_f$age_group2[phe_f$age >= 60 &
                   phe_f$age < 70  ] <- "60-70"
phe_f$age_group2[phe_f$age > 70  ] <- ">70"
phe_f$age_group2[phe_f$age < 60  ] <- "<60"
phe_f$age_group2 |> table()
save(phe_f,phenotype_xena,file = "outdata/LUAD_cln_clean.rds")
