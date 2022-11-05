set.seed(1999)
x1 = rnorm(1000)          
x2 = rnorm(1000)
z = 1 + 2*x1 + 3*x2       
pr = 1/(1+exp(-z))        
y = rbinom(1000,1,pr)
#使用logloss作为训练目标函数
df = data.frame(y=y,x1=x1,x2=x2)
glm.fit=glm( y~x1+x2,data=df,family="binomial")

#下面使用auc作为训练目标函数
library(ROCR)

CalAUC <- function(real,pred){
  rocr.pred=prediction(pred,real)
  rocr.perf=performance(rocr.pred,'auc')
  as.numeric(rocr.perf@y.values)
}

#初始值
beta0=c(1,1,1)

loss <- function(beta){
  z=beta[1]+beta[2]*x1+beta[3]*x2
  pred=1/(1+exp(-z))
  -CalAUC(y,pred)
}

res=optim(beta0,loss,method = "Nelder-Mead",control = list(maxit = 100))
cat("直接用AUC训练:",-res$value)
cat("使用glm函数",CalAUC(y,glm.fit$fitted.values))





# ROC_柳叶刀小鼠标---------------------------------------------------------------------
library(ROCR)
data(ROCR.simple)
pred <- prediction(ROCR.simple$predictions, ROCR.simple$labels) 
#ROCR.simple$predictions为预测标签，ROCR.simple$labels为真实标签
perf <- performance(pred,"tpr","fpr")
auc <- performance(pred,'auc')
auc = unlist(slot(auc,"y.values"))
plot(perf,
     xlim=c(0,1), ylim=c(0,1),col='red', 
     main=paste("ROC curve (", "AUC = ",auc,")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)

