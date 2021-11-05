library(pROC) #绘制ROC曲线
library(e1071)
library(RColorBrewer)
library(pheatmap)

load("hall_exp.rdata")
load("boruta_feature.rdata")
par(mfrow=c(1,2))
col<-brewer.pal(11,"Spectral")
confirm_result<-gsub("`","",confirm_result)
feature_exp<-hall_exp[,confirm_result]
type<-rep("1",dim(hall_exp)[1])
type[hall_exp$class%in%c(1,3)]<-"2"
feature_exp$class<-type

# train_sub = sample(nrow(feature_exp),0.6*nrow(feature_exp))
# # save(train_sub,file = "训练集样本的位置.rdata")
# train_data = feature_exp[train_sub,]
# test_data = feature_exp[-train_sub,]

data <- feature_exp

library(plyr)
CVgroup <- function(k,datasize,seed){
  cvlist <- list()
  set.seed(seed)
  n <- rep(1:k,ceiling(datasize/k))[1:datasize]    #将数据分成K份，并生成的完成数据集n
  temp <- sample(n,datasize)   #把n打乱
  x <- 1:k
  dataseq <- 1:datasize
  cvlist <- lapply(x,function(x) dataseq[temp==x])  #dataseq中随机生成k个随机有序数据列
  return(cvlist)
}
k <- 5
datasize <- nrow(feature_exp)
cvlist <- CVgroup(k = k,datasize = datasize,seed = 100)
cvlist
save(cvlist,file = "5重交叉验证位置.rdata")
#######################SVM
k=5
pred <- data.frame()
for (i in 1:k){
  train <- data[-cvlist[[i]],]  #刚才通过cvgroup生成的函数
  test <- data[cvlist[[i]],]
  model <-svm(as.numeric(class) ~ ., 
              data = train,
              type = 'eps',kernel = 'radial')  #建模，ntree=j 指的树数
  prediction <- predict(model,subset(test,select = -class))   #预测
  foroc<-roc(test$class,as.numeric(prediction))
  # randomtree <- rep(j,length(foroc[["sensitivities"]]))   #随机森林树的数量
  kcross <- rep(i,length(foroc[["sensitivities"]]))   #i是第几次循环交叉，共K次
  temp <- data.frame(kcross,
                     sensitivities=foroc[["sensitivities"]],
                     specificities=foroc[["specificities"]],
                     auc=rep(foroc[["auc"]][1],length(foroc[["sensitivities"]])))#真实值、预测值、随机森林树数、预测组编号捆绑在一起组成新的数据框tenp
  pred <- rbind(pred,temp)   #temp按行和pred合并
  # print(paste("SVM：",j))  #循环至树数j的随机森林模型
  # progress.bar$step() #输出进度条。告知完成了这个任务的百分之几
}
library(ggplot2)
library(RColorBrewer)
mycol<-c("#2E9FDF", "#4C931B", "#E7B800", "#FC4E07","#DD3497")
p1<-ggplot(pred,aes(x = 1-specificities, y = sensitivities, color = as.factor(kcross))) + 
  geom_path(size=1) + 
  theme_bw()+
  labs(x = "1-specificity", y = "sensitivity", title ="ROC curve of SVM") +
  scale_colour_manual(values = mycol)+ ###############"YlGn","Reds","RdPu"
  annotate("text",x=0.7,y=0.50,label=paste("K=1, AUC=",round(pred$auc[pred$kcross==1][1],3),sep = ""),size=4,color=mycol[1])+
  annotate("text",x=0.7,y=0.40,label=paste("K=2, AUC=",round(pred$auc[pred$kcross==2][1],3),sep = ""),size=4,color=mycol[2])+
  annotate("text",x=0.7,y=0.30,label=paste("K=3, AUC=",round(pred$auc[pred$kcross==3][1],3),sep = ""),size=4,color=mycol[3])+
  annotate("text",x=0.7,y=0.20,label=paste("K=4, AUC=",round(pred$auc[pred$kcross==4][1],3),sep = ""),size=4,color=mycol[4])+
  annotate("text",x=0.7,y=0.10,label=paste("K=5, AUC=",round(pred$auc[pred$kcross==5][1],3),sep = ""),size=4,color=mycol[5])





######################决策树分类器
library("rpart")
library("rpart.plot")
library(survival)
pred <- data.frame()
for (i in 1:k){
  train <- data[-cvlist[[i]],]  #刚才通过cvgroup生成的函数
  test <- data[cvlist[[i]],]
  model <-rpart(Surv(as.numeric(class))~.,data = train, method = 'exp')
  # rpart.plot(model_tree,type = 2)
  prediction <- predict(model,subset(test,select = -class))   #预测
  foroc<-roc(test$class,prediction)
  # randomtree <- rep(j,length(foroc[["sensitivities"]]))   #随机森林树的数量
  kcross <- rep(i,length(foroc[["sensitivities"]]))   #i是第几次循环交叉，共K次
  temp <- data.frame(kcross,
                     sensitivities=foroc[["sensitivities"]],
                     specificities=foroc[["specificities"]],
                     auc=rep(foroc[["auc"]][1],length(foroc[["sensitivities"]])))#真实值、预测值、随机森林树数、预测组编号捆绑在一起组成新的数据框tenp
  pred <- rbind(pred,temp)   #temp按行和pred合并
  # print(paste("SVM：",j))  #循环至树数j的随机森林模型
  # progress.bar$step() #输出进度条。告知完成了这个任务的百分之几
}
library(ggplot2)
library(RColorBrewer)
mycol<-c("#2E9FDF", "#4C931B", "#E7B800", "#FC4E07","#DD3497")
p2<-ggplot(pred,aes(x = 1-specificities, y = sensitivities, color = as.factor(kcross))) + 
  geom_path(size=1) + 
  theme_bw()+
  labs(x = "1-specificity", y = "sensitivity", title ="ROC curve of Decision tree") +
  scale_colour_manual(values = mycol)+ ###############"YlGn","Reds","RdPu"
  annotate("text",x=0.7,y=0.50,label=paste("K=1, AUC=",round(pred$auc[pred$kcross==1][1],3),sep = ""),size=4,color=mycol[1])+
  annotate("text",x=0.7,y=0.40,label=paste("K=2, AUC=",round(pred$auc[pred$kcross==2][1],3),sep = ""),size=4,color=mycol[2])+
  annotate("text",x=0.7,y=0.30,label=paste("K=3, AUC=",round(pred$auc[pred$kcross==3][1],3),sep = ""),size=4,color=mycol[3])+
  annotate("text",x=0.7,y=0.20,label=paste("K=4, AUC=",round(pred$auc[pred$kcross==4][1],3),sep = ""),size=4,color=mycol[4])+
  annotate("text",x=0.7,y=0.10,label=paste("K=5, AUC=",round(pred$auc[pred$kcross==5][1],3),sep = ""),size=4,color=mycol[5])
library(cowplot)
plot_grid(p1,p2)
