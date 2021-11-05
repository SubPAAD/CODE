# # install.packages("factoextra")
# library("factoextra")
# load("immun_risk_cor.rdata")
# all_cor_immun<-all_cor_immun[all_cor_immun[,1]<0.05,]
# load("处理后的cibresort_result.rdata")
# result_ciber<-result_ciber[-which(row.names(result_ciber)=="Plasma cells"),]
# result_ciber<-result_ciber[row.names(all_cor_immun),]
# #先求样本之间两两相似性
# result <- dist(t(result_ciber), method = "euclidean")
# #产生层次结构
# result_hc <- hclust(d = result, method = "ward.D2")
# #进行初步展示
# fviz_dend(result_hc, k = 3, 
#           cex = 0.5, 
#           k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
#           color_labels_by_k = TRUE, 
#           rect = TRUE          
# )
# cluster2<-cutree(result_hc,k=3)
# cluster2<-cluster2[all_sample_anno$sample]
# table(names(cluster2)==all_sample_anno$sample)
# table(cluster2,all_sample_anno$cluster)

################筛选特征
# load("差异表达致癌基因.rdata")
load("转录调控网络.rdata")
load("PAAD_symbol_exp.rdata")
load("GSVA_处理后结果.rdata")
# feature<-names(hall_pvalue[hall_pvalue<0.01])
table(hall_tf_tg$pvalue<0.01&hall_tf_tg$cor>0.6)
hall_tf_tg2<-hall_tf_tg[hall_tf_tg$pvalue<0.01&hall_tf_tg$cor>0.6,]
feature<-unique(c(hall_tf_tg2$TF.Name,hall_tf_tg2$TG.Name))

hall_exp<-as.data.frame(t(tumor_exp[row.names(tumor_exp)%in%feature,]))
hall_exp<-hall_exp[all_sample_anno$sample,]
table(row.names(hall_exp)==all_sample_anno$sample)
hall_exp$class<-all_sample_anno$cluster
save(hall_exp,file = "hall_exp.rdata")
###############Boruta赛选特征
library(Boruta)
gse10_boruta<-Boruta(class~.,data=hall_exp,doTrace = 2)
print(gse10_boruta)
plot(gse10_boruta,xlab = "", xaxt = "n")
lz<-lapply(1:ncol(gse10_boruta$ImpHistory),function(i)
  
  gse10_boruta$ImpHistory[is.finite(gse10_boruta$ImpHistory[,i]),i])
names(lz) <- colnames(gse10_boruta$ImpHistory)
Labels <- sort(sapply(lz,median))
axis(side = 1,las=2,labels = names(Labels),
     at = 1:ncol(gse10_boruta$ImpHistory), cex.axis = 0.5)
#########
final_gse10_boruta <- TentativeRoughFix(gse10_boruta)
print(final_gse10_boruta)
plot(final_gse10_boruta,xlab = "", xaxt = "n")
axis(side = 1,las=2,labels = names(Labels),
     at = 1:ncol(final_gse10_boruta$ImpHistory), cex.axis = 0.5)

confirm_result<-getSelectedAttributes(final_gse10_boruta, withTentative = F)
save(confirm_result,file = "boruta_feature.rdata")
