
library(pheatmap)
library(ConsensusClusterPlus)
library(ggplot2)
source('CIBERSORT.R')
mycol=c("#C6147E","#C08229","#4C931B","#329890")
load("PAAD_symbol_exp.rdata")
lm22<-read.csv("LM22.txt",header = T,sep = "\t")
immun_dat<-tumor_exp[row.names(tumor_exp)%in%lm22$Gene.symbol,]
# write.table(immun_dat,"immun_dat.txt",sep = "\t",row.names = T,col.names = T,quote =F)

ciber_dat <- CIBERSORT('LM22.txt','immun_dat.txt', perm = 1000, QN = TRUE)
save(ciber_dat,file = "cibresort_result.rdata")

ciber_dat<-ciber_dat[,1:22]
result_ciber<-apply(ciber_dat,1,function(x){log(x[3]+0.001)/log(x+0.001)})
result_ciber<-t(scale(t(result_ciber[1:22,]),center = T))
result_ciber[result_ciber=="NaN"]<-0
# save(result_ciber,file = "处理后的cibresort_result.rdata")
result_ciber[result_ciber>2]<-2

table(all_sample_anno$sample%in%colnames(result_ciber))
result_ciber<-result_ciber[,all_sample_anno$sample]
table(all_sample_anno$sample==colnames(result_ciber))
result_ciber<-result_ciber[-which(row.names(result_ciber)=="Plasma cells"),]
result_ciber2<-as.data.frame(result_ciber)
colnames(result_ciber2)<-paste("Cluster",all_sample_anno$cluster,sep = "_")

library(reshape2)
all_p<-list()
past_gene_name<-vector()
for (i in 1:dim(result_ciber2)[1]) {
  dat<-result_ciber2[i,]
  dat2<-melt(as.matrix(dat))
  colnames(dat2)<-c("cell","cluster","value")
  fit <- aov(value ~ cluster, data=dat2)
  pvalue<-summary(fit)[[1]][1,5]
  pvalue<-signif(pvalue,3)
  # pvalue<-formatC(pvalue,format = "e",small.interval=2)
  if (pvalue<0.05){
    p<-ggplot(dat2, aes(x=cluster, y=value,color=cluster))+
      # facet_wrap(~Var1,ncol = 5)+
      stat_boxplot(geom = 'errorbar',width=0.15)+
      geom_boxplot()+
      # stat_summary(fun.data = calc_stat, geom="boxplot")+
      geom_jitter(width = 0.2, alpha = 0.5, size=1)+
      # ggtitle(paste("The fraction of ",row.names(result_ciber2)[i])) +
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 30, hjust = 1)) +
      labs(x="", y=paste("The fraction of ",row.names(result_ciber2)[i]))+
      annotate("text",x=2,y=max(dat2$value),label = paste("one-way anova, p =",pvalue))+
      scale_color_manual(values = mycol)
    all_p<-c(all_p,list(p))
    past_gene_name<-c(past_gene_name,row.names(result_ciber2)[i])
  }
  
}
library(cowplot)
plot_grid(all_p[[1]],all_p[[2]],all_p[[3]],all_p[[4]],all_p[[5]],all_p[[6]],ncol = 3)
save(all_p,file = "免疫细胞分布的boxplot.rdata")
#########################
library(ggsci)
library(RColorBrewer)
mycol<-pal_d3("category20")(20)
anno_col<-data.frame(Stage=all_sample_anno$type,
                     Cluster=as.character(all_sample_anno$cluster))
row.names(anno_col)<-all_sample_anno$sample
anno_color<-list(Stage=mycol[6:10],
                 Cluster=mycol[1:4])
names(anno_color[[2]])<-unique(as.character(all_sample_anno$cluster))
names(anno_color[[1]])<-unique(as.character(all_sample_anno$type))

aa<-result_ciber[,order(result_ciber[13,],decreasing = T)]
pheatmap(aa,
         # scale = "row",
         cluster_rows = T,
         cluster_cols = F,
         show_colnames = FALSE,
         show_rownames = T,
         color = c(colorRampPalette(rev(brewer.pal(9,"BrBG")))(100)),
         annotation_col = anno_col,
         # annotation_row = anno_row,
         annotation_colors = anno_color,
         border_color = NA
)
######################
###################免疫原性分析
library(reshape2)
library(ggsignif)
feature1<-c("HLA-DRA","HLA-DRB1","CD74")
feature2<-c("HLA-E","HLA-F","NLRC5")
feature_exp<-tumor_exp[row.names(tumor_exp)%in%c(feature1,feature2),]
feature_exp<-feature_exp[,all_sample_anno$sample]
table(colnames(feature_exp)==all_sample_anno$sample)
colnames(feature_exp)<-paste("Cluster",all_sample_anno$cluster,sep = "_")
feature1_exp<-feature_exp[feature1,]
feature2_exp<-feature_exp[feature2,]

all_p<-list()
for (i in 1:dim(feature_exp)[1]) {
  dat<-melt(as.matrix(feature_exp[i,]))
  fit <- aov(value ~ Var2, data=dat)
  pvalue<-summary(fit)[[1]][1,5]
  pvalue<-signif(pvalue,3)
  p<-ggplot(data=dat, mapping = aes(x = Var2, y = value, color = Var2)) +geom_violin(trim=FALSE,scale = "width")+
    geom_boxplot(width=0.3,position=position_dodge(0.8))+
    theme_bw()+
    geom_jitter(width = 0.2, alpha = 0.5, size=1)+
    theme(legend.position="none",axis.text.x = element_text(angle = 30, hjust = 1))+
    scale_color_manual(values = mycol)+
    xlab("")+ylab("The expression of gene (FPKM)")+ggtitle(row.names(feature_exp)[i])+
    annotate("text",x=2,y=max(dat$value),label = paste("one-way anova, p =",pvalue))
  all_p<-c(all_p,list(p))
}
library(cowplot)
plot_grid(all_p[[1]],all_p[[2]],all_p[[3]],
          all_p[[4]],all_p[[5]],all_p[[6]],ncol = 3)
