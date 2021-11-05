

#####################生成基因名称表达谱
dat_exp<-read.csv("TCGA-PAAD.htseq_fpkm-uq.tsv.gz",header = T,sep = "\t",row.names = 1)
colnames(dat_exp)<-gsub(colnames(dat_exp), pattern = '[.]', replacement = '-')
gene_anno<-read.csv("gene_annoation.txt",sep = "\t",header = T)
table(row.names(dat_exp)%in%gene_anno$id)
gene_anno<-gene_anno[match(row.names(dat_exp),gene_anno$id),]
table(row.names(dat_exp)==gene_anno$id)
dat_exp$symbol<-gene_anno$gene
dat_exp<-aggregate(dat_exp[,1:182],by=list(dat_exp$symbol),FUN = mean)
row.names(dat_exp)<-dat_exp$Group.1
dat_exp<-dat_exp[,-1]
tumor_exp<-dat_exp[,-grep("-11",colnames(dat_exp))]
save(tumor_exp,file = "PAAD_symbol_exp.rdata")

###############
# library(reshape2)
# anno_dat<-read.csv("TCGA-PAAD.GDC_phenotype.tsv.gz",header = T,sep = "\t")
# anno_dat<-anno_dat[-grep("-11",anno_dat$submitter_id.samples),]
# # anno_dat<-anno_dat[!duplicated(anno_dat),]
# anno_dat<-anno_dat[anno_dat$tumor_stage.diagnoses!="" & anno_dat$tumor_stage.diagnoses!="not reported",]
# anno_dat$tumor_stage.diagnoses[grep("stage iv",anno_dat$tumor_stage.diagnoses)]<-"Stage IV"
# anno_dat$tumor_stage.diagnoses[grep("stage iii",anno_dat$tumor_stage.diagnoses)]<-"Stage III"
# anno_dat$tumor_stage.diagnoses[grep("stage ii",anno_dat$tumor_stage.diagnoses)]<-"Stage II"
# anno_dat$tumor_stage.diagnoses[grep("stage i",anno_dat$tumor_stage.diagnoses)]<-"Stage I"
# 
# sample_name<-intersect(colnames(tumor_exp),anno_dat$submitter_id.samples)
# anno_dat2<-anno_dat[match(sample_name,anno_dat$submitter_id.samples),]
# tumor_exp2<-tumor_exp[,sample_name]
# table(colnames(tumor_exp2)==anno_dat2$submitter_id.samples)
# colnames(tumor_exp2)<-anno_dat2$tumor_stage.diagnoses
# all_pvalue<-vector()
# for (k in 1:dim(tumor_exp2)[1]) {
#   plot_dat2<-melt(as.matrix(tumor_exp2[k,]))
#   fit <- aov(value ~ Var2, data=plot_dat2)
#   pvalue<-summary(fit)[[1]][1,5]
#   # pvalue<-signif(pvalue,3)
#   # if(pvalue<0.01){
#   #   pdf(paste("Stage_aov_",row.names(tumor_exp2)[k],".pdf",sep = ""),width = 3.8,height = 4)
#   #   print(ggplot(data = plot_dat2, mapping = aes(x = Var2, y =value , color = Var2)) +
#   #           geom_boxplot()+
#   #           stat_boxplot(geom = 'errorbar',width=0.15)+
#   #           geom_jitter(width = 0.2, alpha = 0.5, size=1)+
#   #           theme_bw()+
#   #           # geom_signif(comparisons = list(c("NO","YES")),step_increase=0.1,map_signif_level = F,test = wilcox.test,color="black")+
#   #           theme(legend.position="none")+
#   #           scale_color_manual(values = pal_d3("category10")(10)[1:4])+
#   #           xlab("")+ylab(paste("The expression of ",row.names(tumor_exp2)[k],sep = ""))+
#   #           # labs(title =colnames(plot_dat2)[i])+
#   #           annotate("text",x=2,y=max(plot_dat2$value),label = paste("one-way anova, p =",pvalue)))
#   #   dev.off()
#   # }
#   all_pvalue<-c(all_pvalue,pvalue)
# }
# diff_dat<-data.frame(symbol=row.names(tumor_exp2),pvalue=all_pvalue)
# save(diff_dat,file = "diff_stage_exp.rdata")
# diff_dat<-diff_dat[which(diff_dat$pvalue<0.05),]
###########################hallmarker分析
library(GSVA)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(ggsignif)
library(RColorBrewer)
source("utils.R")
# diff_exp<-tumor_exp[diff_dat$symbol,]
dd<-apply(tumor_exp,1,function(x){sum(x>0)})
diff_exp<-tumor_exp[dd>10,]
pathways <- gmtPathways("h.all.v7.2.symbols.gmt")
hall_gene<-unique(unlist(pathways))
dat_gsva<-diff_exp[row.names(diff_exp)%in%hall_gene,]
es <- gsva(as.matrix(dat_gsva), pathways,
           min.sz=10, max.sz=500, verbose=TRUE)
library(ConsensusClusterPlus)
# drug_acti = sweep(drug_acti,1, apply(drug_acti,1,median,na.rm=T))
results = ConsensusClusterPlus(as.matrix(es),maxK=10,reps=200,pItem=0.8,pFeature=1,
                               title = "aa",clusterAlg="hc",distance="pearson",seed=10,plot="pdf")
consen_class<-results[[4]][["consensusClass"]]

consen_cor<-results[[4]][["ml"]]
colnames(consen_cor)<-names(consen_class)
row.names(consen_cor)<-names(consen_class)
name2<-names(sort(consen_class,decreasing = F))
consen_cor<-consen_cor[name2,name2]

pheatmap(consen_cor,
         # scale = "row",
         cluster_rows = F,
         cluster_cols = F,
         show_colnames = FALSE,
         show_rownames = F,
         color = colorRampPalette(c("white",brewer.pal(9,"BrBG")[5:7]))(100),
         annotation_col = anno_col,
         # annotation_row = anno_row,
         annotation_colors = anno_color,
         border_color = NA
)


##########################
library(wesanderson) 
library(reshape2)
hall_score<-es
table(colnames(hall_score)==names(consen_class))
colnames(hall_score)<-paste("cluster ",as.character(consen_class),sep = "")
hall_score<-as.data.frame(hall_score)
hall_pvalue<-vector()
all_p<-list()
hall_name<-vector()
# col<-RColorBrewer::brewer.pal(n = 9, name = "PiYG")
for (i in 1:dim(hall_score)[1]) {
  plot_dat2<-melt(as.matrix(hall_score[i,]))
  fit <- aov(value ~ Var2, data=plot_dat2)
  pvalue<-summary(fit)[[1]][1,5]
  pvalue<-signif(pvalue,3)
  # hall_pvalue<-c()
  if(pvalue<0.01){
    p<-ggplot(data = plot_dat2, mapping = aes(x = Var2, y =value , color = Var2)) +
                geom_boxplot()+
                stat_boxplot(geom = 'errorbar',width=0.15)+
                geom_jitter(width = 0.2, alpha = 0.5, size=1)+
                theme_bw()+
                theme(legend.position="none")+
                scale_color_manual(values = mycol[1:4])+
                xlab("")+ylab("The activity score")+
                ggtitle(row.names(hall_score)[i])+
                # labs(title =colnames(plot_dat2)[i])+
                annotate("text",x=2,y=max(plot_dat2$value),label = paste("one-way anova, p =",pvalue))+
      coord_flip()
    all_p<-c(all_p,list(p))
    hall_name<-c(hall_name,row.names(hall_score)[i])
  }
}
names(all_p)<-hall_name
save(all_p,file = "hallmark_boxplot.rdata")
library(cowplot)
plot_grid( all_p[["HALLMARK_HYPOXIA"]],all_p[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]]
           ,ncol = 1)
##########################
library(pheatmap)
all_sample_anno<-data.frame(sample=colnames(es),cluster=as.vector(consen_class))
all_sample_anno$type<-"Un-detected"
all_sample_anno$type[all_sample_anno$sample%in%anno_dat$submitter_id.samples[anno_dat$tumor_stage.diagnoses=="Stage I"]]<-"Stage I"
all_sample_anno$type[all_sample_anno$sample%in%anno_dat$submitter_id.samples[anno_dat$tumor_stage.diagnoses=="Stage II"]]<-"Stage II"
all_sample_anno$type[all_sample_anno$sample%in%anno_dat$submitter_id.samples[anno_dat$tumor_stage.diagnoses=="Stage III"]]<-"Stage III"
all_sample_anno$type[all_sample_anno$sample%in%anno_dat$submitter_id.samples[anno_dat$tumor_stage.diagnoses=="Stage IV"]]<-"Stage IV"
all_sample_anno<-all_sample_anno[order(all_sample_anno$cluster),]
save(es,all_sample_anno,file = "GSVA_处理后结果.rdata")

mycol<-pal_d3("category20")(20)
anno_col<-data.frame(Stage=all_sample_anno$type,
                     Cluster=as.character(all_sample_anno$cluster))
row.names(anno_col)<-all_sample_anno$sample
anno_color<-list(Stage=mycol[6:10],
                 Cluster=mycol[1:4])
names(anno_color[[2]])<-unique(as.character(all_sample_anno$cluster))
names(anno_color[[1]])<-unique(as.character(all_sample_anno$type))
# names(anno_color[[3]])<-unique(as.character(all_sample_anno$cluster))

hall_score<-es[,all_sample_anno$sample]
pheatmap(hall_score,
         # scale = "row",
         cluster_rows = T,
         cluster_cols = F,
         show_colnames = FALSE,
         show_rownames = T,
         color = c(colorRampPalette(rev(brewer.pal(9,"BrBG")))(100)),
         annotation_col = anno_col,
         # annotation_row = anno_row,
         annotation_colors = anno_color,
         border_color = NA,
         cutree_cols = 5
)
