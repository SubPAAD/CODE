library(ggplot2)
library(ggsignif)
library(RColorBrewer)
mycol=c("#C6147E","#C08229","#4C931B","#329890")
load("all_risk_score.rdata")
anno_dat<-read.csv("TCGA-PAAD.GDC_phenotype.tsv.gz",header = T,sep = "\t")
anno_dat<-anno_dat[,c(1,grep("tumor_stage.diagnoses",colnames(anno_dat))
                      ,grep("primary_diagnosis.diagnoses",colnames(anno_dat)),
                      grep("prior_malignancy.diagnoses",colnames(anno_dat)),
                      grep("history_of_neoadjuvant_treatment",colnames(anno_dat)),
                      grep("radiation_therapy",colnames(anno_dat)),
                      grep("tissue_or_organ_of_origin.diagnoses",colnames(anno_dat)),
                      grep("tobacco_smoking_history",colnames(anno_dat)))] 

anno_dat<-anno_dat[-grep("-11",anno_dat$submitter_id.samples),]
# anno_dat<-anno_dat[!duplicated(anno_dat),]
anno_dat<-anno_dat[anno_dat$tumor_stage.diagnoses!="" & anno_dat$tumor_stage.diagnoses!="not reported",]
anno_dat$tumor_stage.diagnoses[grep("stage iv",anno_dat$tumor_stage.diagnoses)]<-"Stage IV"
anno_dat$tumor_stage.diagnoses[grep("stage iii",anno_dat$tumor_stage.diagnoses)]<-"Stage III"
anno_dat$tumor_stage.diagnoses[grep("stage ii",anno_dat$tumor_stage.diagnoses)]<-"Stage II"
anno_dat$tumor_stage.diagnoses[grep("stage i",anno_dat$tumor_stage.diagnoses)]<-"Stage I"

# anno_dat$submitter_id.samples<-unlist(lapply(strsplit(anno_dat$submitter_id.samples,split = "-"),function(x){
#   paste(x[1:3],collapse  = "-")
# }))
# anno_dat<-anno_dat[!duplicated(anno_dat),]
samp_inte<-intersect(anno_dat$submitter_id.samples,all_risk_score$sample)

row.names(anno_dat)<-anno_dat$submitter_id.samples
row.names(all_risk_score)<-all_risk_score$sample
anno_dat<-anno_dat[samp_inte,]
all_risk_score<-all_risk_score[samp_inte,]
table(all_risk_score$sample==anno_dat$submitter_id.samples)

dat<-cbind(anno_dat,all_risk_score$total_risk_score)
colnames(dat)[9]<-"risk_score"

fit <- aov(risk_score ~ tumor_stage.diagnoses, data=dat)
pvalue<-summary(fit)[[1]][1,5]
pvalue<-signif(pvalue,3)
p1<-ggplot(data = dat, mapping = aes(x = dat[,2], y = dat[,9], fill = dat[,2])) +
  geom_boxplot()+
  # geom_jitter(width = 0.2, alpha = 0.5, size=1)+
  theme_bw()+
  # geom_signif(comparisons = list(c("Stage I","Stage II"),
  #                                c("Stage I","Stage III"),
  #                                c("Stage I","Stage IV"),
  #                                c("Stage II","Stage IV")
  # ),step_increase=0.1,map_signif_level = F,test = wilcox.test,color="black")+
  theme(legend.position="none")+
  scale_fill_manual(values = mycol)+
  xlab("")+ylab("Risk score")+
  labs(title =colnames(dat)[2])+
  annotate("text",x=2,y=max(dat$risk_score),label = paste("one-way anova, p =",pvalue))

dat$tissue_or_organ_of_origin.diagnoses<-as.character(dat$tissue_or_organ_of_origin.diagnoses)
a<-table(anno_dat$tissue_or_organ_of_origin.diagnoses)
b<-names(a)[a<10]
dat2<-dat[!dat$tissue_or_organ_of_origin.diagnoses%in%b,]
dat2<-dat2[!is.na(dat2$tissue_or_organ_of_origin.diagnoses),]
fit <- aov(risk_score ~ tissue_or_organ_of_origin.diagnoses, data=dat2)
pvalue<-summary(fit)[[1]][1,5]
pvalue<-signif(pvalue,3)
p2<-ggplot(data = dat2, mapping = aes(x = dat2[,7], y = dat2[,9], fill = dat2[,7])) +
  geom_boxplot()+
  theme_bw()+
  geom_signif(comparisons = list(c("Body of pancreas","Head of pancreas"),
                                 c("Tail of pancreas","Pancreas, NOS"),
                                 c("Body of pancreas","Pancreas, NOS")
  ),step_increase=0.1,map_signif_level = F,test = wilcox.test,color="black")+
  theme(legend.position="none")+
  scale_fill_manual(values = mycol)+
  xlab("")+ylab("Risk score")+
  labs(title =colnames(dat2)[7])+
  annotate("text",x=2,y=max(dat$risk_score),label = paste("one-way anova, p =",pvalue))



dat$primary_diagnosis.diagnoses<-as.character(dat$primary_diagnosis.diagnoses)
a<-table(dat$primary_diagnosis.diagnoses)
b<-names(a)[a<6]
dat3<-dat[!dat$primary_diagnosis.diagnoses%in%b,]
dat3<-dat3[!is.na(dat3$primary_diagnosis.diagnoses),]
fit <- aov(risk_score ~ primary_diagnosis.diagnoses, data=dat3)
pvalue<-summary(fit)[[1]][1,5]
pvalue<-signif(pvalue,3)
p3<-ggplot(data = dat3, mapping = aes(x = dat3[,3], y = dat3[,9], fill = dat3[,3])) +
  geom_boxplot()+
  theme_bw()+
  # geom_signif(comparisons = list(c("Infiltrating duct carcinoma, NOS","Adenocarcinoma, NOS"),
  #                                c("Adenocarcinoma, NOS","Neuroendocrine carcinoma, NOS"),
  #                                c("Infiltrating duct carcinoma, NOS","Neuroendocrine carcinoma, NOS")),step_increase=0.1,map_signif_level = F,test = wilcox.test,color="black")+
  theme(legend.position="none")+
  scale_fill_manual(values = mycol)+
  xlab("")+ylab("Risk score")+
  labs(title =colnames(dat3)[3])+
  annotate("text",x=2,y=max(dat$risk_score),label = paste("one-way anova, p =",pvalue))



dat4<-dat[!dat$radiation_therapy=="",]
p4<-ggplot(data = dat4, mapping = aes(x = dat4[,6], y = dat4[,9], fill = dat4[,6])) +
  geom_boxplot()+
  theme_bw()+
  geom_signif(comparisons = list(c("NO","YES")),step_increase=0.1,map_signif_level = F,test = wilcox.test,color="black")+
  theme(legend.position="none")+
  scale_fill_manual(values = mycol[c(1,4)])+
  xlab("")+ylab("Risk score")+
  labs(title =colnames(dat3)[6])

# dat5<-dat[!dat$prior_malignancy.diagnoses=="",]
# p5<-ggplot(data = dat5, mapping = aes(x = dat5[,4], y = dat5[,9], fill = dat5[,4])) +
#   geom_boxplot()+
#   theme_bw()+
#   geom_signif(comparisons = list(c("no","yes")),step_increase=0.1,map_signif_level = F,test = wilcox.test,color="black")+
#   theme(legend.position="none")+
#   scale_fill_manual(values = mycol[c(2,9)])+
#   xlab("")+ylab("Risk score")+
#   labs(title =colnames(dat3)[5])

library(cowplot)
plot_grid(p1,p2,p3,p4,ncol = 2)
