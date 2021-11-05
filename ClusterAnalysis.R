

library(survival)
library('survminer')


survival_cancer1<-read.csv("TCGA-PAAD.survival.tsv",header = T,sep = "\t")
survival_cancer1<-survival_cancer1[!survival_cancer1$OS.time==0,]
smap_nam<-intersect(survival_cancer1$sample,all_sample_anno$sample)
survival_cancer<-survival_cancer1[match(smap_nam,survival_cancer1$sample),c(1,4,2)]
colnames(survival_cancer)<-c("sample","OS.m","status")
cluster_id<-all_sample_anno[match(smap_nam,all_sample_anno$sample),2]
table(survival_cancer$sample==smap_nam)
survival_cancer$type<-cluster_id

fit_km <- survfit(Surv(OS.m, status) ~ type, data = survival_cancer)     
ggsurvplot(fit_km, conf.int = F,pval = T,legend.title="Sample_type",
           risk.table = T, 
           palette = brewer.pal(11,"BrBG")[c(2,4,7,9)])
