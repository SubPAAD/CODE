

library(reshape2)
load("PAAD_symbol_exp.rdata")
dd<-apply(tumor_exp,1,function(x){sum(x>0)})
diff_exp<-tumor_exp[dd>10,]
source("utils.R")
pathways <- gmtPathways("h.all.v7.2.symbols.gmt")
load("hallmark_boxplot.rdata")
pathways<-pathways[names(all_p)]
hall_gene<-unique(unlist(pathways))

dat_hall<-diff_exp[row.names(diff_exp)%in%hall_gene,]

load("GSVA_处理后结果.rdata")
dat_hall<-dat_hall[,all_sample_anno$sample]
table(colnames(dat_hall)==all_sample_anno$sample)
colnames(dat_hall)<-paste("Cluster ",as.character(all_sample_anno$cluster),sep = "")
dat_hall<-as.data.frame(dat_hall)
hall_pvalue<-vector()
# hall_name<-vector()
for (i in 1:dim(dat_hall)[1]) {
  plot_dat2<-melt(as.matrix(dat_hall[i,]))
  fit <- aov(value ~ Var2, data=plot_dat2)
  pvalue<-summary(fit)[[1]][1,5]
  pvalue<-signif(pvalue,3)
  hall_pvalue<-c(hall_pvalue,pvalue)
}
names(hall_pvalue)<-row.names(dat_hall)
# save(hall_pvalue,file = "差异表达致癌基因.rdata")
tf_dat<-read.csv("Homo_sapiens_TF.txt",sep = "\t",header = T)
table(tf_dat$Symbol%in%names(hall_pvalue[hall_pvalue<0.01]))
hall_tf<-tf_dat$Symbol[tf_dat$Symbol%in%names(hall_pvalue[hall_pvalue<0.01])]

load("TF_TG_interaction.rdata")
hall_tf_tg<-tf_tg[tf_tg$TF.Name%in%hall_tf,]
hall_tf_tg<-hall_tf_tg[hall_tf_tg$TG.Name%in%names(hall_pvalue[hall_pvalue<0.01]),]
hall_tf_tg<-hall_tf_tg[!hall_tf_tg$TG.Name%in%hall_tf,]

aa<-vector()
bb<-vector()
for (j in 1:dim(hall_tf_tg)[1]) {
  aa<-c(aa,cor.test(as.vector(t(dat_hall[hall_tf_tg[j,1],])),as.vector(t(dat_hall[hall_tf_tg[j,2],])))[["p.value"]])
  
  bb<-c(bb,cor.test(as.vector(t(dat_hall[hall_tf_tg[j,1],])),as.vector(t(dat_hall[hall_tf_tg[j,2],])))[["estimate"]][["cor"]])
}
hall_tf_tg$pvalue<-aa
hall_tf_tg$cor<-bb
save(hall_tf_tg,file = "转录调控网络.rdata")
table(hall_tf_tg$pvalue<0.01&hall_tf_tg$cor>0.6)
hall_tf_tg2<-hall_tf_tg[hall_tf_tg$pvalue<0.01&hall_tf_tg$cor>0.6,]

all_inter<-hall_tf_tg2[,c(1,2)]
all_inter<-all_inter[!duplicated(all_inter),]
all_inter_anno<-data.frame(symbol=c(unique(all_inter$TF.Name),unique(all_inter$TG.Name)),
                           type=c(rep("TF",length(unique(all_inter$TF.Name))),
                                  rep("TG",length(unique(all_inter$TG.Name)))))
write.table(all_inter,"all_inter.txt",sep = "\t",quote = F,row.names = F)
write.table(all_inter_anno,"all_inter_anno.txt",sep = "\t",quote = F,row.names = F)
