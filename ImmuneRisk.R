
load("处理后的cibresort_result.rdata")
load("all_risk_score.rdata")
result_ciber<-result_ciber[,all_risk_score$sample]
table(colnames(result_ciber)==all_risk_score$sample)
result_ciber<-result_ciber[-which(row.names(result_ciber)=="Plasma cells"),]
all_cor_immun<-vector()
for (i in 1:dim(result_ciber)[1]) {
  cor_dat<-cor.test(as.vector(result_ciber[i,]),all_risk_score$total_risk_score)
  pvalue<-cor_dat[["p.value"]]
  cor_value<-cor_dat[["estimate"]][["cor"]]
  all_cor_immun<-rbind(all_cor_immun,c(pvalue,cor_value))
}
row.names(all_cor_immun)<-row.names(result_ciber)
colnames(all_cor_immun)<-c("p_value","cor")
# save(all_cor_immun,file = "immun_risk_cor.rdata")
all_cor_immun<-all_cor_immun[all_cor_immun[,1]<0.01,]
aa<-all_cor_immun
aa[,1]<-aa[,2]
pheatmap(aa,
         # scale = "row",
         cluster_rows = T,
         cluster_cols = T,
         show_colnames = FALSE,
         show_rownames = T,
         color = c(colorRampPalette(rev(brewer.pal(9,"PiYG"))[5:9])(70)),
         cutree_cols = 2
)

library(survival)
library(survminer)
for (i in 1:dim(all_cor_immun)[1]) {
  high_low_test <- (result_ciber[row.names(all_cor_immun)[i],] > median(result_ciber[row.names(all_cor_immun)[i],]))
  high_low_test[high_low_test == TRUE] <- 'high'
  high_low_test[high_low_test == FALSE] <- 'low'
  all_risk_score$type <- high_low_test
  fit_km2 <- survfit(Surv(OS.m, status) ~ type, data = all_risk_score)     
  ggsurvplot(fit_km2, conf.int = F,pval = T,legend.title=paste("The fraction of ",row.names(all_cor_immun)[i],sep = ""),
             legend.labs=c("High", 
                           "Low"), risk.table = F, 
             palette = c("#C6147E",
                         "#329890"))
  
  dat<-data.frame(macro=as.vector(result_ciber[row.names(all_cor_immun)[i],]),
                  risk=risk_score_dat$total_risk_score,type="aa")
  ggplot(dat,aes(macro,risk,color=type))+
    geom_point(aes(color=type),size=1)+
    geom_smooth(method = lm)+
    scale_color_manual(values = pal_d3("category20")(20)[11])+
    # annotate(
    #   geom = "text", x = xrng[1], y = yrng[2]-2, 
    #   label = cc1, hjust = 0, vjust = 1, size = 3
    # )+
    # annotate(
    #   geom = "text", x = xrng[1]+0.5, y = yrng[2]-6, 
    #   label = cc2, hjust = 0, vjust = 1, size = 3
    # )+
    theme_bw()+
    labs(x="The fraction of Macrophages M0",y="Risk score")
  
}