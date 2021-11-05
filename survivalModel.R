library(ggplot2)
library(RColorBrewer)

load("转录调控网络.rdata")
load("PAAD_symbol_exp.rdata")
table(hall_tf_tg$pvalue<0.01&hall_tf_tg$cor>0.6)
hall_tf_tg2<-hall_tf_tg[hall_tf_tg$pvalue<0.01&hall_tf_tg$cor>0.6,]
feature<-c(unique(hall_tf_tg2$TF.Name),unique(hall_tf_tg2$TG.Name))

#################单因素cox回归

survival_cancer1<-read.csv("TCGA-PAAD.survival.tsv",header = T,sep = "\t")
survival_cancer1<-survival_cancer1[!survival_cancer1$OS.time==0,]
smap_nam<-intersect(survival_cancer1$sample,colnames(tumor_exp))
survival_cancer<-survival_cancer1[match(smap_nam,survival_cancer1$sample),c(1,4,2)]
dat_exp<-tumor_exp[,smap_nam]
colnames(survival_cancer)<-c("sample","OS.m","status")
feature_exp<-as.data.frame(t(dat_exp[feature,]))
table(row.names(feature_exp)==survival_cancer$sample)

survival_cancer2<-cbind(survival_cancer,feature_exp)
colnames(survival_cancer2)<-gsub(colnames(survival_cancer2), pattern = '-', replacement = '_')
uni_cox_in_bulk<-function(gene_list, survival_info_df){
  library('survival')
  gene_list <- gsub(gene_list, pattern = '-', replacement = '_')
  uni_cox <- function(single_gene){
    formula <- as.formula(paste0('Surv(OS.m, status)~', single_gene))
    surv_uni_cox <- summary(coxph(formula, data = survival_cancer2))
    ph_hypothesis_p <- cox.zph(coxph(formula, data = survival_cancer2))$table[1,3]
    if (surv_uni_cox$coefficients[,5]<0.05 & ph_hypothesis_p>0.05){  #get the pvalue
      single_cox_report<- data.frame('uni_cox_sig_genes'=single_gene,
                                     'beta'=surv_uni_cox$coefficients[,1],
                                     'Hazard_Ratio'=exp(surv_uni_cox$coefficients[,1]),
                                     'z_pvalue'=surv_uni_cox$coefficients[,5],
                                     'Wald_pvalue'=as.numeric(surv_uni_cox$waldtest[3]),
                                     'Likelihood_pvalue'=as.numeric(surv_uni_cox$logtest[3]))
      single_cox_report
    }
  }
  uni_cox_list <-lapply(gene_list, uni_cox)
  do.call(rbind, uni_cox_list)
}
uni_cox_df<- uni_cox_in_bulk(gene_list =colnames(survival_cancer2)[4:dim(survival_cancer2)[2]], survival_info_df = survival_cancer2)
####################拉索模型筛选特征
library(glmnet)
x <- as.matrix(survival_cancer2[,uni_cox_df$uni_cox_sig_genes])
y <- survival_cancer2[,c('OS.m', 'status')]
names(y) <- c('time', 'status')
y$time <- as.double(y$time)
y$status <- as.double(y$status)
y <- as.matrix(survival::Surv(y$time, y$status))
lasso <- glmnet(x, y, family='cox', type.measure = 'deviance')
par(mfrow=c(2,1))
plot(lasso,xvar = "lambda",col= colorRampPalette(brewer.pal(n = 11, "BrBG"))(40))
lasso_fit <- cv.glmnet(x, y, family='cox', type.measure = 'deviance')
plot(lasso_fit)
coefficient <- coef(lasso_fit, s=lasso_fit$lambda.min)
Active.Index <- which(as.numeric(coefficient) != 0)
active.coefficients <- as.numeric(coefficient)[Active.Index]
sig_gene_multi_cox <- rownames(coefficient)[Active.Index]
save(sig_gene_multi_cox,file = "拉索特征筛选结果.rdata")


########################cox多因素回归模型
train_site<-sample(1:dim(feature_exp)[1],0.6*dim(feature_exp)[1])
# save(train_site,file = "train_site.rdata")
train_dat<-feature_exp[train_site,]
test_dat<-feature_exp[-train_site,]
survival_cancer2<-cbind(survival_cancer[train_site,],train_dat)
colnames(survival_cancer2)<-gsub(colnames(survival_cancer2), pattern = '-', replacement = '_')
library(survminer)
library(survival)
formula_for_multivariate <-as.formula(paste0('Surv(OS.m, status)~', paste(sig_gene_multi_cox, sep = '', collapse = '+')))
multi_variate_cox <-coxph(formula_for_multivariate, data = survival_cancer2)
#check if variances are supported by PH hypothesis.
ph_hypo_multi <- cox.zph(multi_variate_cox)
# ggcoxzph(ph_hypo_multi)
#The last row of the table records the test results on the GLOBAL model. Delete it.
ph_hypo_table <- ph_hypo_multi$table[-nrow(ph_hypo_multi$table),]
#Remove variances not supported by ph hypothesis and perform the 2nd regression.
formula_for_multivariate <- as.formula(paste0('Surv(OS.m, status)~', paste(rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05], sep = '', collapse = '+')))
multi_variate_cox_2 <- coxph(formula_for_multivariate, data = survival_cancer2)
ggforest(model = multi_variate_cox_2, data = survival_cancer2, main = 'Hazard ratios of candidate genes', fontsize = 1)
############列线图
mycol<-pal_d3("category10")(10)
library(rms)
dd<-datadist(survival_cancer2)
options(datadist="dd")
cox_report<-summary(multi_variate_cox_2)$coefficients
cox_report<-cox_report[cox_report[,5]<0.05,]
candidate_genes_for_cox2 <- c(row.names(cox_report))

formula_for_multivariate2<-as.formula(paste0('Surv(OS.m, status)~', paste(candidate_genes_for_cox2, sep = '', collapse = '+')))
f2 <- cph(formula_for_multivariate,data = survival_cancer2,x=T,y=T,surv = T)

library(regplot)

pbcox<-coxph(formula_for_multivariate, data = survival_cancer2)
regplot(pbcox,
        #对观测2的六个指标在列线图上进行计分展示
        # observation=survival_cancer2[3,], #也可以不展示
        
        #预测3年和5年的死亡风险，此处单位是day
        failtime = c(180,365), 
        prfail = T, #cox回归中需要TRUE
        showP = T, #是否展示统计学差异
        droplines = F,#观测2示例计分是否画线
        # colors = mycol, #用前面自己定义的颜色
        rank="sd", #根据统计学差异的显著性进行变量的排序
        interval="confidence",
        points = T,
        plots=c("density","boxes")) #展示观测的可信区间
#########验证

cal1 <- calibrate(f2, cmethod='KM', method="boot", 
                  u=180, m=33, B=150)
cal2 <- calibrate(f2, cmethod='KM', method="boot", 
                  u=365, m=33, B=150)
# cal3 <- calibrate(f2, cmethod='KM', method="boot", 
#                   u=545, m=25, B=150)
#####半年和1年联合绘制
plot(cal1,lwd = 2,lty = 1,errbar.col = "#C6147E",
     bty = "l", #只画左边和下边框
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = "#C6147E",
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.2)

plot(cal2,lwd = 2,lty = 1,errbar.col = "#4C931B",
     xlim = c(0,1),ylim= c(0,1),col = "#4C931B",add = T)

abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

# plot(cal3,lwd = 2,lty = 1,errbar.col = "grey",
#      xlim = c(0,1),ylim= c(0,1),col = "grey",add = T)

legend("topleft", #图例的位置
       legend = c("6-months","12-months"), #图例文字
       col =c("#C6147E","#4C931B",mycol[4]), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.2,#图例字体大小
       bty = "n")#不显示图例边框
dev.off()
#################ROC曲线图
riskscore <- function(survival_cancer_df, candidate_genes_for_cox, cox_report) {
  all_risk_score<-vector()
  library('dplyr')
  risk_score_table <- survival_cancer_df[,candidate_genes_for_cox]
  for(each_sig_gene in colnames(risk_score_table)){
    risk_score <- risk_score_table[,each_sig_gene]*cox_report[each_sig_gene,1]
    all_risk_score<-cbind(all_risk_score,risk_score)
  }
  total_risk_score<-rowSums(all_risk_score)
  risk_score_table <- cbind(risk_score_table, total_risk_score,survival_cancer_df[,c('sample','OS.m','status')])
  risk_score_table <- risk_score_table[,c('sample','OS.m','status', candidate_genes_for_cox, 'total_risk_score')]
  return(risk_score_table)
}
cox_report<-summary(multi_variate_cox_2)$coefficients
cox_report<-cox_report[cox_report[,5]<0.05,]
candidate_genes_for_cox2 <- c(row.names(cox_report))
risk_score_table_multi_cox2 <- riskscore(survival_cancer_df=survival_cancer2, 
                                         candidate_genes_for_cox=candidate_genes_for_cox2, 
                                         cox_report=cox_report)
# save(risk_score_table_multi_cox2,file = "all_risk_score.rdata")

multi_ROC <- function(time_vector, risk_score_table){
  library('survivalROC')
  single_ROC <- function(single_time){
    for_ROC <- survivalROC(Stime = risk_score_table$OS.m,
                           status = risk_score_table$status,
                           marker = risk_score_table$total_risk_score,
                           predict.time = single_time, method = 'KM')
    data.frame('True_positive'=for_ROC$TP, 'False_positive'=for_ROC$FP, 
               'Cut_values'=for_ROC$cut.values, 'Time_point'=rep(single_time, length(for_ROC$TP)),
               'AUC'=rep(for_ROC$AUC, length(for_ROC$TP)))
  }
  multi_ROC_list<- lapply(time_vector, single_ROC)
  do.call(rbind, multi_ROC_list)
}
#We evaluate 11 AUCs between 3-5 years.
for_multi_ROC<- multi_ROC(time_vector = c(365*seq(0.5,1,0.1)), risk_score_table = risk_score_table_multi_cox2)
library(ggplot2)
library(ggsci)
ggplot(for_multi_ROC,aes(x = False_positive, y = True_positive, color = as.factor(Time_point))) + 
  geom_path(size=1) + 
  theme_bw()+
  # theme(panel.grid.major = element_blank(),   
  #       panel.grid.minor = element_blank(),
  #       #不显示网格
  #       axis.line = element_line(colour = "black",size=0.5)
  # )+
  labs(x = "False positive rate", y = "Ture positive rate", title ="ROC curve") +
  scale_colour_manual(values = brewer.pal(9,"RdPu")[2:7])+ ###############"YlGn","Reds","RdPu"
  annotate("text",x = 0.75, y = 0.15,
           label = paste("AUC max = ", round(max(for_multi_ROC$AUC),2), '\n', 'AUC max time = ',for_multi_ROC$Time_point[which(for_multi_ROC$AUC==max(for_multi_ROC$AUC))[1]], 'days', sep = ''))
#################绘制生存曲线
cut.off <-median(risk_score_table_multi_cox2$total_risk_score)
high_low <- (risk_score_table_multi_cox2$total_risk_score > cut.off)
high_low[high_low == TRUE] <- 'high'
high_low[high_low == FALSE] <- 'low'
risk_score_table_multi_cox2$high_low <- high_low
save(cut.off,cox_report,risk_score_table_multi_cox2,file = "PAAD_cox_result.rdata")

pbcox<-coxph(formula_for_multivariate2, data = survival_cancer2)
library(ggrisk)
ggrisk(pbcox,color.A = c(low = "#20837B", high = "#B47626"),
       color.C = c(low = "#20837B", median = "white", high = "#B47626"),
       color.B = c(code.0 = "#20837B", code.1 = "#B47626"),
       cutoff.value = "median",size.points = 1)

library('survminer')
fit_km <- survfit(Surv(OS.m, status) ~high_low, data = risk_score_table_multi_cox2)     
ggsurvplot(fit_km, conf.int = F,pval = T,legend.title="Total risk score",
           legend.labs=c(paste0('>',as.character(round(cut.off,2))), paste0('<=',as.character(round(cut.off,2)))), risk.table = F, 
           palette = c("#20837B",
                       "#B47626"))
#################测试集验证
# save(candidate_genes_for_cox2,cox_report,cut.off,file = "生存分析辅助数据.rdata")
vali_dat<-as.data.frame(t(test_dat))
sur_dat<-survival_cancer
sur_dat<-sur_dat[sur_dat$sample%in%colnames(vali_dat),]
sur_dat<-sur_dat[match(colnames(vali_dat),sur_dat$sample),]
# row.names(vali_dat)<-gsub(row.names(vali_dat), pattern = '-', replacement = '_')
vali_candidate<-as.data.frame(t(vali_dat[candidate_genes_for_cox2,]))
table(row.names(vali_candidate)==sur_dat$sample)
survival_cancer_test<-cbind(sur_dat,vali_candidate)

pbcox_test<-coxph(formula_for_multivariate2, data = survival_cancer_test)
cox_report2<-summary(pbcox_test)$coefficients
# ggrisk(pbcox_test,color.A = c(low = "#83C2E3", high = "#F1E7C1"),
#        color.C = c(low = "#83C2E3", median = "white", high = "#F1E7C1"),
#        color.B = c(code.0 = "#83C2E3", code.1 = "#F1E7C1"),
#        cutoff.value = "median",size.points = 1)
ggrisk(pbcox_test,color.A = c(low = "#20837B", high = "#B47626"),
       color.C = c(low = "#20837B", median = "white", high = "#B47626"),
       color.B = c(code.0 = "#20837B", code.1 = "#B47626"),
       cutoff.value = "median",size.points = 1)

risk_score_test <- riskscore(survival_cancer_test, candidate_genes_for_cox2, cox_report)

high_low_test <- (risk_score_test$total_risk_score > median(risk_score_test$total_risk_score))
high_low_test[high_low_test == TRUE] <- 'high'
high_low_test[high_low_test == FALSE] <- 'low'
risk_score_test$type <- high_low_test

# risk_score_test$status[which(risk_score_test$OS.m > AUC_max_time)] <- 0
# risk_score_test$OS.m[which(risk_score_test$OS.m > AUC_max_time)] <- AUC_max_time
fit_km2 <- survfit(Surv(OS.m, status) ~ type, data = risk_score_test)     
ggsurvplot(fit_km2, conf.int = F,pval = T,legend.title="Validation total risk score",
           legend.labs=c(paste0('>',as.character(round(cut.off,2))), paste0('<=',as.character(round(cut.off,2)))), risk.table = F, 
           palette = c("#20837B",
                       "#B47626"))

