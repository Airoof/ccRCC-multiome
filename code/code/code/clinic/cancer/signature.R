
setwd('E:\\work\\RCC')
library(Seurat)
rcc.combined<-readRDS('data/combined/rcc_combined_4_73477.rds')
meta<-readRDS('data/combined/meta.rds')
meta<-meta[colnames(rcc.combined),]
rcc.combined$subtype_cancer<-meta$subtype_cancer

####deal tcga data
data<-readRDS('data/tcga/data_rownamer_tpm.rds')
patient<-substring(colnames(data),14,14)
data<-data[,patient==0]
patient<-substring(colnames(data),16,16)
data<-data[,patient=='A']

#####tumor cell subtype signature
##(1) all cell signature
type1<-FindMarkers(rcc.combined,ident.1 = 'type 1',group.by = 'subtype_cancer',only.pos = T)
type1<-subset(type1,p_val_adj<0.05)
type1<-type1[order(type1$avg_log2FC,decreasing = T),]
sig1<-rownames(type1)[1:50]

type2<-FindMarkers(rcc.combined,ident.1 = 'type 2',group.by = 'subtype_cancer',only.pos = T)
type2<-subset(type2,p_val_adj<0.05)
type2<-type2[order(type2$avg_log2FC,decreasing = T),]
sig2<-rownames(type2)[1:10]

type3<-FindMarkers(rcc.combined,ident.1 = 'type 3',group.by = 'subtype_cancer',only.pos = T)
type3<-subset(type3,p_val_adj<0.05)
type3<-type3[order(type3$avg_log2FC,decreasing = T),]
sig3<-rownames(type3)[1:10]

type4<-FindMarkers(rcc.combined,ident.1 = 'type 4',group.by = 'subtype_cancer',only.pos = T)
type4<-subset(type4,p_val_adj<0.05)
type4<-type4[order(type4$avg_log2FC,decreasing = T),]
sig4<-rownames(type4)[1:50]

#ssgsea
library(GSVA)
sigs<-list(sig1,sig2,sig3,sig4)
gsva_matrix <- gsva(as.matrix(data), sigs, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
rownames(gsva_matrix)<-c('type_1_signature','type_2_signature','type_3_signature','type_4_signature')

#os analysis
surdata<-read.delim('data/tcga/TCGA-KIRC.survival.tsv',row.names = 1)
rownames(surdata)<-gsub('-','.',rownames(surdata))
sample<-intersect(rownames(surdata),colnames(data))
surdata<-surdata[sample,]
gsva_matrix<-t(gsva_matrix[,sample])

surdata<-cbind(surdata,gsva_matrix)

surdata$level<-ifelse(surdata$type_2_signature>median(surdata$type_2_signature),'high(n=261)','low(n=262)')
library("survival")
library("survminer")
sfit <- survfit(Surv(OS.time/365, OS)~level, data=surdata)
print(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE,palette = c("#Ff3a56", "#7b6afd"))

##########top 50 signature
gsva_matrix<-readRDS('duilie/cancersubtype/gsva_matrix_allcell_marker.rds')
ph<-readRDS('data/tcga/ph.rds')
rownames(ph)<-gsub('-','.',rownames(ph))

patient<-intersect(rownames(ph),colnames(gsva_matrix))
ph<-ph[patient,]
gsva_matrix<-gsva_matrix[,patient]

ph<-cbind(ph,t(gsva_matrix))

###trans
ph$level<-ifelse(ph$type_2_signature>median(ph$type_2_signature),'high','low')
ph$trans<-ph$primary_lymph_node_presentation_assessment
library(ggpubr)
ggboxplot(subset(ph,trans %in% c('YES',"NO")),'trans','type_4_signature',color = 'trans')+
  stat_compare_means(aes(group=trans), label = "p.signif")+
  theme(axis.text.x = element_text(angle = 60,vjust = 1,size = 10,hjust = 1))

ggplot(data = subset(ph,trans %in% c('YES',"NO")),aes(x=trans,y=type_2_signature))+
  geom_violin(trim = T,aes(fill=trans))+    #填充颜色，不剪尾
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),geom="pointrange", color = "red")+  #加均值及标准差
  theme_minimal()+    #背景
  geom_jitter(shape = 16, size = 1, position = position_jitter(0.1))+
  scale_fill_manual(values=c("#cb6c67", "#b883e2"))+
  stat_compare_means(aes(group=trans), label = "p.format")

###stage
ph$level<-ifelse(ph$type_2_signature>median(ph$type_2_signature),'high','low')
ph$stage<-ifelse(ph$pathologic_T %in% c('T1','T1a','T1b','T2','T2a','T2b'),'early','late')
library(ggpubr)
ggboxplot(ph,'stage','type_4_signature',color = 'stage')+
  stat_compare_means(aes(group=stage), label = "p.signif")+
  theme(axis.text.x = element_text(angle = 60,vjust = 1,size = 10,hjust = 1))

ggplot(data = ph,aes(x=stage,y=type_3_signature))+
  geom_violin(trim = T,aes(fill=stage))+    #填充颜色，不剪尾
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),geom="pointrange", color = "red")+  #加均值及标准差
  theme_minimal()+    #背景
  geom_jitter(shape = 16, size = 1, position = position_jitter(0.1))+
  scale_fill_manual(values=c("#cb6c67", "#b883e2"))+
  stat_compare_means(aes(group=stage), label = "p.format")

###stage
ph$level<-ifelse(ph$type_3_signature>median(ph$type_3_signature),'high','low')
ph$stage<-ph$pathologic_T
ph$stage<-substring(ph$stage,1,2)
ph$stage<-factor(ph$stage,levels = c('T1','T2','T3','T4'))
ggboxplot(ph,'stage','type_4_signature',color = 'stage')+
  stat_compare_means(aes(group=stage), label = "p.signif")+
  theme(axis.text.x = element_text(angle = 60,vjust = 1,size = 10,hjust = 1))

ggplot(data = ph,aes(x=stage,y=type_2_signature))+
  geom_violin(trim = T,aes(fill=stage))+    #填充颜色，不剪尾
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),geom="pointrange", color = "red")+  #加均值及标准差
  theme_minimal()+    #背景
  geom_jitter(shape = 16, size = 1, position = position_jitter(0.1))+
  #scale_fill_manual(values=c("#cb6c67", "#b883e2"))+
  stat_compare_means(aes(group=stage), label = "p.format")


#####stage sur
ph<-readRDS('data/tcga/ph.rds')
surdata<-read.delim('data/tcga/TCGA-KIRC.survival.tsv',row.names = 1)
sample<-intersect(rownames(ph),rownames(surdata))
surdata<-surdata[sample,]
ph<-ph[sample,]
plot_data<-cbind(surdata,ph)
#filter
patient<-substring(rownames(plot_data),14,14)
plot_data<-plot_data[patient==0,]
patient<-substring(rownames(plot_data),16,16)
plot_data<-plot_data[patient=='A',]

plot_data$stage<-ifelse(plot_data$pathologic_T %in% c('T1','T1a','T1b','T2','T2a','T2b'),'early(n=338)','late(n=192)')
library("survival")
library("survminer")
sfit <- survfit(Surv(OS.time/365, OS)~stage, data=plot_data)
print(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE,palette = c("#Ff3a56", "#7b6afd"))

plot_data$stage<-plot_data$pathologic_T
plot_data$stage<-substring(plot_data$stage,1,2)
#plot_data$stage<-factor(plot_data$stage,levels = c('T1','T2','T3','T4'))
plot_data$stage<-ifelse(plot_data$stage=='T1','T1(n=270)',plot_data$stage)
plot_data$stage<-ifelse(plot_data$stage=='T2','T2(n=68)',plot_data$stage)
plot_data$stage<-ifelse(plot_data$stage=='T3','T3(n=181)',plot_data$stage)
plot_data$stage<-ifelse(plot_data$stage=='T4','T4(n=11)',plot_data$stage)
library("survival")
library("survminer")
sfit <- survfit(Surv(OS.time/365, OS)~stage, data=plot_data)
print(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE)
