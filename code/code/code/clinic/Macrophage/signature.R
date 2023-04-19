
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
tam_ccl3<-FindMarkers(rcc.combined,ident.1 = 'Macro1',group.by = 'subtype_cancer',only.pos = T)
tam_ccl3<-subset(tam_ccl3,p_val_adj<0.05)
tam_ccl3<-tam_ccl3[order(tam_ccl3$avg_log2FC,decreasing = T),]
sig_tam_ccl3<-rownames(tam_ccl3)[1:50]

tam_cd163<-FindMarkers(rcc.combined,ident.1 = 'Macro2',group.by = 'subtype_cancer',only.pos = T)
tam_cd163<-subset(tam_cd163,p_val_adj<0.05)
tam_cd163<-tam_cd163[order(tam_cd163$avg_log2FC,decreasing = T),]
sig_tam_cd163<-rownames(tam_cd163)[1:50]

tam_cd74<-FindMarkers(rcc.combined,ident.1 = 'Macro3',group.by = 'subtype_cancer',only.pos = T)
tam_cd74<-subset(tam_cd74,p_val_adj<0.05)
tam_cd74<-tam_cd74[order(tam_cd74$avg_log2FC,decreasing = T),]
sig_tam_cd74<-rownames(tam_cd74)[1:50]

tam_hif1a<-FindMarkers(rcc.combined,ident.1 = 'Macro0',group.by = 'subtype_cancer',only.pos = T)
tam_hif1a<-subset(tam_hif1a,p_val_adj<0.05)
tam_hif1a<-tam_hif1a[order(tam_hif1a$avg_log2FC,decreasing = T),]
sig_tam_hif1a<-rownames(tam_hif1a)[1:50]

tam_mki67<-FindMarkers(rcc.combined,ident.1 = 'Macrophage(proliferating)',group.by = 'subtype_cancer',only.pos = T)
tam_mki67<-subset(tam_mki67,p_val_adj<0.05)
tam_mki67<-tam_mki67[order(tam_mki67$avg_log2FC,decreasing = T),]
sig_tam_mki67<-rownames(tam_mki67)[1:50]

#ssgsea
library(GSVA)
sigs<-list(sig_tam_ccl3,sig_tam_cd163,sig_tam_cd74,sig_tam_hif1a,sig_tam_mki67)
gsva_matrix <- gsva(as.matrix(data), sigs, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
rownames(gsva_matrix)<-c('TAM_CCL3_signature','TAM_CD163_signature','TAM_CD74_signature','TAM_HIF1A_signature','TAM_MKI67_signature')
saveRDS(gsva_matrix,'duilie/macro//gsva_matrix_allcell_50.rds')

#os analysis
surdata<-read.delim('data/tcga/TCGA-KIRC.survival.tsv',row.names = 1)
rownames(surdata)<-gsub('-','.',rownames(surdata))
sample<-intersect(rownames(surdata),colnames(data))
surdata<-surdata[sample,]
gsva_matrix<-t(gsva_matrix[,sample])

surdata<-cbind(surdata,gsva_matrix)

surdata$level<-ifelse(surdata$TAM_HIF1A_signature>quantile(surdata$TAM_HIF1A_signature)[4],'high',surdata$level)
surdata$level<-ifelse(surdata$TAM_HIF1A_signature<quantile(surdata$TAM_HIF1A_signature)[2],'low',surdata$level)
library("survival")
library("survminer")
sfit <- survfit(Surv(OS.time/365, OS)~level, data=surdata)
print(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE)


##########top 50 signature
gsva_matrix<-readRDS('duilie/macro/gsva_matrix_allcell_50.rds')
ph<-readRDS('data/tcga/ph.rds')
rownames(ph)<-gsub('-','.',rownames(ph))

patient<-intersect(rownames(ph),colnames(gsva_matrix))
ph<-ph[patient,]
gsva_matrix<-gsva_matrix[,patient]

ph<-cbind(ph,t(gsva_matrix))

###trans
ph$level<-ifelse(ph$TAM_MKI67_signature>median(ph$TAM_MKI67_signature),'high','low')
ph$trans<-ph$primary_lymph_node_presentation_assessment
library(ggpubr)
ggboxplot(subset(ph,trans %in% c('YES',"NO")),'trans','TAM_MKI67_signature',color = 'trans')+
  stat_compare_means(aes(group=trans), label = "p.signif")+
  theme(axis.text.x = element_text(angle = 60,vjust = 1,size = 10,hjust = 1))

ggplot(data = subset(ph,trans %in% c('YES','NO')),aes(x=trans,y=TAM_HIF1A_signature))+
  geom_violin(trim = T,aes(fill=trans)) +   
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),geom="pointrange", color = "red")+  
  theme_minimal()+
  stat_compare_means(aes(group=trans), label = "p.signif")+
  theme(panel.grid =element_blank())+
  theme(axis.line = element_line(size=0.5, colour = "black"))   
  # scale_fill_manual(values=c("#FB9A99FF", "#B2DF8AFF"))

###stage
ph$level<-ifelse(ph$TAM_MKI67_signature>median(ph$TAM_MKI67_signature),'high','low')
ph$stage<-ifelse(ph$pathologic_T %in% c('T1','T1a','T1b','T2','T2a','T2b'),'early','late')
library(ggpubr)
ggboxplot(ph,'stage','TAM_MKI67_signature',color = 'stage')+
  stat_compare_means(aes(group=stage), label = "p.signif")+
  theme(axis.text.x = element_text(angle = 60,vjust = 1,size = 10,hjust = 1))

###stage
ph$level<-ifelse(ph$TAM_MKI67_signature>median(ph$TAM_MKI67_signature),'high','low')
ph$stage<-ph$pathologic_T
ph$stage<-substring(ph$stage,1,2)
ph$stage<-factor(ph$stage,levels = c('T1','T2','T3','T4'))
ggboxplot(ph,'stage','TAM_MKI67_signature',color = 'stage')+
  stat_compare_means(aes(group=stage), label = "p.signif")+
  theme(axis.text.x = element_text(angle = 60,vjust = 1,size = 10,hjust = 1))
