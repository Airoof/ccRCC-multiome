
setwd('E:\\work\\RCC')

###
#rcc.m
library(Seurat)
rcc.m<-readRDS('immune/macrophage/rcc.m_2840.rds')
m1m2<-subset(rcc.m,celltype.refined %in% c('Macro0','Macro1','Macro2','Macro3'))
data<-m1m2@assays$RNA@data

#pathway
library(readxl)
mm <- read_excel("C:/Users/txw/Downloads/1-s2.0-S0092867420316135-mmc3.xlsx", 
                 skip = 1)
m1<-data.frame(na.omit(mm[,1]))[,1]
m2<-data.frame(na.omit(mm[,2]))[,1]
m3<-data.frame(na.omit(mm[,3]))[,1]
m4<-data.frame(na.omit(mm[,4]))[,1]

gene<-list()
gene[[colnames(mm)[1]]]<-m1
gene[[colnames(mm)[2]]]<-m2
gene[[colnames(mm)[3]]]<-m3
gene[[colnames(mm)[4]]]<-m4

##aucell
setwd('pathway/macro/')
library(AUCell)
cells_rankings <- AUCell_buildRankings(data, nCores=1, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(gene, cells_rankings)
saveRDS(cells_AUC,'cells_AUC.rds')

res<-cells_AUC@assays@data@listData[["AUC"]]
res<-data.frame(t(res))
res$cluster<-m1m2$seurat_clusters

m1score=c()
m2score=c()
proscore=c()
antiscore=c()

for(j in c(0,1,2,3)){
  t=subset(res,cluster==j)
  m1score<-c(m1score,mean(t[,1]))
  m2score<-c(m2score,mean(t[,2]))
  proscore<-c(proscore,mean(t[,3]))
  antiscore<-c(antiscore,mean(t[,4]))
}
df2<-data.frame(m1score=m1score,m2score=m2score,proscore=proscore,antiscore=antiscore)

library(ggplot2)
ggplot(res,aes(x=M1_Polarization,y=M2_Polarization,color=cluster))+geom_point()+theme_bw()+geom_point(data=df2,aes(x=m1score,y=m2score),color=c('red','green','blue','purple'),size=5,shape=17)
ggplot(res,aes(y=Anti_inflammatory,x=Pro_inflammatory,color=cluster))+geom_point()+geom_point(data=df2,aes(x=antiscore,y=proscore),color=c('red','green','blue','purple'),size=5,,shape=17)+theme_bw()


#############################################
##gsva
##mmc4
rcc.m<-readRDS('immune/macrophage/rcc.m_2840.rds')
m1m2<-subset(rcc.m,celltype.refined %in% c('TAM-HIF1A','TAM-CCL3','TAM-CD163','TAM-CD74'))

library(GSVA)
library(readxl)
mmc4_c <- read_excel("E:/work/RCC/pathway/macro/mmc4.xlsx", 
                     sheet = "C", skip = 1)
mmc4_c<-data.frame(mmc4_c)
m1<-as.character(mmc4_c[1,3:215])
m1<-m1[!is.na(m1)]
m2<-as.character(mmc4_c[2,3:215])
m2<-m2[!is.na(m2)]
m3<-as.character(mmc4_c[3,3:215])
m3<-m3[!is.na(m3)]
m4<-as.character(mmc4_c[4,3:215])
m4<-m4[!is.na(m4)]
m5<-as.character(mmc4_c[5,3:215])
m5<-m5[!is.na(m5)]
m6<-as.character(mmc4_c[6,3:215])
m6<-m6[!is.na(m6)]
m7<-as.character(mmc4_c[7,3:215])
m7<-m7[!is.na(m7)]
m8<-as.character(mmc4_c[8,3:215])
m8<-m8[!is.na(m8)]
m9<-as.character(mmc4_c[9,3:215])
m9<-m9[!is.na(m9)]

mmc4_b <- read_excel("E:/work/RCC/pathway/macro/mmc4.xlsx", 
                     sheet = "B", skip = 1)
mmc4_b<-data.frame(mmc4_b)
m10<-as.character(mmc4_b[1,2:73])
m10<-m10[!is.na(m10)]
m11<-as.character(mmc4_b[2,2:73])
m11<-m11[!is.na(m11)]

gene<-list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11)
data<-rcc.m@assays$RNA@data
res<-gsva(data,gene,method='gsva')
rownames(res)<-c(mmc4_c$Signature,mmc4_b$Signature)[1:11]
saveRDS(res,'pathway/macro/gsva_res_mmc4.rds')

#plot
res<-res[,colnames(m1m2)]
res<-data.frame(t(res))
res$cluster<-m1m2$celltype.refined

m1s=c()
m2s=c()
m3s=c()
m4s=c()
m5s=c()
m6s=c()
m7s=c()
m8s=c()
m9s=c()
m10s=c()
m11s=c()
p<-c()

for(j in c('TAM-HIF1A','TAM-CCL3','TAM-CD163','TAM-CD74')){
  t=subset(res,cluster==j)
  m1s<-c(m1s,mean(t[,1]))
  m2s<-c(m2s,mean(t[,2]))
  m3s<-c(m3s,mean(t[,3]))
  m4s<-c(m4s,mean(t[,4]))
  m5s<-c(m5s,mean(t[,5]))
  m6s<-c(m6s,mean(t[,6]))
  m7s<-c(m7s,mean(t[,7]))
  m8s<-c(m8s,mean(t[,8]))
  m9s<-c(m9s,mean(t[,9]))
  m10s<-c(m10s,mean(t[,10]))
  m11s<-c(m11s,mean(t[,11]))
}

for(i in 1:11){
  kruskal.test(as.formula(paste0(colnames(res)[i],'~cluster')), res)->t
  p<-c(p,t$p.value)
}
df2<-data.frame(IFNg_stimmed_MDM=m1s,IFNg_Response=m2s,Res_to_Type_I_IFN=m3s,
                Antigen_Proc_and_Pres=m4s,Complement=m5s,Proteasome=m6s,
                IL6_STAT3_Signaling=m7s,Fc_Receptor_Signaling=m8s,TGF_Beta_Signaling=m9s,
                M1_Curated=m10s,M2_Curated=m11s)
df2[5,]<-p

ggplot(res,aes(x=M1.Curated,y=M2.Curated,color=cluster))+geom_point(size=2)+theme_bw()+geom_point(data=df2,aes(x=M1_Curated,y=M2_Curated),color=c('red','green','#0092C7','purple'),size=8,shape=17)

####plot
library(pheatmap)
rownames(df2)<-c('TAM-HIF1A','TAM-CCL3','TAM-CD163','TAM-CD74')
pheatmap(t(df2), cluster_rows = F,cluster_cols = F,show_colnames = T)

#single cell
tmp<-res
anno_col<-tmp[,11:12]
anno_col<-anno_col[order(anno_col$cluster),]
cell_name<-rownames(anno_col)
anno_col<-data.frame(anno_col[,2])
rownames(anno_col)<-cell_name
colnames(anno_col)<-'cluster'

for(i in 1:11){
  tmp[,i]<-(tmp[,i]-min(tmp[,i]))/(max(tmp[,i])-min(tmp[,i]))
}
pheatmap(t(tmp[rownames(anno_col),-12]),cluster_rows = F,cluster_cols = F,show_colnames = F,annotation_col = anno_col,scale = 'none')

#####vlnplot
library(ggplot2)
library(paletteer)
paletteer_d("RColorBrewer::Paired")
ggplot(data = res,aes(x=cluster,y=M2.Curated))+
  geom_violin(trim = T,aes(fill=cluster))+    
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),geom="pointrange", color = "red")+ 
  theme_minimal()+    #±³¾°
  scale_fill_manual(values=c("#FB9A99FF", "#B2DF8AFF", "#A6CEE3FF",'#CAB2D6FF'))

####dotplot
df2<-data.frame(IFNg_stimmed_MDM=m1s,IFNg_Response=m2s,Res_to_Type_I_IFN=m3s,
                Antigen_Proc_and_Pres=m4s,Complement=m5s,Proteasome=m6s,
                IL6_STAT3_Signaling=m7s,Fc_Receptor_Signaling=m8s,TGF_Beta_Signaling=m9s,
                M1_Curated=m10s,M2_Curated=m11s)
ggplot(res,aes(x=M1.Curated,y=M2.Curated,color=cluster))+geom_point(size=2)+theme_bw()+scale_color_manual(values = c("#F8766D", '#7CAE00', "#00BFC4","#8167F5"))+
  geom_point(data=df2,aes(x=M1_Curated,y=M2_Curated),color=c("red", 'green', "blue","purple"),size=6,shape=17)

#density plot
ggplot(data = tmp,aes(x=M2.Curated,fill=cluster))+geom_density(alpha=0.3)+
  theme_minimal()
  #scale_fill_manual(values=c("#FB9A99FF", "#B2DF8AFF", "#A6CEE3FF",'#CAB2D6FF'))

#rigde plot 
library(ggridges)
library(sciplot)
tmp$cluster<-factor(tmp$cluster,levels = c('TAM-HIF1A','TAM-CCL3','TAM-CD163','TAM-CD74'))
ggplot(tmp,aes(x = M2.Curated,y = cluster,fill = cluster))+
  geom_density_ridges(scale = 2,size = 0.75,alpha = 0.75)+
  theme_bw()


##########ÈÈÍ¼
library(pheatmap)
anno_row<-data.frame(matrix(0,ncol = 2,nrow = 43))
anno_row$X2<-c( 'PPARG', 'THBS1', 'VEGFA','CD44','HES1','SPP1',
                'CTSB', 'CTSL', 'FABP5', 'GPNMB','LGALS3', 'PLA2G7',
                'CCL3', 'CCL3L1','CCL4L2','KLF6', 'EGR3',
                'CD163', 'CD163L1','MRC1', 'IGF1', 'F13A1','HNMT','MAF','MS4A4A',
                'CD74', 'FCGR3A','HLA-DPA1','MT-CO2','C1QA','CD68', 'FCGR1A', 'FCGR2A', 'APOL6','CD14',
                'MKI67', 'TOP2A', 'CCNA2', 'CDK1','HIST1H4C','HMGN2', 'STMN1', 'TYMS')
anno_row$X1<-c(rep('Angiogenic signature',6),rep('Lipid-related genes',6),rep('Chemokines',5),rep('M2 signature',8),rep('M1 signature',10),rep('Proliferation',8))

macro<-subset(rcc.m,celltype.refined %in% c('TAM-HIF1A','TAM-CCL3','TAM-CD163','TAM-CD74','TAM-MKI67'))
data_macro<-macro@assays$RNA@data
data_macro<-data_macro[anno_row$X2,]
df<-data.frame(matrix(0,ncol = 5,nrow = 43))
rownames(df)<-anno_row$X2
colnames(df)<-c('TAM-HIF1A','TAM-CCL3','TAM-CD163','TAM-CD74','TAM-MKI67')
for(i in c('TAM-HIF1A','TAM-CCL3','TAM-CD163','TAM-CD74','TAM-MKI67')){
  cells<-colnames(subset(macro,celltype.refined==i))
  tmp_data<-data_macro[,cells]
  df[,i]<-rowMeans(tmp_data)
}

rownames(anno_row)<-anno_row$X2
anno_row<-data.frame(signature=anno_row$X1)
rownames(anno_row)<-rownames(data_macro)
anno_row$signature<-factor(anno_row$signature,levels = c('Angiogenic signature','Lipid-related genes','Chemokines','M2 signature','M1 signature','Proliferation'),ordered = T)

pheatmap(df,cluster_rows = F,cluster_cols = F,scale = 'row',annotation_row = anno_row,
         legend_labels = c('Angiogenic signature','Lipid-related genes','Chemokines','M2 signature','M1 signature','Proliferation'))
