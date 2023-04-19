
setwd('E:\\work\\RCC\\immune\\macrophage')
library(Seurat)
library(ggplot2)
rcc.combined<-readRDS('/public/home/xwtang/work/RCC/data/intergrated/rcc_combined_4.rds')
rcc.m<-subset(rcc.combined,celltype.refined=='Myeloid cell')
DefaultAssay(rcc.m)<-'integrated'

#plot
pdf('/public/home/xwtang/work/RCC/immune_cell/macrophage/pic/test.pdf',width = 8,height = 6)
FeaturePlot(rcc.m,c('LRP2','CUBN','CD247','PTPRC'))
VlnPlot(rcc.m,c('PTPRC','CD247','EPCAM'))
DimPlot(rcc.m,label = T)
dev.off()

rcc.m<-FindVariableFeatures(rcc.m,nFeatures=2000)
rcc.m <- ScaleData(rcc.m, verbose = FALSE)
rcc.m <- RunPCA(rcc.m, npcs = 30, verbose = FALSE)
rcc.m <- RunUMAP(rcc.m, reduction = "pca", dims = 1:15)
rcc.m <- FindNeighbors(rcc.m, reduction = "pca", dims = 1:15) 
rcc.m <- FindClusters(rcc.m, resolution = 1)

#plot
pdf('/public/home/xwtang/work/RCC/immune_cell/macrophage/pic/m_cluster.pdf',width = 8,height = 6)
DimPlot(rcc.m,label = T)
dev.off()

pdf('/public/home/xwtang/work/RCC/immune_cell/macrophage/pic/m_patient.pdf',width = 8,height = 6)
DimPlot(rcc.m,group.by = 'orig.ident')
dev.off()

pdf('/public/home/xwtang/work/RCC/immune_cell/T/pic/test.pdf',width = 8,height = 6)
FeaturePlot(rcc.t,'CD4')
dev.off()

###
cg=c('HLA-DPA1','APOE','ZNF331','SLC16A10','CD163L1','F13A1','NLGN1','RORA','CD247','TOX','CCL3','CCL3L1','CLEC9A','NEGR1','TOP2A','MKI67','FCN1','RIPOR2')
pdf('/public/home/xwtang/work/RCC/immune_cell/T/pic/all_marker.pdf')
DotPlot(rcc.m, features = unique(cg),
        assay='RNA'  )  + coord_flip()+ RotatedAxis() 
dev.off()



########
#0 M2  cd206(MRC1)/cd163
#1 M2
#2 M2
#3 M2
#4 M1
#5 M2
#6 Macrophage-EPI
#7 Macrophage-EPI
#8 M1 cd206-(MCR1)/cd163-/FCGR3A+
#9 Macrophage-EPI
#10 CD16+ monocyte cd16(fcgr3a)/cd68/cd52/s100a8/fcn1
#11 DC CLEC9A/adam19/clnk/flt3  mhc分子高表达
#12 TAM(proliferative)

#0 Macro0
#1 Macro1
#2 Macro2
#3 Macro3
#4 Macro2
#5 Macrophage-EPC
#6 Macrophage-EPC
#7 Macro0
#8 Macro2
#9 DC CLEC9A/adam19/clnk/flt3
#10 Macrophage-CD247+ cell
#11 CD16+ monocyte cd16(fcgr3a)/cd68/cd52/s100a8/fcn1
#12 Macrophage(proliferating)


type<-c('0'='Macro0','1'='Macro1','2'='Macro2','3'='Macro3','4'='Macro2','5'='Macrophage-EPC',
        '6'='Macrophage-EPC','7'='Macro0','8'='Macro2','9'='DC','10'='Macrophage-CD247+ cell',
        '11'='CD16+ monocyte','12'='Macrophage(proliferating)')
rcc.m[['celltype.refined']] = unname(type[rcc.m@meta.data$seurat_clusters])
saveRDS(rcc.m,'E:\\work\\RCC\\immune\\macrophage\\rcc.m.rds')


###plot
DefaultAssay(rcc.m)<-'RNA'
m<-c('HLA-DRB1','C1QB','APOE','HSPA1B','CD163','SEPP1','IFITM1','IFITM3','IFI27','FLT1','SPARC','RGS5')
VlnPlot(rcc.m,m)
FeaturePlot(rcc.m,m)

m1<-c('FCGR3A','FCGR2A','FCGR1A','CD68','CD80','CD86','CXCL10','CCR7','IL7R')
FeaturePlot(rcc.m,m1)
VlnPlot(rcc.m,m1)

m2<-c('MRC1','CD163','CCL18','CD209','IL10','ABL2')
VlnPlot(rcc.m,m2)

tn<-c('CCR7','LEF1','SELL','TCF7','CD27','CD28','S1PR1','CD247','CD8A')
mono<-c('FCN1','APOBEC3A','THBS1')
hlahi<-c('C1QA','C1QB','APOE','TREM2')
hlaint<-c('CD163','CSF1R','HSPA1A','HSPA1B','SEPP1','BAG3')
isg<-c('RGS5','IGFBP7','CCL5','LST1','EMCN','AQP1','MT2A')
m<-c('ITGAM','ITGAX','CD80','CD86','MRC1','ARNTL')
m<-c('HLA-DRA','CD14','FCGR3A','S100A8','FCER1A','IRF8','CLEC10A','TPSAB1')

FeaturePlot(rcc.m,M)

####
FeaturePlot(rcc.m,c('CA9','LRP2','CUBN','CD163','CD68','FCGR3A','CD52','S100A8','CLEC9A','NCAM1','CD4','CD8A','PECAM1','FLT1'))

rcc.m$celltype.refined<-factor(rcc.m$celltype.refined,levels = c('DC','Monocyte','M1','M2','Macrophage-Cancer cell','Multilet'),ordered = T)
cg=c('CLEC9A','CD52','S100A8','CD68','FCGR3A','CD163','CA9','LRP2','CUBN','NCAM1','CD247','PECAM1','FLT1')
DotPlot(rcc.m, features = unique(cg),
        assay='RNA',group.by = 'celltype.refined'  )  + coord_flip()+ RotatedAxis() 

FeaturePlot(rcc.m,c('CA9','LRP2','CUBN','CD163','CD68','FCGR3A','CD52','S100A8','CLEC9A'))
VlnPlot(rcc.m,c('LRP2','CUBN','CD163','FCGR3A','MKI67','CD52','FCN1','CLEC9A'),group.by = 'celltype.refined',pt.size = 0)

########
#https://www.sciencedirect.com/science/article/pii/S0092867420316135?via%3Dihub#figs2
#Single-cell landscape of the ecosystem in early-relapse hepatocellular carcinoma 
DefaultAssay(rcc.m)<-'RNA'
library(GSVA)
library(readxl)
mm <- read_excel("E:\\work\\RCC\\pathway\\macro\\1-s2.0-S0092867420316135-mmc3.xlsx", 
                                             skip = 1)
m1<-data.frame(na.omit(mm[,1]))[,1]
m2<-data.frame(na.omit(mm[,2]))[,1]
m3<-data.frame(na.omit(mm[,3]))[,1]
m4<-data.frame(na.omit(mm[,4]))[,1]

gene<-list(m1,m2,m3,m4)
data<-rcc.m@assays$RNA@data
res<-gsva(data,gene,method='ssgsea')
rownames(res)<-colnames(mm)[1:4]
saveRDS(res,'gsva_res_RNA.rds')
#saveRDS(res,'gsva_res.rds')

#plot
m1m2<-subset(rcc.m,celltype.refined %in% c('Macro0','Macro1','Macro2','Macro3'))
res<-data.frame(t(data.frame(res)))
# res$type<-'M2'
# res[colnames(m1),'type']<-'M1'
res$cluster<-rcc.m$celltype.refined

m1score=c()
m2score=c()
proscore=c()
antiscore=c()

for(j in c('Macro0','Macro1','Macro2','Macro3')){
  t=subset(res,cluster==j)
  m1score<-c(m1score,mean(t[,1]))
  m2score<-c(m2score,mean(t[,2]))
  proscore<-c(proscore,mean(t[,3]))
  antiscore<-c(antiscore,mean(t[,4]))
}
df2<-data.frame(m1score=m1score,m2score=m2score,proscore=proscore,antiscore=antiscore)


ggplot(res[gsub('-','.',colnames(m1m2)),],aes(x=M1_Polarization,y=M2_Polarization,color=cluster))+geom_point()+theme_bw()+geom_point(data=df2,aes(x=m1score,y=m2score),color=c('red','green','blue','purple'),size=5,shape=17)
ggplot(res[gsub('-','.',colnames(m1m2)),],aes(y=Anti_inflammatory,x=Pro_inflammatory,color=cluster))+geom_point()+geom_point(data=df2,aes(x=antiscore,y=proscore),color=c('red','green','blue','purple'),size=5,,shape=17)+theme_bw()


#####热图 file:///C:/Users/txw/Downloads/pnas.2103240118.sapp.pdf
library(pheatmap)

M1m<-c('FCGR1A','FCGR2A','FCGR3A','ITGAX','SOCS1','CXCL10')
M2m<-c('MRC1','F13A1','CD163','STAB1','CHID1','TGM2','FN1','MSR1','TGFBI','OLR1')
data<-data.frame(m1m2@assays$RNA@data)[c(M1m,M2m),]
data = sweep(data,1, apply(data,1,median,na.rm=T))
#data<-data.frame(t(scale(t(data))))
#data[data>3]<-3

anno_col<-m1m2@meta.data
anno_col<-anno_col[order(anno_col$seurat_clusters),]
#row<-rownames(anno_col)
anno_col<-anno_col[,c(6,22)]
#rownames(anno_col)<-row
rownames(anno_col)<-gsub('-','\\.',rownames(anno_col))
#colnames(anno_col)<-'cluster'

anno_row<-data.frame(c(rep('M1 markers',6),rep('M2 markers',10)))
rownames(anno_row)<-c(M1m,M2m)
colnames(anno_row)<-'type'

data<-data[,rownames(anno_col)]
pheatmap(data, annotation_col = anno_col, annotation_row = anno_row,cluster_rows = F,cluster_cols = F,
         show_colnames = F)



#################
library(GSVA)
library(clusterProfiler)
## 读入hallmarks gene set
hallmarks <- read.gmt("E:\\work\\RCC\\pathway\\macro\\c2.cp.kegg.v7.5.1.symbols.gmt")
# 需要网络
kegg_list = split(hallmarks$gene, hallmarks$term)
data<-rcc.m@assays$RNA@data
res<-gsva(expr=data,kegg_list,kcdf="Gaussian",method='gsva')
saveRDS(res,'kegg_res.rds')


###########
Idents(m1m2)<-'celltype.refined'
#vlnplot
m1<-c('CD68','FCGR1A','FCGR2A','FCGR3A','ITGAX')
m2<-c('CD163','MRC1','F13A1','STAB1','TGFBI')
VlnPlot(m1m2,m1,pt.size = 0)
VlnPlot(m1m2,m2,pt.size = 0)

macro3<-c('C1QA','C1QB','C1QC','TREM2','APOE')
VlnPlot(m1m2,macro3,pt.size = 0)

mhc<-c('HLA-A','HLA-B','HLA-C','HLA-DMA','HLA-DMB','HLA-DOA','HLA-DPA1','HLA-DPB1','HLA-DQA1','HLA-DQA2','HLA-DQB1',
       'HLA-DRB1')
VlnPlot(m1m2,mhc,pt.size = 0)


###############
###enrichment
library(org.Hs.eg.db)
library(clusterProfiler)
library(patchwork)
library(tidyverse)
ego_ALL<-enrichGO(gene = row.names(y),
                  OrgDb = 'org.Hs.eg.db',
                  keyType = 'SYMBOL',
                  ont = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.05)
#kegg分析
options(clusterProfiler.download.method = "wininet")
genelist<-bitr(rownames(macro0),
               fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = 'org.Hs.eg.db')
genelist <- pull(genelist,ENTREZID)               
ekegg <- enrichKEGG(gene = genelist, organism = 'hsa')


##################
###plot
m1m2<-subset(rcc.m,celltype.refined %in% c('Macro0','Macro1','Macro2','Macro3'))
df<-data.frame(matrix(0,ncol = 4,nrow = 6))
colnames(df)<-c('Macro0','Macro1','Macro2','Macro3')
rownames(df)<-c('ccRCCb','ccRCCc','Normala','Normald','ccRCCa','ccRCCd')
for(i in c('Macro0','Macro1','Macro2','Macro3')){
  t=subset(m1m2,celltype.refined==i)
  df[,i]=as.character(table(t$orig.ident))
}

data<-data.frame(celltype=c(rep('Macro0',6),rep('Macro1',6),rep('Macro2',6),rep('Macro3',6)),
                 sample=rep(c('ccRCCa','ccRCCb','ccRCCc','ccRCCd','Normala','Normald'),4),
                 value=rep(0,24))
for(i in 1:24){
  data[i,3]=df[data[i,2],data[i,1]]
}             
data[,3]<-as.numeric(data[,3])
ggplot( data, aes( x = celltype, y = value, fill = sample))+
  geom_bar(stat="identity",position="stack", color="black", width=0.7,size=0.25)

#stage
data$stage<-rep(c('T1a','T3a','T1a','T3a','Normal','Normal'),4)
ggplot( data, aes( x = celltype, y = value, fill = stage))+
  geom_bar(stat="identity",position="fill", color="black", width=0.7,size=0.25)


################################
#rcc.m delete
library(Hmisc)
cellname<-readRDS("E:\\work\\RCC\\data\\cellname.rds")
rcc.m<-rcc.m[,cellname]
rcc.m<-subset(rcc.m,celltype.refined %nin% c('Macrophage-CD247+ cell','Macrophage-EPC'))

DefaultAssay(rcc.m)<-'integrated'
rcc.m<-RunPCA(rcc.m,features = VariableFeatures(rcc.m))
rcc.m<-RunUMAP(rcc.m,dims = 1:30)
saveRDS(rcc.m,'rcc.m_2840.rds')




