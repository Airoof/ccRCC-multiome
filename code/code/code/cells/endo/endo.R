
setwd('E:\\work\\RCC')

library(Seurat)
rcc.combined<-readRDS('data/combined/rcc_combined.rds')
endo<-subset(rcc.combined,celltype.refined=='ENDO')

DefaultAssay(endo)<-'integrated'
DimPlot(endo,label = T)

endo<-FindVariableFeatures(endo,nFeatures=2000)
endo <- ScaleData(endo, verbose = FALSE)
endo <- RunPCA(endo, npcs = 30, verbose = FALSE) 
endo <- RunUMAP(endo, reduction = "pca", dims = 1:15) 
endo <- FindNeighbors(endo, reduction = "pca", dims = 1:15)
endo <- FindClusters(endo, resolution =1 )
DimPlot(endo,label = T)
endo$tissue<-ifelse(endo$orig.ident %in% c('normala','normald'),'normal','tumor')
DimPlot(endo,group.by = 'tissue')

DefaultAssay(endo)<-'RNA'
##https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8214680/
gene_pnas<-c('CLDN5','AQP1','PLVAP','PDGFB','RGS5','ITGA8','POSTN')
FeaturePlot(endo,gene_pnas)
VlnPlot(endo,gene_pnas,pt.size = 0)

##https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8969307/
gene_gb<-c('SOST','NDUFA4L2','GJA5','ACKR1','CXCR4')
FeaturePlot(endo,gene_gb)
VlnPlot(endo,gene_gb,pt.size = 0)

VlnPlot(endo,c('nFeature_RNA','nCount_RNA'))
####
#0 PLVAP peritubular capillary endothelial cells PCEC
#1 PLVAP peritubular capillary endothelial cells PCEC
#2 sost/crhbp/ehd3 GCEC glomerular  capillary endothelial cells
#3 PLVAP peritubular capillary endothelial cells PCEC
#4 ackr1/selp venular endothelial cells VEC
#5 ENDO-EPC
#6 ENDO-EPC
#7 ENDO-EPC
#8 ENDO-PTPRC+ cell
#9 cldn5/SERPINE2/sox17  AEC arteriolar endothelial cells 
#10 ENDO-PTPRC+ cell
#11 ENDO(proliferating)

type<-c('0'='PCEC','1'='PCEC','2'='GCEC','3'='PCEC','4'='VEC',
        '5'='ENDO-EPC','6'='ENDO-EPC','7'='ENDO-EPC','8'='ENDO-PTPRC_cell','9'=
          'AEC','10'='ENDO-PTPRC_cell','11'='ENDO(proliferating)')
endo[['celltype']] = unname(type[endo@meta.data$seurat_clusters])
saveRDS(endo,'cells/endo/endo.rds')
