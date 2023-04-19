
setwd('E:\\work\\RCC')

library(Seurat)
library(Signac)
set.seed(1)

#read atac data
atac<-readRDS('/public/home/xwtang/work/RCC/atac/data/multi_afterMACS2_intersected.rds')
rcca<-readRDS('/public/home/xwtang/work/RCC/data/rcca/rcca_2.rds')
rcca<-RenameCells(rcca,add.cell.id='RCCa')
rcca<-rcca[,colnames(atac)]
rccd<-readRDS('/public/home/xwtang/work/RCC/data/rccd/rccd_2.rds')
rccd<-RenameCells(rccd,add.cell.id='RCCd')
rccd<-rccd[,colnames(atac)]
hsma<-readRDS('/public/home/xwtang/work/RCC/data/hsma/hsma_2.rds')
hsma<-RenameCells(hsma,add.cell.id='HSMa')
hsma<-hsma[,colnames(atac)]
hsmb<-readRDS('/public/home/xwtang/work/RCC/data/hsmb/hsmb_2.rds')
hsmb<-RenameCells(hsmb,add.cell.id='HSMb')
hsmb<-hsmb[,colnames(atac)]
normala<-readRDS('/public/home/xwtang/work/RCC/data/normala/normala_2.rds')
normala<-RenameCells(normala,add.cell.id='Normala')
normala<-normala[,colnames(atac)]
normald<-readRDS('/public/home/xwtang/work/RCC/data/normald/normald_2.rds')
normald<-RenameCells(normald,add.cell.id='Normald')
normald<-normald[,colnames(atac)]

sclist<-list(rcca,rccd,hsma,hsmb,normala,normald)
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = sclist)
rcc.anchors <- FindIntegrationAnchors(object.list = sclist, anchor.features = features)
# this command creates an 'integrated' data assay
rcc.combined <- IntegrateData(anchorset = rcc.anchors)
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(rcc.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
rcc.combined <- ScaleData(rcc.combined, verbose = FALSE)
rcc.combined <- RunPCA(rcc.combined, npcs = 30, verbose = FALSE)
rcc.combined <- RunUMAP(rcc.combined, reduction = "pca", dims = 1:30)
rcc.combined <- FindNeighbors(rcc.combined, reduction = "pca", dims = 1:30)
rcc.combined <- FindClusters(rcc.combined, resolution = 1)
saveRDS(rcc.combined,'data/combined/rcc_combined.rds')

###########
library(Seurat)
library(ggplot2)
rcc.combined<-readRDS('data/combined/rcc_combined.rds')
##
meta1<-readRDS('combined/tumor_combined_2_meta.rds')
meta2<-readRDS('combined/normal_combined_2_meta.rds')
rcc.combined$celltype<-c(meta1$celltype,meta2$celltype)
type<-c('0'='EPC','1'='EPC','2'='EPC','3'='EPC','4'='EPC',
        '5'='EPC','6'='Myeloid cell','7'='LOH','8'='EPC','9'=
          'ENDO-PT','10'='ENDO','11'='LOH','12'=
          'LOH-PT','13'='Lymphocyte','14'='CNT','15'=
          'DCT','16'='LOH','17'='EPC','18'=
          'MES','19'='EPC','20'='ICA','21'='ICB','22'='Lymphocyte','23'='ICA','24'='PODO','25'='ENDO','26'='Mast','27'='Plasma cell')
rcc.combined[['celltype.refined']] = unname(type[rcc.combined@meta.data$seurat_clusters])
saveRDS(rcc.combined,'data/combined/rcc_combined_2.rds')
#plot
DimPlot(rcc.combined,label = T)
cg=c('PTPRC','CD247','UBE2C','TOP2A','MS4A1','JCHAIN','CD163','CD68','FCGR3A','KIT','SLC18A2','PECAM1','FLT1','EMCN','CFH','PDGFRB','NPHS1','NPHS2','AQP6','ATP6V0D2','SLC26A4','SLC12A3','SLC8A1','SLC12A1','CUBN','LRP2','SLC22A7','SLC22A8','CA9')
DotPlot(rcc.combined, features = unique(cg),
        assay='RNA',group.by = 'celltype.refined'  )  + coord_flip() + RotatedAxis() 

