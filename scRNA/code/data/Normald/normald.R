
setwd('E:\\work\\RCC\\data\\normald')

set.seed(1)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)
library(BiocGenerics)

counts_t <- Read10X_h5("outs/filtered_feature_bc_matrix.h5")
# create a Seurat object containing the RNA adata
normald <- CreateSeuratObject(
  counts = counts_t$`Gene Expression`,
  assay = "RNA",project = "normald",min.cells=3,min.features=200
)
#使用MT-作为线粒体基因集
normald[["percent.mt"]]<-PercentageFeatureSet(normald,pattern = "^MT-")
#qc指标可视化
VlnPlot(normald,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
normald<-subset(normald,subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt<15 & nCount_RNA>1000 & nCount_RNA<50000)

#数据归一化   
normald<-NormalizeData(normald)
#高变异基因选择
normald<-FindVariableFeatures(normald,selection.method = "vst",nfeatures = 2000)

#缩放数据
all.genes<-rownames(normald)
normald<-ScaleData(normald,features = all.genes)
#降维
normald<-RunPCA(normald,features = VariableFeatures(object = normald))
ElbowPlot(normald)
#计算最邻近距离
normald<-FindNeighbors(normald,dims=1:20)
#聚类，包含设置下游聚类的“间隔尺度”的分辨率参数resolution，增加值会导致更多的聚类
normald<-FindClusters(normald,resolution = 0.5)
normald<-RunUMAP(normald,dims=1:20) #saveRDS(normald,"normald.rds")
DimPlot(normald,reduction = "umap",label = T)

FeaturePlot(normald,c('SLC44A5','VCAM1'))

##注释
#0  PST proximal convoluted tubule CUBN/LRP2/SLC34A1/SLC5A12/ALDOB SLC22A8	Proximal convoluted tubule
#1  PCT proximal convoluted tubule CUBN/LRP2/SLC34A1/SLC5A12/ALDOB SLC22A7	Proximal straight tubule
#2  LOH loop of Henle  slc12a1
#3  DCT distal convoluted tubule  slc12a3
#4  DCT-CT distal convoluted tubule connecting tubule slc8a1
#5  LOH-PCT
#6  LOH loop of Henle  slc12a1 CLCNKA 升肢细段
#7  PST proximal convoluted tubule CUBN/LRP2/SLC34A1/SLC5A12/ALDOB  SLC22A8	Proximal convoluted tubule
#8  LOH-PCT
#9  ENDO endothelial kdr/pecam1/PLVAP/flt1
#10 CD-PC collecting duct principal cell  aqp2
#11 CD-ICA collecting duct intercalated cell A  aqp6/atp6v0d2
#12 SLC44A5+VCAM1+ EPI
#13 DCT distal convoluted tubule  slc12a3 
#14 MES mesangial cell pdgfrb  
#15 PEC parietal epithelial cells  cfh
#16 DC-ICB collecting duct intercalated cell B  slc26a4/atp6v0d2
#17 DCT-CT distal convoluted tubule connecting tubule slc8a1
#18 PODO podocyte  nphs1/nphs2
#19 LEU ptprc

type<-c('0'='PST','1'='PCT','2'='LOH','3'='DCT','4'='DCT-CT',
        '5'='LOH-PT','6'='LOH','7'='PST','8'='LOH-PT','9'=
          'ENDO','10'='CD-PC','11'='DC-ICA','12'=
          'SLC44A5+VCAM1+ EPI','13'='DCT','14'='MES','15'=
          'PEC','16'='DC-ICB','17'='DCT-CT','18'=
          'PODO','19'='LEU')
normald[['celltype']] = unname(type[normald@meta.data$seurat_clusters])
markers<-FindAllMarkers(normald)
saveRDS(markers,'markers.rds')
saveRDS(normald,'normald.rds')

####
cg1=c('SLC34A1','LRP2','HAVCR1','CFH','SLC12A1','SLC12A3','SLC8A1','AQP2','SLC26A7','SLC26A4','NPHS2','EMCN','PIEZO2','COL1A2','PTPRC')
p1 <- DotPlot(normald, features = unique(cg1),
              assay='RNA'  )  + coord_flip()
p1

cg2=c('CUBN','LRP2','SLC34A1','SLC5A12','SLC5A2','ALDOB','CFH','SLC12A1','SLC12A3','SLC12A2','SLC8A1','AQP2','AQP6','SLC26A4','ATP6V0D2','NPHS1','NPHS2','PECAM1','FLT1','PDGFRB','PTPRC')
p2 <- DotPlot(normald, features = unique(cg2),
              assay='RNA'  )  + coord_flip()
p2 

#######
normald$celltype<-factor(normald$celltype,levels = c('CD-PC','DC-ICA','DC-ICB','DCT','DCT-CT','LOH','LOH-PT','PCT','PST',
                                                      'SLC44A5+VCAM1+ EPI','ENDO','MES','PEC','PODO','LEU'),ordered = T)
cg=c('AQP2','AQP6','ATP6V0D2','SLC26A4','SLC12A3','SLC12A2','SLC8A1','SLC12A1','CUBN','LRP2','SLC34A1','SLC5A12','ALDOB','SLC22A8','SLC22A7',
     'SLC44A5','VCAM1','PECAM1','FLT1','EMCN','PDGFRB','CFH','NPHS1','NPHS2','PTPRC')
DotPlot(normald, features = unique(cg),
        assay='RNA',group.by = 'celltype'  )  + coord_flip()+ RotatedAxis() 
################细胞注释###########
library(SingleR)
library(tibble)
hpca.se<-HumanPrimaryCellAtlasData()
bpe.se<-BlueprintEncodeData()
dic.se<-DatabaseImmuneCellExpressionData()
obj<-normald
obj@meta.data$cell.type<-Idents(obj)
test_1<-as.SingleCellExperiment(obj)
anno<-SingleR(test=test_1,ref = list(HP=hpca.se,BP=bpe.se),labels = list(hpca.se$label.fine,bpe.se$label.fine),method = "cluster",cluster=test_1$cell.type)
anno$cluster<-rownames(anno)
fin<-anno %>% dplyr::tbl_df() %>% dplyr::select(cluster,labels)
new.cluster.ids<-fin$labels
names(new.cluster.ids)<-levels(obj)
obj<-RenameIdents(obj,new.cluster.ids)
DimPlot(obj,reduction = "umap",label = TRUE,pt.size = 0.5)+NoLegend()

