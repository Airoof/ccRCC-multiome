
setwd('E:\\work\\RCC\\data\\hsma')

set.seed(1)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)
library(BiocGenerics)

counts_t <- Read10X_h5("outs/filtered_feature_bc_matrix.h5")
# create a Seurat object containing the RNA adata
hsma <- CreateSeuratObject(
  counts = counts_t$`Gene Expression`,
  assay = "RNA",project = "hsma",min.cells=3,min.features=200
)
#使用MT-作为线粒体基因集
hsma[["percent.mt"]]<-PercentageFeatureSet(hsma,pattern = "^MT-")
#qc指标可视化
VlnPlot(hsma,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
hsma<-subset(hsma,subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt<15 & nCount_RNA>1000 & nCount_RNA<50000)

#数据归一化   
hsma<-NormalizeData(hsma)
#高变异基因选择
hsma<-FindVariableFeatures(hsma,selection.method = "vst",nfeatures = 2000)

#缩放数据
all.genes<-rownames(hsma)
hsma<-ScaleData(hsma,features = all.genes)
#降维
hsma<-RunPCA(hsma,features = VariableFeatures(object = hsma))
ElbowPlot(hsma)
#计算最邻近距离
hsma<-FindNeighbors(hsma,dims=1:20)
#聚类，包含设置下游聚类的“间隔尺度”的分辨率参数resolution，增加值会导致更多的聚类
hsma<-FindClusters(hsma,resolution = 0.5)
hsma<-RunUMAP(hsma,dims=1:20) #saveRDS(hsma,"hsma.rds")
DimPlot(hsma,reduction = "umap",label = T)

##注释
#0  cancer cell ca9
#1  cancer cell ca9
#2  cancer cell ca9
#3  cancer cell ca9
#4  cancer cell ca9
#5  endothelial kdr/pecam1/PLVAP
#6  macrophage  cd68/cd163 
#7  T cell
#8  Mesangial cell  PDGFRB
#9  cancer cell ca9
#10 cancer cell ca9
#11 mast  kit/slc18a2/tpsb2


type<-c('0'='cancer cell','1'='cancer cell','2'='cancer cell','3'='cancer cell','4'='cancer cell','5'='endothelial','6'='macrophage','7'='T cell','8'='mesangial cell','9'=
          'cancer cell','10'='cancer cell','11'='mast')
hsma[['celltype']] = unname(type[hsma@meta.data$seurat_clusters])
markers<-FindAllMarkers(hsma)
saveRDS(markers,'markers.rds')
saveRDS(hsma,'hsma.rds')
################细胞注释###########
library(SingleR)
library(tibble)
hpca.se<-HumanPrimaryCellAtlasData()
bpe.se<-BlueprintEncodeData()
dic.se<-DatabaseImmuneCellExpressionData()
obj<-hsma
obj@meta.data$cell.type<-Idents(obj)
test_1<-as.SingleCellExperiment(obj)
anno<-SingleR(test=test_1,ref = list(HP=hpca.se,BP=bpe.se),labels = list(hpca.se$label.fine,bpe.se$label.fine),method = "cluster",cluster=test_1$cell.type)
anno$cluster<-rownames(anno)
fin<-anno %>% dplyr::tbl_df() %>% dplyr::select(cluster,labels)
new.cluster.ids<-fin$labels
names(new.cluster.ids)<-levels(obj)
obj<-RenameIdents(obj,new.cluster.ids)
DimPlot(obj,reduction = "umap",label = TRUE,pt.size = 0.5)+NoLegend()
