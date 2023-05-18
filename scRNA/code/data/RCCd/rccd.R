setwd('E:\\work\\RCC\\data\\RCCd')

set.seed(1)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)
library(BiocGenerics)

counts_t <- Read10X_h5("outs/filtered_feature_bc_matrix.h5")
# create a Seurat object containing the RNA adata
rccd <- CreateSeuratObject(
  counts = counts_t$`Gene Expression`,
  assay = "RNA",project = "RCCd",min.cells=3,min.features=200
)
#ʹ��MT-��Ϊ���������
rccd[["percent.mt"]]<-PercentageFeatureSet(rccd,pattern = "^MT-")
#qcָ����ӻ�
VlnPlot(rccd,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
rccd<-subset(rccd,subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt<15 & nCount_RNA>1000 & nCount_RNA<50000)

#���ݹ�һ��   
rccd<-NormalizeData(rccd)
#�߱������ѡ��
rccd<-FindVariableFeatures(rccd,selection.method = "vst",nfeatures = 2000)

#��������
all.genes<-rownames(rccd)
rccd<-ScaleData(rccd,features = all.genes)
#��ά
rccd<-RunPCA(rccd,features = VariableFeatures(object = rccd))
ElbowPlot(rccd)
#�������ڽ�����
rccd<-FindNeighbors(rccd,dims=1:20)
#���࣬�����������ξ���ġ�����߶ȡ��ķֱ��ʲ���resolution������ֵ�ᵼ�¸���ľ���
rccd<-FindClusters(rccd,resolution = 0.5)
rccd<-RunUMAP(rccd,dims=1:20) #saveRDS(rccd,"rccd.rds")
DimPlot(rccd,reduction = "umap",label = T)

##ע��
#0  cancer cell ca9
#1  cancer cell ca9
#2  cancer cell ca9
#3  cancer cell ca9
#4  cancer cell ca9
#5  macrophage cd68/cd163/hla-dra
#6  macrophage cd68/cd163/hla-dra
#7  cancer cell ca9
#8  Mesangial cell  PDGFRB
#9  endothelial kdr/pecam1/PLVAP
#10 mast  kit/slc18a2/tpsb2
#11 CD8T cd3d/cd52/trac 



type<-c('0'='cancer cell','1'='cancer cell','2'='cancer cell','3'='cancer cell','4'='cancer cell','5'='macrophage','6'='macrophage','7'='cancer cell','8'='mesangial cell','9'=
          'endothelial','10'='mast','11'='CD8T')
rccd[['celltype']] = unname(type[rccd@meta.data$seurat_clusters])
markers<-FindAllMarkers(rccd)
saveRDS(markers,'markers.rds')
saveRDS(rccd,'rccd.rds')

################ϸ��ע��###########
library(SingleR)
library(tibble)
hpca.se<-HumanPrimaryCellAtlasData()
bpe.se<-BlueprintEncodeData()
dic.se<-DatabaseImmuneCellExpressionData()
obj<-rccd
obj@meta.data$cell.type<-Idents(obj)
test_1<-as.SingleCellExperiment(obj)
anno<-SingleR(test=test_1,ref = list(HP=hpca.se,BP=bpe.se),labels = list(hpca.se$label.fine,bpe.se$label.fine),method = "cluster",cluster=test_1$cell.type)
anno$cluster<-rownames(anno)
fin<-anno %>% dplyr::tbl_df() %>% dplyr::select(cluster,labels)
new.cluster.ids<-fin$labels
names(new.cluster.ids)<-levels(obj)
obj<-RenameIdents(obj,new.cluster.ids)
DimPlot(obj,reduction = "umap",label = TRUE,pt.size = 0.5)+NoLegend()