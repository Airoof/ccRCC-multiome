
setwd('E:\\work\\RCC\\data\\hsmb')

set.seed(1)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)
library(BiocGenerics)

counts_t <- Read10X_h5("outs/filtered_feature_bc_matrix.h5")
# create a Seurat object containing the RNA adata
hsmb <- CreateSeuratObject(
  counts = counts_t$`Gene Expression`,
  assay = "RNA",project = "hsmb",min.cells=3,min.features=200
)
#使用MT-作为线粒体基因集
hsmb[["percent.mt"]]<-PercentageFeatureSet(hsmb,pattern = "^MT-")
#qc指标可视化
VlnPlot(hsmb,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
hsmb<-subset(hsmb,subset =nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt<15 & nCount_RNA>1000 & nCount_RNA<50000)

#数据归一化   
hsmb<-NormalizeData(hsmb)
#高变异基因选择
hsmb<-FindVariableFeatures(hsmb,selection.method = "vst",nfeatures = 2000)

#缩放数据
all.genes<-rownames(hsmb)
hsmb<-ScaleData(hsmb,features = all.genes)
#降维
hsmb<-RunPCA(hsmb,features = VariableFeatures(object = hsmb))
ElbowPlot(hsmb)
#计算最邻近距离
hsmb<-FindNeighbors(hsmb,dims=1:20)
#聚类，包含设置下游聚类的“间隔尺度”的分辨率参数resolution，增加值会导致更多的聚类
hsmb<-FindClusters(hsmb,resolution = 0.5)
hsmb<-RunUMAP(hsmb,dims=1:20) #saveRDS(hsmb,"hsmb.rds")
DimPlot(hsmb,reduction = "umap",label = T)

##注释
#0  cancer cell ca9
#1  cancer cell ca9
#2  
#3  cancer cell ca9
#4  macrophage  cd68/cd163 
#5  T cell cd247
#6  endothelial kdr/pecam1/PLVAP
#7  cancer cell ca9
#8  cancer cell ca9
#9  CD8T  cd8a
#10 Mesangial cell  PDGFRB
#11 cancer cell ca9 (proliferation)

cg=c('CA9','SLC31A2','CD68','CD163','CD247','CD8A','CD4','GNLY','GZMB','KIT','SLC18A2','PTPRC','PECAM1','PDGFRB','NPHS1','NPHS2','EMCN','CLEC4C','JCHAIN','IRF8','UBE2C','TOP2A',
     'CTSS','CD19','MS4A1','FCGR3A','NCAM1','KDR','FLT1','VCAM1','CFH','NPHS1','NPHS2')
DotPlot(hsmb, features = unique(cg),
        assay='RNA'  )  + coord_flip()+ RotatedAxis() 

type<-c('0'='cancer cell','1'='cancer cell','2'='cancer cell','3'='cancer cell','4'='macrophage','5'='T cell','6'='endothelial',
        '7'='cancer cell','8'='cancer cell','9'='T cell','10'='mesangial cell','11'='cancer cell')
hsmb[['celltype']] = unname(type[hsmb@meta.data$seurat_clusters])
markers<-FindAllMarkers(hsmb)
saveRDS(markers,'markers.rds')
saveRDS(hsmb,'hsmb.rds')
################细胞注释###########
library(SingleR)
library(tibble)
hpca.se<-HumanPrimaryCellAtlasData()
bpe.se<-BlueprintEncodeData()
#dic.se<-DatabaseImmuneCellExpressionData()
obj<-hsmb
obj@meta.data$cell.type<-Idents(obj)
test_1<-as.SingleCellExperiment(obj)
anno<-SingleR(test=test_1,ref = list(HP=hpca.se,BP=bpe.se),labels = list(hpca.se$label.fine,bpe.se$label.fine),method = "cluster",cluster=test_1$cell.type)
anno$cluster<-rownames(anno)
fin<-anno %>% dplyr::tbl_df() %>% dplyr::select(cluster,labels)
new.cluster.ids<-fin$labels
names(new.cluster.ids)<-levels(obj)
obj<-RenameIdents(obj,new.cluster.ids)
DimPlot(obj,reduction = "umap",label = TRUE,pt.size = 0.5)+NoLegend()


###cluster11
####富集分析
marker11<-subset(markers,cluster==11)
colnames(marker11)[7]<-'SYMBOL'
#gsea
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

df_id<-bitr(marker11$SYMBOL, #转换的列是df数据框中的SYMBOL列
            fromType = "SYMBOL",#需要转换ID类型
            toType = "ENTREZID",#转换成的ID类型
            OrgDb = "org.Hs.eg.db")#对应的物种，小鼠的是org.Mm.eg.db

df_all<-merge(marker11,df_id,by="SYMBOL",all=F)#使用merge合并

df_all_sort <- df_all[order(df_all$avg_log2FC, decreasing = T),]#先按照logFC降序排序
gene_fc = df_all_sort$avg_log2FC #把foldchange按照从大到小提取出来
names(gene_fc) <- df_all_sort$ENTREZID #给上面提取的foldchange加上对应上ENTREZID

KEGG <- gseKEGG(gene_fc, organism = "hsa") #具体参数在下面

GO <- gseGO(
  gene_fc, #gene_fc
  ont = "BP",# "BP"、"MF"和"CC"或"ALL"
  OrgDb = org.Hs.eg.db,#人类注释基因
  keyType = "ENTREZID",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",#p值校正方法
)
