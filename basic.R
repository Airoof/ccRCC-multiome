library(Seurat)
library(data.table)
data <- Read10X(data.dir = "./")

pbmc <- CreateSeuratObject(counts = data, project = "sample1",min.cells = 3)

pbmc2 <- pbmc

pbmc <- merge(pbmc, pbmc2)

pbmc$percent.mt <- PercentageFeatureSet(pbmc, pattern = "^MT-")

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
singlets <- as.vector(as.matrix(fread('singlets.txt',header = F)))
wt <- pbmc[,singlets]

pbmc <- subset(pbmc, subset = nFeature_RNA >= 500 & nCount_RNA <= 100000 &
                 percent.mt <= 20 )

# Matrix::writeMM(pbmc@assays$RNA@counts,'matrix.mtx')
# cat(colnames(pbmc),file='barcodes.tsv',sep = '\n')
# library(scater)
# library(SingleCellExperiment)
# library(SingleR)
# hpca.se <- HumanPrimaryCellAtlasData()
# expe <- SingleCellExperiment(pbmc@assays$RNA@counts)
# expe@assays@data@listData$counts <- as.matrix(expe@assays@data@listData[[1]])
# 
# expe <- logNormCounts(expe)
# 
# pred <- SingleR(test = expe, ref = hpca.se, labels = hpca.se$label.main)
# pbmc[['slabel']] <- pred$first.labels

pbmc <- SCTransform(pbmc,return.only.var.genes = F)

pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

pbmc <- ScaleData(pbmc)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
DimPlot(pbmc)
ElbowPlot(pbmc,50)
dim_num <- 20
pbmc <- FindNeighbors(pbmc, dims = 1:dim_num)
pbmc <- FindClusters(pbmc, resolution = 0.6)
pbmc <- RunUMAP(pbmc,dims = 1:dim_num)
pbmc <- RunUMAP(pbmc,dims = 1:dim_num, umap.method = 'uwot')
pbmc <- RunTSNE(pbmc,dims = 1:dim_num)
DimPlot(pbmc, reduction = "umap") 
DimPlot(pbmc, reduction = "umap",label = T) 
markers <- FindMarkers(pbmc, ident.1 = c(2, 5))
DimPlot(pbmc, reduction = "pca",label = T)
DimPlot(pbmc, reduction = "umap",label = T,group.by = 'orig.ident')
DimPlot(pbmc, reduction = "umap",group.by = 'celltype')
rat1neuron <- pbmc
FeaturePlot(rat1neuron, c("Kcnn2"))
FeaturePlot(rat2neuron, c("Kcnn2"))
write.csv(markers, 'neuron_markers.csv')
pbmc <- pbmc[[pbmc$DF.classifications_0.25_0.25_125 == 'Singlet']]
DimPlot(pbmcorig, reduction = "umap",label = T,group.by = 'orig.ident')
DimPlot(pbmc, reduction = "tsne",label = T)
DimPlot(pbmc, reduction = "umap",label = T,group.by = 'celltype')
DimPlot(pbmc, reduction = "umap",label = T,group.by = 'panglaotype')
DimPlot(pbmc, reduction = "tsne",label = T,group.by = 'celltype')
FeaturePlot(pbmc,features = c('CD4'))
saveRDS(pbmc, 'gastric_1.rds')
pbmc.markers <- FindAllMarkers(pbmc)
pbmcsub <- subset(pbmc,subset = orig.ident == 'ko')

pbmc[['celltype']] <- 'None'

VlnPlot(pbmc, features = c("CD14"))

library(magrittr)
library(dplyr)
topn <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(pbmc, features = topn$gene) + NoLegend()