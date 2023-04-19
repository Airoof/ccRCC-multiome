
setwd('E:\\work\\RCC\\data\\hsma')

#######
hsma<-readRDS('hsma_1.rds')
library(DoubletFinder)
#这是一个测试最佳参数的过程，运行速度慢
sweep.res.list_hsma <- paramSweep_v3(hsma, PCs = 1:30, sct = FALSE)
#head(sweep.res.list_hsma)
sweep.stats_hsma <- summarizeSweep(sweep.res.list_hsma, GT = FALSE)
bcmvn_hsma <- find.pK(sweep.stats_hsma) #可以看到最佳参数的点
## 所以最佳的参数是：
mpK<-as.numeric(as.vector(bcmvn_hsma$pK[which.max(bcmvn_hsma$BCmetric)]))

##Homotypic Doublet Proportion Estimate -------------------------------------
annotations <- hsma@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)  
#DoubletRate = ncol(hsma)*8*1e-6 #按每增加1000个细胞，双细胞比率增加千分之8来计算
DoubletRate = 0.2
#估计同源双细胞比例，根据modelHomotypic()中的参数人为混合双细胞。这里是从seurat_clusters中来混双细胞 
nExp_poi <- round(DoubletRate*length(hsma$seurat_clusters))  #最好提供celltype，而不是seurat_clusters。
# 计算双细胞比例
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## 最后，使用确定好的参数鉴定Doublets. Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
#hsma <- doubletFinder_v3(hsma, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
hsma <- doubletFinder_v3(hsma, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)
# 使用nExp = nExp_poi和nExp = nExp_poi.adj,分别进行doublets鉴定，以便后续确定哪些细胞是Doublet-High Confidience
DimPlot(hsma, reduction = "umap", group.by ="DF.classifications_0.25_0.005_1488")

saveRDS(hsma,'hsma_2.rds')
