
setwd('E:\\work\\RCC\\data\\normald')

#######
normald<-readRDS('normald_1.rds')
library(DoubletFinder)
#����һ��������Ѳ����Ĺ��̣������ٶ���
sweep.res.list_normald <- paramSweep_v3(normald, PCs = 1:30, sct = FALSE)
#head(sweep.res.list_normald)
sweep.stats_normald <- summarizeSweep(sweep.res.list_normald, GT = FALSE)
bcmvn_normald <- find.pK(sweep.stats_normald) #���Կ�����Ѳ����ĵ�
## ������ѵĲ����ǣ�
mpK<-as.numeric(as.vector(bcmvn_normald$pK[which.max(bcmvn_normald$BCmetric)]))

##Homotypic Doublet Proportion Estimate -------------------------------------
annotations <- normald@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)  
#DoubletRate = ncol(normald)*8*1e-6 #��ÿ����1000��ϸ����˫ϸ����������ǧ��֮8������
DoubletRate = 0.2
#����ͬԴ˫ϸ������������modelHomotypic()�еĲ�����Ϊ���˫ϸ���������Ǵ�seurat_clusters������˫ϸ�� 
nExp_poi <- round(DoubletRate*length(normald$seurat_clusters))  #����ṩcelltype��������seurat_clusters��
# ����˫ϸ������
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## ���ʹ��ȷ���õĲ�������Doublets. Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
#normald <- doubletFinder_v3(normald, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
normald <- doubletFinder_v3(normald, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)
# ʹ��nExp = nExp_poi��nExp = nExp_poi.adj,�ֱ����doublets�������Ա����ȷ����Щϸ����Doublet-High Confidience
DimPlot(normald, reduction = "umap", group.by ="DF.classifications_0.25_0.005_3657")

saveRDS(normald,'normald_2.rds')