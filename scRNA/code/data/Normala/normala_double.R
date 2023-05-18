
setwd('E:\\work\\RCC\\data\\normala')

#######
normala<-readRDS('normala_1.rds')
library(DoubletFinder)
#����һ��������Ѳ����Ĺ��̣������ٶ���
sweep.res.list_normala <- paramSweep_v3(normala, PCs = 1:30, sct = FALSE)
#head(sweep.res.list_normala)
sweep.stats_normala <- summarizeSweep(sweep.res.list_normala, GT = FALSE)
bcmvn_normala <- find.pK(sweep.stats_normala) #���Կ�����Ѳ����ĵ�
## ������ѵĲ����ǣ�
mpK<-as.numeric(as.vector(bcmvn_normala$pK[which.max(bcmvn_normala$BCmetric)]))

##Homotypic Doublet Proportion Estimate -------------------------------------
annotations <- normala@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)  
#DoubletRate = ncol(normala)*8*1e-6 #��ÿ����1000��ϸ����˫ϸ����������ǧ��֮8������
DoubletRate = 0.2
#����ͬԴ˫ϸ������������modelHomotypic()�еĲ�����Ϊ���˫ϸ���������Ǵ�seurat_clusters������˫ϸ�� 
nExp_poi <- round(DoubletRate*length(normala$seurat_clusters))  #����ṩcelltype��������seurat_clusters��
# ����˫ϸ������
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## ���ʹ��ȷ���õĲ�������Doublets. Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
#normala <- doubletFinder_v3(normala, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
normala <- doubletFinder_v3(normala, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)
# ʹ��nExp = nExp_poi��nExp = nExp_poi.adj,�ֱ����doublets�������Ա����ȷ����Щϸ����Doublet-High Confidience
DimPlot(normala, reduction = "umap", group.by ="DF.classifications_0.25_0.005_3142")

saveRDS(normala,'normala_2.rds')