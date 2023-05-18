
setwd('E:\\work\\RCC\\data\\hsma')

#######
hsma<-readRDS('hsma_1.rds')
library(DoubletFinder)
#����һ��������Ѳ����Ĺ��̣������ٶ���
sweep.res.list_hsma <- paramSweep_v3(hsma, PCs = 1:30, sct = FALSE)
#head(sweep.res.list_hsma)
sweep.stats_hsma <- summarizeSweep(sweep.res.list_hsma, GT = FALSE)
bcmvn_hsma <- find.pK(sweep.stats_hsma) #���Կ�����Ѳ����ĵ�
## ������ѵĲ����ǣ�
mpK<-as.numeric(as.vector(bcmvn_hsma$pK[which.max(bcmvn_hsma$BCmetric)]))

##Homotypic Doublet Proportion Estimate -------------------------------------
annotations <- hsma@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)  
#DoubletRate = ncol(hsma)*8*1e-6 #��ÿ����1000��ϸ����˫ϸ����������ǧ��֮8������
DoubletRate = 0.2
#����ͬԴ˫ϸ������������modelHomotypic()�еĲ�����Ϊ���˫ϸ���������Ǵ�seurat_clusters������˫ϸ�� 
nExp_poi <- round(DoubletRate*length(hsma$seurat_clusters))  #����ṩcelltype��������seurat_clusters��
# ����˫ϸ������
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## ���ʹ��ȷ���õĲ�������Doublets. Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
#hsma <- doubletFinder_v3(hsma, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
hsma <- doubletFinder_v3(hsma, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)
# ʹ��nExp = nExp_poi��nExp = nExp_poi.adj,�ֱ����doublets�������Ա����ȷ����Щϸ����Doublet-High Confidience
DimPlot(hsma, reduction = "umap", group.by ="DF.classifications_0.25_0.005_1488")

saveRDS(hsma,'hsma_2.rds')