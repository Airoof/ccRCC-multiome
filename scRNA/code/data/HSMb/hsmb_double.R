
setwd('E:\\work\\RCC\\data\\hsmb')

#######
hsmb<-readRDS('hsmb_1.rds')
library(DoubletFinder)
#����һ��������Ѳ����Ĺ��̣������ٶ���
sweep.res.list_hsmb <- paramSweep_v3(hsmb, PCs = 1:30, sct = FALSE)
#head(sweep.res.list_hsmb)
sweep.stats_hsmb <- summarizeSweep(sweep.res.list_hsmb, GT = FALSE)
bcmvn_hsmb <- find.pK(sweep.stats_hsmb) #���Կ�����Ѳ����ĵ�
## ������ѵĲ����ǣ�
mpK<-as.numeric(as.vector(bcmvn_hsmb$pK[which.max(bcmvn_hsmb$BCmetric)]))

##Homotypic Doublet Proportion Estimate -------------------------------------
annotations <- hsmb@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)  
#DoubletRate = ncol(hsmb)*8*1e-6 #��ÿ����1000��ϸ����˫ϸ����������ǧ��֮8������
DoubletRate = 0.2
#����ͬԴ˫ϸ������������modelHomotypic()�еĲ�����Ϊ���˫ϸ���������Ǵ�seurat_clusters������˫ϸ�� 
nExp_poi <- round(DoubletRate*length(hsmb$seurat_clusters))  #����ṩcelltype��������seurat_clusters��
# ����˫ϸ������
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## ���ʹ��ȷ���õĲ�������Doublets. Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
#hsmb <- doubletFinder_v3(hsmb, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
hsmb <- doubletFinder_v3(hsmb, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)
# ʹ��nExp = nExp_poi��nExp = nExp_poi.adj,�ֱ����doublets�������Ա����ȷ����Щϸ����Doublet-High Confidience
DimPlot(hsmb, reduction = "umap", group.by ="DF.classifications_0.25_0.005_2186")

saveRDS(hsmb,'hsmb_2.rds')