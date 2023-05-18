
setwd('E:\\work\\RCC\\data\\rccd')

#######
rccd<-readRDS('rccd_1.rds')
library(DoubletFinder)
#����һ��������Ѳ����Ĺ��̣������ٶ���
sweep.res.list_rccd <- paramSweep_v3(rccd, PCs = 1:30, sct = FALSE)
#head(sweep.res.list_rccd)
sweep.stats_rccd <- summarizeSweep(sweep.res.list_rccd, GT = FALSE)
bcmvn_rccd <- find.pK(sweep.stats_rccd) #���Կ�����Ѳ����ĵ�
## ������ѵĲ����ǣ�
mpK<-as.numeric(as.vector(bcmvn_rccd$pK[which.max(bcmvn_rccd$BCmetric)]))

##Homotypic Doublet Proportion Estimate -------------------------------------
annotations <- rccd@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)  
#DoubletRate = ncol(rccd)*8*1e-6 #��ÿ����1000��ϸ����˫ϸ����������ǧ��֮8������
DoubletRate = 0.2
#����ͬԴ˫ϸ������������modelHomotypic()�еĲ�����Ϊ���˫ϸ���������Ǵ�seurat_clusters������˫ϸ�� 
nExp_poi <- round(DoubletRate*length(rccd$seurat_clusters))  #����ṩcelltype��������seurat_clusters��
# ����˫ϸ������
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## ���ʹ��ȷ���õĲ�������Doublets. Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
#rccd <- doubletFinder_v3(rccd, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
rccd <- doubletFinder_v3(rccd, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)
# ʹ��nExp = nExp_poi��nExp = nExp_poi.adj,�ֱ����doublets�������Ա����ȷ����Щϸ����Doublet-High Confidience
DimPlot(rccd, reduction = "umap", group.by ="DF.classifications_0.25_0.005_2958")

saveRDS(rccd,'rccd_2.rds')