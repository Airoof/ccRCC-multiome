
setwd('E:\\work\\RCC\\data\\RCCa')

set.seed(1)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)
library(BiocGenerics)

counts_t <- Read10X_h5("outs/filtered_feature_bc_matrix.h5")
# create a Seurat object containing the RNA adata
rcca <- CreateSeuratObject(
  counts = counts_t$`Gene Expression`,
  assay = "RNA",project = "RCCa",min.cells=3,min.features=200
)
#ʹ��MT-��Ϊ���������
rcca[["percent.mt"]]<-PercentageFeatureSet(rcca,pattern = "^MT-")
#qcָ����ӻ�
VlnPlot(rcca,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
rcca<-subset(rcca,subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt<15 & nCount_RNA>1000 & nCount_RNA<50000)

# #���ݹ�һ��   
# rcca<-NormalizeData(rcca)
# #�߱������ѡ��
# rcca<-FindVariableFeatures(rcca,selection.method = "vst",nfeatures = 2000)
# 
# #��������
# all.genes<-rownames(rcca)
# rcca<-ScaleData(rcca,features = all.genes)

#��ά
rcca<-RunPCA(rcca,features = VariableFeatures(object = rcca))
ElbowPlot(rcca)
#�������ڽ�����
rcca<-FindNeighbors(rcca,dims=1:20)
#���࣬�����������ξ���ġ�����߶ȡ��ķֱ��ʲ���resolution������ֵ�ᵼ�¸���ľ���
rcca<-FindClusters(rcca,resolution = 0.5)
rcca<-RunUMAP(rcca,dims=1:20) #saveRDS(rcca,"rcca.rds")
DimPlot(rcca,reduction = "umap",label = T)

##ע��
#0  cancer cell ca9
#1  cancer cell ca9
#2  cancer cell ca9
#3  macrophage cd68/cd163/hla-dra
#4  endothelial kdr/pecam1/PLVAP
#5  cancer cell ca9
#6  T cell
#7  cancer cell ca9
#8  cancer cell ca9
#9  cancer cell ca9
#10 cancer cell ca9
#11 epithelial
#12 NK gnly/gzmb

type<-c('0'='cancer cell','1'='cancer cell','2'='cancer cell','3'='macrophage','4'='endothelial','5'='cancer cell','6'='T cell','7'='cancer cell','8'='cancer cell','9'=
        'cancer cell','10'='cancer cell','11'='epithelial','12'='NK')
rcca[['celltype']] = unname(type[rcca@meta.data$seurat_clusters])
markers<-FindAllMarkers(rcca)
saveRDS(markers,'markers.rds')
saveRDS(rcca,'rcca_3.rds')
################ϸ��ע��###########
library(SingleR)
library(tibble)
hpca.se<-HumanPrimaryCellAtlasData()
bpe.se<-BlueprintEncodeData()
dic.se<-DatabaseImmuneCellExpressionData()
obj<-rcca
obj@meta.data$cell.type<-Idents(obj)
test_1<-as.SingleCellExperiment(obj)
anno<-SingleR(test=test_1,ref = list(HP=hpca.se,BP=bpe.se),labels = list(hpca.se$label.main,bpe.se$label.main),method = "cluster",cluster=test_1$cell.type)
anno$cluster<-rownames(anno)
fin<-anno %>% dplyr::tbl_df() %>% dplyr::select(cluster,labels)
new.cluster.ids<-fin$labels
names(new.cluster.ids)<-levels(obj)
obj<-RenameIdents(obj,new.cluster.ids)
DimPlot(obj,reduction = "umap",label = TRUE,pt.size = 0.5)+NoLegend()


###
markers_rcca<-readRDS('markers.rds')
getdiff<-function(i){
  t<-subset(markers_rcca,cluster==i)
  t<-subset(t,p_val_adj<0.05)
  colnames(t)[7]<-'SYMBOL'
  return(t)
}
for(i in c(1,2,3,8,9)){
  diff<-getdiff(i)
  
  
  ####��������
  #gsea
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  df_id<-bitr(diff$SYMBOL, #ת��������df���ݿ��е�SYMBOL��
              fromType = "SYMBOL",#��Ҫת��ID����
              toType = "ENTREZID",#ת���ɵ�ID����
              OrgDb = "org.Hs.eg.db")#��Ӧ�����֣�С�����org.Mm.eg.db
  
  df_all<-merge(diff,df_id,by="SYMBOL",all=F)#ʹ��merge�ϲ�
  
  df_all_sort <- df_all[order(df_all$avg_log2FC, decreasing = T),]#�Ȱ���logFC��������
  gene_fc = df_all_sort$avg_log2FC #��foldchange���մӴ�С��ȡ����
  names(gene_fc) <- df_all_sort$ENTREZID #��������ȡ��foldchange���϶�Ӧ��ENTREZID
  
  KEGG <- gseKEGG(gene_fc, organism = "hsa") #�������������
  
  GO <- gseGO(
    gene_fc, #gene_fc
    ont = "BP",# "BP"��"MF"��"CC"��"ALL"
    OrgDb = org.Hs.eg.db,#����ע�ͻ���
    keyType = "ENTREZID",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",#pֵУ������
  )
  if(nrow(GO@result)>0){write.csv(GO@result,paste0('E:\\work\\RCC\\marker\\rcca\\GO_',as.character(i),'.csv'))}
  if(nrow(KEGG@result)>0){write.csv(KEGG@result,paste0('E:\\work\\RCC\\marker\\rcca\\KEGG_',as.character(i),'.csv'))}
}

########
#sct
# run sctransform
rcca<-readRDS('rcca_2.rds')
rcca <- SCTransform(rcca, vars.to.regress = "percent.mt", verbose = FALSE)