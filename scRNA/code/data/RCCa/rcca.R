
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
#使用MT-作为线粒体基因集
rcca[["percent.mt"]]<-PercentageFeatureSet(rcca,pattern = "^MT-")
#qc指标可视化
VlnPlot(rcca,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
rcca<-subset(rcca,subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt<15 & nCount_RNA>1000 & nCount_RNA<50000)

# #数据归一化   
# rcca<-NormalizeData(rcca)
# #高变异基因选择
# rcca<-FindVariableFeatures(rcca,selection.method = "vst",nfeatures = 2000)
# 
# #缩放数据
# all.genes<-rownames(rcca)
# rcca<-ScaleData(rcca,features = all.genes)

#降维
rcca<-RunPCA(rcca,features = VariableFeatures(object = rcca))
ElbowPlot(rcca)
#计算最邻近距离
rcca<-FindNeighbors(rcca,dims=1:20)
#聚类，包含设置下游聚类的“间隔尺度”的分辨率参数resolution，增加值会导致更多的聚类
rcca<-FindClusters(rcca,resolution = 0.5)
rcca<-RunUMAP(rcca,dims=1:20) #saveRDS(rcca,"rcca.rds")
DimPlot(rcca,reduction = "umap",label = T)

##注释
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
################细胞注释###########
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
  
  
  ####富集分析
  #gsea
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  df_id<-bitr(diff$SYMBOL, #转换的列是df数据框中的SYMBOL列
              fromType = "SYMBOL",#需要转换ID类型
              toType = "ENTREZID",#转换成的ID类型
              OrgDb = "org.Hs.eg.db")#对应的物种，小鼠的是org.Mm.eg.db
  
  df_all<-merge(diff,df_id,by="SYMBOL",all=F)#使用merge合并
  
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
  if(nrow(GO@result)>0){write.csv(GO@result,paste0('E:\\work\\RCC\\marker\\rcca\\GO_',as.character(i),'.csv'))}
  if(nrow(KEGG@result)>0){write.csv(KEGG@result,paste0('E:\\work\\RCC\\marker\\rcca\\KEGG_',as.character(i),'.csv'))}
}

########
#sct
# run sctransform
rcca<-readRDS('rcca_2.rds')
rcca <- SCTransform(rcca, vars.to.regress = "percent.mt", verbose = FALSE)
