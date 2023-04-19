
setwd('E:\\work\\RCC')

library(Seurat)
normal.combined<-readRDS('data/combined/normal_combined_4.rds')
pt<-subset(normal.combined,celltype %in% c('PT','PROM1_PT'))
DefaultAssay(pt)<-'integrated'

#plot
DimPlot(pt,label = T)

pt<-FindVariableFeatures(pt,nFeatures=2000)
pt <- ScaleData(pt, verbose = FALSE)
pt <- RunPCA(pt, npcs = 30, verbose = FALSE)
pt <- RunUMAP(pt, reduction = "pca", dims = 1:15)
pt <- FindNeighbors(pt, reduction = "pca", dims = 1:15) 
pt <- FindClusters(pt, resolution = 0.5)
DimPlot(pt,label = T)

##
DefaultAssay(pt)<-'RNA'
markers<-FindAllMarkers(pt,only.pos = T)
markers<-subset(markers,p_val_adj<0.05)

################
#0 	man1c1/slc36a2/slc22a8   pct-S1
#1  slc5a8/slc5a1/SLC7A13/slc22a7 pst
#2  man1c1/slc36a2/slc22a8   pct-S1
#3  slc22a7/slc22a8  pct-pst
#4  slc5a8/slc5a1/SLC7A13/slc22a7 pst
#5  pt-progenitor(resting)
#6  pt-progenitor(activated)
#7  habp2/runx1/cdh11  pct-s2
#8  pt-progenitor(pre-activated)
#9  

pt$celltype.refined<-pt$celltype
pt$celltype.refined<-ifelse(pt$seurat_clusters%in%c(0,2),'PT-S1',pt$celltype.refined)
pt$celltype.refined<-ifelse(pt$seurat_clusters%in%c(1,4),'PT-S2',pt$celltype.refined)
pt$celltype.refined<-ifelse(pt$seurat_clusters==3,'PT-S1/S2',pt$celltype.refined)
pt$celltype.refined<-ifelse(pt$seurat_clusters==5,'PROM1_PT-S1',pt$celltype.refined)
pt$celltype.refined<-ifelse(pt$seurat_clusters==6,'PROM1_PT-S3',pt$celltype.refined)
pt$celltype.refined<-ifelse(pt$seurat_clusters==7,'PT-S3',pt$celltype.refined)
pt$celltype.refined<-ifelse(pt$seurat_clusters==8,'PROM1_PT-S2',pt$celltype.refined)

#select cell
pt9=subset(pt,seurat_clusters==9)
p=DimPlot(pt9)
cell=CellSelector(p)
pt$celltype.refined<-ifelse(colnames(pt)%in% cell,'PT-S1/S2',pt$celltype.refined)
cell=CellSelector(p)
pt$celltype.refined<-ifelse(colnames(pt)%in% cell,'PT-S2',pt$celltype.refined)
cell=CellSelector(p)
pt$celltype.refined<-ifelse(colnames(pt)%in% cell,'PROM1_PT-S1',pt$celltype.refined)
cell=CellSelector(p)
pt$celltype.refined<-ifelse(colnames(pt)%in% cell,'PT-S1',pt$celltype.refined)

saveRDS(pt,'immune/PT/pt.rds')
###enrichment
library(org.Hs.eg.db)
library(clusterProfiler)
library(patchwork)
library(tidyverse)
marker<-FindMarkers(pt,ident.1 = 6)
marker<-subset(marker,p_val_adj<0.05)
marker$SYMBOL<-rownames(marker)
ego_ALL<-enrichGO(gene = row.names(marker),
                  OrgDb = 'org.Hs.eg.db',
                  keyType = 'SYMBOL',
                  ont = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.01,
                  qvalueCutoff = 0.05)
#kegg分析
options(clusterProfiler.download.method = "wininet")
genelist<-bitr(rownames(marker8),
               fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = 'org.Hs.eg.db')
genelist <- pull(genelist,ENTREZID)               
ekegg <- enrichKEGG(gene = genelist, organism = 'hsa')

##gsea
df_id<-bitr(marker$SYMBOL, #转换的列是df数据框中的SYMBOL列
            fromType = "SYMBOL",#需要转换ID类型
            toType = "ENTREZID",#转换成的ID类型
            OrgDb = "org.Hs.eg.db")#对应的物种，小鼠的是org.Mm.eg.db
df_all<-merge(marker,df_id,by="SYMBOL",all=F)#使用merge合并

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
ridgeplot(GO,10)

##############
pt<-readRDS('immune/PT/pt.rds')
cellname<-readRDS('data/cellname.rds')
pt<-pt[,cellname]

saveRDS(pt,'immune/PT/pt_9899.rds')

#########
cg<-c('SLC34A1','SLC22A8','CYP3A5','PCDH15','SORCS1','SLC5A1','SLC5A11','SLC22A7','SLC7A13','DGKB','FMO5',
      'RUNX1','AFF2','CACNA1E','ST18','TGFBR2','SLC12A1','SLIT3','NRG1','ROBO2','DCC','AFF2','VCAM1','SLC9A9','NRG3','NRG1')
pdf('test.pdf',width = 20,height = 20)
VlnPlot(pt,unique(cg))
dev.off()

pdf('test.pdf',width = 20,height = 20)
FeaturePlot(pt,unique(cg))
dev.off()

#######enrichment

markers<-FindAllMarkers(pt,only.pos = T)
setwd('E:\\work\\RCC\\marker\\pt')
getdiff<-function(i){
  t<-subset(markers,cluster==i)
  t<-subset(t,p_val_adj<0.05)
  colnames(t)[7]<-'SYMBOL'
  return(t)
}
for(i in unique(markers$cluster)[4:7]){
  diff<-getdiff(i)
  filter_deg<-subset(diff,avg_log2FC>0.3)
  
  ####富集分析
  #gsea
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  library(patchwork)
  library(tidyverse)
  ego_ALL<-enrichGO(gene = row.names(filter_deg),
                    OrgDb = 'org.Hs.eg.db',
                    keyType = 'SYMBOL',
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.05)
  #kegg分析
  options(clusterProfiler.download.method = "wininet")
  genelist<-bitr(rownames(filter_deg),
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = 'org.Hs.eg.db')
  genelist <- pull(genelist,ENTREZID)               
  ekegg <- enrichKEGG(gene = genelist, organism = 'hsa')
  
  if(nrow(ego_ALL@result)>0){write.csv(ego_ALL@result,paste0('E:\\work\\RCC\\marker\\pt\\enrich_GO_',as.character(i),'.csv'))}
  if(nrow(ekegg@result)>0){write.csv(ekegg@result,paste0('E:\\work\\RCC\\marker\\pt\\enrich_KEGG_',as.character(i),'.csv'))}
}

###########pt atac umap
ptatac<-readRDS('immune/PT/ptCell_multi.rds')
umap<-ptatac@reductions[["umap"]]@cell.embeddings
saveRDS(umap,'immune/PT/umap.rds')

