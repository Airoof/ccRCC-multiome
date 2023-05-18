setwd('E:\\work\\RCC\\data\\normala')

set.seed(1)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)
library(BiocGenerics)

counts_t <- Read10X_h5("outs/filtered_feature_bc_matrix.h5")
# create a Seurat object containing the RNA adata
normala <- CreateSeuratObject(
  counts = counts_t$`Gene Expression`,
  assay = "RNA",project = "normala",min.cells=3,min.features=200
)
#使用MT-作为线粒体基因集
normala[["percent.mt"]]<-PercentageFeatureSet(normala,pattern = "^MT-")
#qc指标可视化
VlnPlot(normala,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
normala<-subset(normala,subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt<15 & nCount_RNA>1000 & nCount_RNA<50000)

#数据归一化   
normala<-NormalizeData(normala)
#高变异基因选择
normala<-FindVariableFeatures(normala,selection.method = "vst",nfeatures = 2000)

#缩放数据
all.genes<-rownames(normala)
normala<-ScaleData(normala,features = all.genes)
#降维
normala<-RunPCA(normala,features = VariableFeatures(object = normala))
ElbowPlot(normala)
#计算最邻近距离
normala<-FindNeighbors(normala,dims=1:20)
#聚类，包含设置下游聚类的“间隔尺度”的分辨率参数resolution，增加值会导致更多的聚类
normala<-FindClusters(normala,resolution = 0.5)
normala<-RunUMAP(normala,dims=1:20) #saveRDS(normala,"normala.rds")
DimPlot(normala,reduction = "umap",label = T)

##注释
#0  LOH loop of Henle  slc12a1 Clcnka 升肢细段
#1  LOH loop of Henle  slc12a1
#2  PCT proximal convoluted tubule CUBN/LRP2/SLC34A1/SLC5A12/ALDOB 
#3  LOH-PCT
#4  SLC44A5+VCAM1+ EPI proximal convoluted tubule
#5  CD-PC collecting duct principal cell  aqp2
#6  PST  SLC22A7	Proximal straight tubule
#7  SLC44A5+ EPI proximal convoluted tubule
#8  DCT distal convoluted tubule  slc12a3 
#9  CD-ICA collecting duct intercalated cell A  aqp6/atp6v0d2
#10 ENDO endothelial kdr/pecam1/PLVAP/flt1
#11 FIB ECMN
#12 MES mesangial cell pdgfrb
#13 LEU ptprc
#14 DC-ICB collecting duct intercalated cell B  slc26a4/atp6v0d2
#15 PODO podocyte  nphs1/nphs2

type<-c('0'='LOH','1'='LOH','2'='PCT','3'='LOH-PT','4'='SLC44A5+VCAM1+ EPI','5'='CD-PC','6'='PST',
        '7'='SLC44A5+ EPI','8'='DCT','9'=
          'DC-ICA','10'='ENDO','11'='FIB','12'='MES',
        '13'='LEU','14'='DC-ICB','15'='PODO')
normala[['celltype']] = unname(type[normala@meta.data$seurat_clusters])
markers<-FindAllMarkers(normala)
saveRDS(markers,'markers.rds')
saveRDS(normala,'normala.rds')

####
cg1=c('SLC34A1','LRP2','HAVCR1','CFH','SLC12A1','SLC12A3','SLC8A1','AQP2','SLC26A7','SLC26A4','NPHS2','EMCN','PIEZO2','COL1A2','PTPRC')
p1 <- DotPlot(normala, features = unique(cg1),
             assay='RNA'  )  + coord_flip()
p1

cg2=c('CUBN','LRP2','SLC34A1','SLC5A12','SLC5A2','ALDOB','CFH','SLC12A1','SLC12A3','SLC12A2','SLC8A1','AQP2','AQP6','SLC26A4','ATP6V0D2','NPHS1','NPHS2','PECAM1','FLT1','PDGFRB','PTPRC')
p2 <- DotPlot(normala, features = unique(cg2),
             assay='RNA'  )  + coord_flip()
p2 

#######
normala$celltype<-factor(normala$celltype,levels = c('CD-PC','DC-ICA','DC-ICB','DCT','DCT-CT','LOH','LOH-PT','PCT','PST','SLC44A5+ EPI',
                                                     'SLC44A5+VCAM1+ EPI','ENDO','FIB','MES','PODO','LEU'),ordered = T)
cg=c('AQP2','AQP6','ATP6V0D2','SLC26A4','SLC12A3','SLC12A2','SLC8A1','SLC12A1','CUBN','LRP2','SLC34A1','SLC5A12','ALDOB','SLC22A8','SLC22A7',
     'SLC44A5','VCAM1','PECAM1','FLT1','EMCN','PDGFRB','CFH','NPHS1','NPHS2','PTPRC')
DotPlot(normala, features = unique(cg),
        assay='RNA' ,group.by = 'celltype' )  + coord_flip()+ RotatedAxis()

################细胞注释###########
library(SingleR)
library(tibble)
hpca.se<-HumanPrimaryCellAtlasData()
bpe.se<-BlueprintEncodeData()
dic.se<-DatabaseImmuneCellExpressionData()
obj<-normala
obj@meta.data$cell.type<-Idents(obj)
test_1<-as.SingleCellExperiment(obj)
anno<-SingleR(test=test_1,ref = list(HP=hpca.se,BP=bpe.se),labels = list(hpca.se$label.fine,bpe.se$label.fine),method = "cluster",cluster=test_1$cell.type)
anno$cluster<-rownames(anno)
fin<-anno %>% dplyr::tbl_df() %>% dplyr::select(cluster,labels)
new.cluster.ids<-fin$labels
names(new.cluster.ids)<-levels(obj)
obj<-RenameIdents(obj,new.cluster.ids)
DimPlot(obj,reduction = "umap",label = TRUE,pt.size = 0.5)+NoLegend()


###cluster
####富集分析
marker6<-subset(markers,cluster==6)
colnames(marker6)[7]<-'SYMBOL'
#gsea
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

df_id<-bitr(marker6$SYMBOL, #转换的列是df数据框中的SYMBOL列
            fromType = "SYMBOL",#需要转换ID类型
            toType = "ENTREZID",#转换成的ID类型
            OrgDb = "org.Hs.eg.db")#对应的物种，小鼠的是org.Mm.eg.db

df_all<-merge(marker6,df_id,by="SYMBOL",all=F)#使用merge合并

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
