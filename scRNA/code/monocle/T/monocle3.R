
setwd('E:\\work\\RCC')
library(Hmisc)
library(monocle3)
library(Seurat)

rcc.t<-readRDS('monocle/T/rcc.t_1666.rds')
rcc.t$celltype.refined<-ifelse(rcc.t$celltype.refined=='CD8+ Trm','CD8+ Tmem',rcc.t$celltype.refined)
rcc.t$celltype.refined<-ifelse(rcc.t$celltype.refined=='CD4+ Tn','CD4+ Tconv',rcc.t$celltype.refined)
rcc.t$celltype.refined<-ifelse(rcc.t$celltype.refined=='CD8+ Temra','CD8+ Teff',rcc.t$celltype.refined)
rcc.t$celltype.refined<-ifelse(rcc.t$celltype.refined=='Treg','CD4+ Treg',rcc.t$celltype.refined)

monot<-subset(rcc.t,celltype.refined %in% c('CD8+ Teff','CD8+ Tmem','CD8+ Tex'))
monot$celltype.refined<-factor(monot$celltype.refined,levels = c('CD8+ Teff','CD8+ Tmem','CD8+ Tex'))

data<-monot@assays$RNA@counts
cell_metadata<-monot@meta.data
gene_annotation<-data.frame(gene_short_name=rownames(data))
rownames(gene_annotation)<-rownames(data)
cds<-new_cell_data_set(data,
                       cell_metadata = cell_metadata,
                       gene_metadata = gene_annotation)
cds<-preprocess_cds(cds,num_dim = 50)
plot_pc_variance_explained(cds)
cds<-reduce_dimension(cds)
cds <- cluster_cells(cds) 
plot_cells(cds)

#seuratµ¼Èëumap
cds.embed<-cds@int_colData$reducedDims$UMAP
embed<-Embeddings(monot,reduction = 'umap')
embed<-embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP<-embed

#constrcuct trajectory
cds<-learn_graph(cds)
cds<-order_cells(cds)

#plot
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=T,
           label_branch_points=T,
           graph_label_size=3,
           group_label_size=4)
pseudotime <- pseudotime(cds, reduction_method = 'UMAP') 

saveRDS(cds,'E:/work/RCC/monocle/T/cds.rds')
