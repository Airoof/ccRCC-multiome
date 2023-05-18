
setwd('E:\\work\\RCC')
library(Seurat)
library(monocle3)

pt<-readRDS('cells/PT/pt_9899.rds')
data<-pt@assays$RNA@counts
cell_metadata<-pt@meta.data
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
embed<-Embeddings(pt,reduction = 'umap')
embed<-embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP<-embed

#constrcuct trajectory
cds<-learn_graph(cds,use_partition = F)
cds<-order_cells(cds)

#plot
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=F,
           label_branch_points=F,
           graph_label_size=3,
           group_label_size=4,
           show_trajectory_graph = F)
pseudotime <- pseudotime(cds, reduction_method = 'UMAP') 

saveRDS(cds,'monocle/pt/cds.rds')
