

runUMAP <- function(object, assay){
  DefaultAssay(object) <- assay
  object <- NormalizeData(object)
  if(assay == 'RNA'){
    normalized <- expm1(object@assays$RNA@data)
  }
  else{
    normalized <- expm1(object@assays$ATAC@data)
  }
  hvf.info <- FindVariableFeatures(normalized)
  
  hvf.info <- hvf.info[which(x = hvf.info[, 1, drop = TRUE] != 0), ]
  hvf.info <- hvf.info[order(hvf.info$vst.variance.standardized, decreasing = TRUE), , drop = FALSE]
  
  nfeatures <- as.integer(nrow(hvf.info) * 0.2)
  VariableFeatures(object) <- head(x = rownames(x = hvf.info), n = nfeatures)
  scaled <- normalized[VariableFeatures(object), ]
  rm <- rowMeans(scaled)
  scaled <- apply(scaled, 2, function(x){x - rm })
  rv <- apply(scaled, 1, var)
  rstd <- sqrt(rv)
  scaled <- apply(scaled, 2, function(x){x / rstd })
  
  pca <- RunPCA(scaled)
  object@reductions[[paste0(assay, '_pca')]] <- pca
  object <- RunUMAP(object, dims = 1:20, reduction = paste0(assay, '_pca'), reduction.name = paste0(assay, '_umap'), reduction.key = paste0(assay, 'UMAP_') )
  return(object)
}
