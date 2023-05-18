library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)

set.seed(1234)

counts <- Read10X("./")
fragpath <- "../atac_fragments.tsv.gz"
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

# create a Seurat object containing the RNA adata
pbmc <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

DefaultAssay(pbmc) <- "ATAC"
# compute nucleosome signal score per cell
pbmc <- NucleosomeSignal(object = pbmc)
pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = pbmc, group.by = 'nucleosome_group', region = 'chr1-1-10000000')

# compute TSS enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc)

# add blacklist ratio and fraction of reads in peaks
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments


pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)
DepthCor(pbmc)

features <- VariableFeatures(pbmc)
data.use <- GetAssayData(
  object = pbmc[['ATAC']],
  slot = "data"
)[features, ]

pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
DimPlot(object = pbmc, label = TRUE) + NoLegend()
saveRDS(pbmc, 'rcca.rds')

ppb <- InsertionBias(pbmc[['ATAC']], BSgenome.Hsapiens.UCSC.hg38)

chr1 <- getSeq(BSgenome.Hsapiens.UCSC.hg38, 'chr1')

DefaultAssay(pbmc) <- "ATAC"
length(VariableFeatures(pbmc))
annotation <- Annotation(pbmc[['ATAC']])

pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(pbmc, group.by = 'high.tss') + NoLegend()

VlnPlot(
  object = pbmc,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg19)

# extract position frequency matrices for the motifs
pwm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

seqlengths(BSgenome.Hsapiens.UCSC.hg19)

print_jaspar <- function(x){
  xname <- x@name
  xname <- gsub('::', '_', xname)
  filename <- paste0(xname, '.jaspar')
  cat(paste0(">", x@ID, "\t", x@name, '\n'), file = filename)
  cat(paste0('A [ ',  paste0(x@profileMatrix[1, ], collapse = '\t'), ' ]\n'), file = filename, append = T)
  cat(paste0('C [ ',  paste0(x@profileMatrix[2, ], collapse = '\t'), ' ]\n'), file = filename, append = T)
  cat(paste0('G [ ',  paste0(x@profileMatrix[3, ], collapse = '\t'), ' ]\n'), file = filename, append = T)
  cat(paste0('T [ ',  paste0(x@profileMatrix[4, ], collapse = '\t'), ' ]\n'), file = filename, append = T)
}

nn <- sapply(namess, function(xname){
  xname <- strsplit(xname, split = '\\(')[[1]][1]
  xname <- gsub('::', '|', xname)
  xname
})

sapply(pwm, print_jaspar)

print_jaspar(pwm[[1]])

namess <- sapply(pwm, function(x){
  cat(A [])
})

pwm <- getMatrixSet(
  x = JASPAR2020,
  opts = list()
)

ppp <- as.data.table(pwm)

cancer <- readRDS('rcca.rds')

# add motif information
DefaultAssay(pbmc) <- 'ATAC'
cancer <- AddMotifs(pbmc, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pwm)

bg_freqs <- c(0.1, 0.4, 0.25, 0.25)
matrix(log2(bg(pwm[[1]])) - log2(bg_freqs), nrow = 4,
       ncol = length(pwm[[1]]),
       byrow = FALSE)
asa <- RegionStats(pbmc, genome = BSgenome.Hsapiens.UCSC.hg38)
pbmc@assays$ATAC@bias

genome_freq <- Biostrings::oligonucleotideFrequency(
  x = Biostrings::getSeq(x = BSgenome.Hsapiens.UCSC.hg38, "chr1"),
  width = 6
)
reads <- MultiGetReadsInRegion(
  object = pbmc,
  region = 'chr1-1-249250621',
)
