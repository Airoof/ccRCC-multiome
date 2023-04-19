
setwd('/public/home/xwtang/work/RCC/archr/data/combined2')
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(parallel)
set.seed(1)

addArchRThreads(threads = 30) 
addArchRGenome("hg38")

inputFiles<-c('/public/home/xwtang/work/RCC/archr/fragment/rcca/atac_fragments.tsv.gz','/public/home/xwtang/work/RCC/archr/fragment/rccd/atac_fragments.tsv.gz',
              '/public/home/xwtang/work/RCC/archr/fragment/hsma/atac_fragments.tsv.gz','/public/home/xwtang/work/RCC/archr/fragment/hsmb/atac_fragments.tsv.gz',
              '/public/home/xwtang/work/RCC/archr/fragment/normala/atac_fragments.tsv.gz','/public/home/xwtang/work/RCC/archr/fragment/normald/atac_fragments.tsv.gz')


ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = c('RCCa','RCCd','HSMa','HSMb','Normala','Normald'),
  minTSS = 1, #Dont set this too high because you can always increase later
  minFrags = 1000,
  maxFrags =  50000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

projHeme1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "/public/home/xwtang/work/RCC/archr/data/combined2/comeout",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

getAvailableMatrices(projHeme1)

#meta
meta<-readRDS('/public/home/xwtang/work/RCC/data/intergrated/meta.rds')
rownames(meta)<-gsub('_','#',rownames(meta))
intecell<-intersect(rownames(meta),rownames(projHeme1))

meta<-meta[intecell,]
projHeme1<-projHeme1[intecell,]
projHeme1$celltype<-meta$celltype
projHeme1$celltype.refined<-meta$celltype.refined
projHeme1$tissue<-meta$tissue

#plot
df <- getCellColData(projHeme1, select = c("log10(nFrags)", "TSSEnrichment"))
pdf('/public/home/xwtang/work/RCC/archr/pic/qcplot.pdf')
ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
dev.off()

pdf('/public/home/xwtang/work/RCC/archr/pic/ridgeplot.pdf')
plotGroups(
  ArchRProj = projHeme1, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "ridges"
)
dev.off()

projHeme1 <- addIterativeLSI(
  ArchRProj = projHeme1,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30
)

#Subtle batch effect,add interations
projHeme1 <- addIterativeLSI(
  ArchRProj = projHeme1,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI2", 
  iterations = 4, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.1, 0.2, 0.4), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 15000, 
  dimsToUse = 1:30
)

#strong batch effect, use harmony
projHeme1 <- addHarmony(
  ArchRProj = projHeme1,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample"
)

#cluster
projHeme1 <- addClusters(
  input = projHeme1,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8
)

#umap
projHeme1 <- addUMAP(
  ArchRProj = projHeme1, 
  reducedDims = "Harmony", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

p1 <- plotEmbedding(ArchRProj = projHeme1, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = projHeme1, colorBy = "cellColData", name = "Clusters", embedding = "UMAP") 

#plot
pdf('/public/home/xwtang/work/RCC/archr/pic/cluster.pdf',width = 10,height = 6)
ggAlignPlots(p1, p2, type = "h")
dev.off()

p1 <- plotEmbedding(ArchRProj = projHeme1, colorBy = "cellColData", name = "celltype", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = projHeme1, colorBy = "cellColData", name = "celltype.refined", embedding = "UMAP")
#plot
pdf('/public/home/xwtang/work/RCC/archr/pic/celltype.pdf',width = 10,height = 6)
ggAlignPlots(p1, p2, type = "h")
dev.off()

saveArchRProject(ArchRProj = projHeme1, outputDirectory = "/public/home/xwtang/work/RCC/archr/data/combined2/Save-ProjHeme1", load = FALSE)


#readrds
projHeme2<-readRDS('/public/home/xwtang/work/RCC/archr/data/combined2/Save-ProjHeme1/Save-ArchR-Project.rds')

#feature plot
pdf('/public/home/xwtang/work/RCC/archr/pic/cd247.pdf',width = 6,height =6 )
plotEmbedding(
  ArchRProj = projHeme2, 
  colorBy = "GeneScoreMatrix", 
  name = c('CD247','PTPRC'), 
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)
dev.off()

#magic
projHeme2 <- addImputeWeights(projHeme2)
#plot
pdf('/public/home/xwtang/work/RCC/archr/pic/cd247.pdf',width = 6,height =6 )
plotEmbedding(
  ArchRProj = projHeme2, 
  colorBy = "GeneScoreMatrix", 
  name = c('CD247','PTPRC'), 
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = getImputeWeights(projHeme2)
)
dev.off()

#track plot
markerGenes  <- c(
  "CD247", #Early Progenitor
  "PTPRC", #Erythroid
  "PECAM1", "PDGFRB", #B-Cell Trajectory
  "CFH", #Monocytes
  "CA9", "VCAM1", "CUBN", "LRP2" #TCells
)
p <- plotBrowserTrack(
  ArchRProj = tumor, 
  groupBy = "celltype.refined", 
  geneSymbol = 'AMH', 
  upstream = 1000,
  downstream = 5000
)

plotPDF(plotList = p, 
        name = "Plot-Tracks-AMH-20220812-tumor.pdf", 
        ArchRProj = projHeme2, 
        addDOC = FALSE, width = 5, height = 5)

#pseudo bulk
addArchRThreads(threads = 1) 
projHeme2 <- addGroupCoverages(ArchRProj = projHeme2, groupBy = "celltype.refined")

###macs2
pathToMacs2<-findMacs2()
projHeme2 <- addReproduciblePeakSet(
  ArchRProj = projHeme2, 
  groupBy = "celltype.refined", 
  pathToMacs2 = pathToMacs2
)
addArchRThreads(threads = 10) 
projHeme2<-addPeakMatrix(projHeme2)

saveArchRProject(ArchRProj = projHeme2, outputDirectory = "/public/home/xwtang/work/RCC/archr/data/combined2/Save-ProjHeme2", load = FALSE)
