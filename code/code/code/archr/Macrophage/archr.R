
setwd('/public/home/xwtang/work/RCC/archr/trajectory/macrophage')
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(parallel)
set.seed(1)

addArchRThreads(threads = 10) 
addArchRGenome("hg38")

##multi cell
multi<-readRDS('/public/home/xwtang/work/RCC/archr/data/combined2/Save-ProjHeme2/Save-ArchR-Project.rds')
##macro
rcc.m<-readRDS('/public/home/xwtang/work/RCC/immune_cell/macrophage/data/rcc.m_2840.rds')
cellname<-colnames(rcc.m)
cellname<-gsub('_','#',cellname)

macro<-multi[cellname,]

#celltype
macro$celltype.refined<-rcc.m$celltype.refined

#cluster
macro <- addClusters(
  input = macro,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "Clusters_macro",
  resolution = 0.8
)

#umap
macro <- addUMAP(
  ArchRProj = macro, 
  reducedDims = "Harmony", 
  name = "UMAP_macro", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)
#magic
macro <- addImputeWeights(macro)

#umap
umap_m<-readRDS('/public/home/xwtang/work/RCC/archr/trajectory/macrophage/umap.rds')
rownames(umap_m)<-gsub('_','#',rownames(umap_m))
umap_m<-umap_m[cellname,]
macro@embeddings$UMAP_macro$df[rownames(umap_m),]<-umap_m
macro@embeddings$UMAP$df[rownames(umap_m),]<-umap_m

#####trajectory
trajectory <- c("Macro3", "Macro1", "Macro2")
trajectory

macro <- addTrajectory(
  ArchRProj = macro, 
  name = "macro_trajectory", 
  groupBy = "celltype.refined",
  trajectory = trajectory, 
  embedding = "UMAP_macro", 
  force = TRUE
)

#heatmap
trajMM  <- getTrajectory(ArchRProj = macro, name = "macro_trajectory", useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
trajGSM <- getTrajectory(ArchRProj = macro, name = "macro_trajectory", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))
trajGIM <- getTrajectory(ArchRProj = macro, name = "macro_trajectory", useMatrix = "GeneIntegrationMatrix", log2Norm = FALSE)
p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"))
trajPM  <- getTrajectory(ArchRProj = macro, name = "macro_trajectory", useMatrix = "PeakMatrix", log2Norm = TRUE)
p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))
plotPDF(p1 ,p2,p4,name = "Plot-LymphoidU-Traj-Heatmaps.pdf", ArchRProj = t, addDOC = FALSE, width = 6, height = 8)


# Integrative pseudo-time analyses
corGSM_MM <- correlateTrajectories(trajGSM, trajMM)
corGSM_MM[[1]]$matchname1
corGSM_MM[[1]]

trajGSM2 <- trajGSM[corGSM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGSM_MM[[1]]$name2, ]

trajCombined <- trajGSM2
assay(trajCombined,withDimnames = F) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))

combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)

rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))
ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2,  pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)
plotPDF(ht1 ,ht2,name = "Plot-LymphoidU-Traj-tcell-Heatmaps.pdf", ArchRProj = t, addDOC = FALSE, width = 6, height = 8)
saveArchRProject(ArchRProj = t, outputDirectory = "t", load = FALSE)
