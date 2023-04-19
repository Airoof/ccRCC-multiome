
setwd('/public/home/xwtang/work/RCC/archr/trajectory/tcell')
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(parallel)
set.seed(1)

addArchRThreads(threads = 10) 
addArchRGenome("hg38")

#color
tmp<-c('#006EB0','#2DACD0','#8AD5C7','#D1EECC','#E1F4DD')
names(tmp)<-c(1,2,3,4,5)
ArchRPalettes[['lzj']]<-tmp
tmp_pal <- tmp
tmp_palOut <- colorRampPalette(tmp_pal)(250)

##multi cell
multi<-readRDS('/public/home/xwtang/work/RCC/archr/data/combined2/Save-ProjHeme2/Save-ArchR-Project.rds')
#umap cell
umap_t<-readRDS('/public/home/xwtang/work/RCC/archr/trajectory/tcell/umap.rds')
rownames(umap_t)<-gsub('_','#',rownames(umap_t))

t<-subsetArchRProject(
  ArchRProj = multi,
  cells = rownames(umap_t),
  outputDirectory ='/public/home/xwtang/work/RCC/archr/trajectory/tcell/archr' ,
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = TRUE
)

#subset, lymphocyte
rcc.t<-readRDS('/public/home/xwtang/work/RCC/immune_cell/T/data/rcc.t.rds')
tname<-gsub('_','#',colnames(rcc.t))
tcell<-intersect(rownames(umap_t),tname)
tcell<-gsub('#','_',tcell)
rcc.t<-rcc.t[,tcell]
t$celltype.refined<-rcc.t$celltype.refined

#cluster
t <- addClusters(
  input = t,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "Clusters_t",
  resolution = 0.8
)

#umap
t <- addUMAP(
  ArchRProj = t, 
  reducedDims = "Harmony", 
  name = "UMAP_t", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)
#magic
t <- addImputeWeights(t)

umap_t<-readRDS('/public/home/xwtang/work/RCC/archr/trajectory/tcell/umap.rds')
rownames(umap_t)<-gsub('_','#',rownames(umap_t))
t@embeddings$UMAP_t$df[rownames(umap_t),]<-umap_t
t@embeddings$UMAP$df[rownames(umap_t),]<-umap_t

#####trajectory
trajectory <- c("CD8+ Temra", "CD8+ Trm", "CD8+ Tex")
trajectory

t <- addTrajectory(
  ArchRProj = t, 
  name = "cd8t_trajectory", 
  groupBy = "celltype.refined",
  trajectory = trajectory, 
  embedding = "UMAP_t", 
  force = TRUE
)

#heatmap
trajMM  <- getTrajectory(ArchRProj = t, name = "cd8t_trajectory", useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = tmp_palOut)
trajGSM <- getTrajectory(ArchRProj = t, name = "cd8t_trajectory", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- plotTrajectoryHeatmap(trajGSM,  pal = tmp_palOut)
trajGIM <- getTrajectory(ArchRProj = t, name = "cd8t_trajectory", useMatrix = "GeneIntegrationMatrix", log2Norm = FALSE)
p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"))
trajPM  <- getTrajectory(ArchRProj = t, name = "cd8t_trajectory", useMatrix = "PeakMatrix", log2Norm = TRUE)
p4 <- plotTrajectoryHeatmap(trajPM, pal = tmp_palOut)
plotPDF(p1 ,p2,p4,name = "Plot-LymphoidU-Traj-Heatmaps_20220906.pdf", ArchRProj = t, addDOC = FALSE, width = 6, height = 8)


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
plotPDF(ht1 ,ht2,name = "Plot-LymphoidU-Traj-tcell-Heatmaps_20220713.pdf", ArchRProj = t, addDOC = FALSE, width = 6, height = 8)
saveArchRProject(ArchRProj = t, outputDirectory = "t", load = FALSE)

