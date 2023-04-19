
setwd('E:\\work\\RCC\\cellchat\\alltype')
library(CellChat)

cellchat<-readRDS("E:\\work\\RCC\\cellchat\\alltype\\cellchat.rds")

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions",vertex.label.cex=0.5)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",vertex.label.cex=0.5)

########
mat <- cellchat@net$weight
par(mfrow =c(3,3),xpd=T)

for (i in 1:nrow(mat)){
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat),title.name = rownames(mat)[i],vertex.label.cex=0.5)
}

######
cellchat@netP$pathways
pathways.show<-c('VEGF')

# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object

#Compute the contribution of each ligand-receendoor pair to the overall signaling pathway and 
#visualize cell-cell communication mediated by a single ligand-receendoor pair
netAnalysis_contribution(cellchat, signaling = pathways.show)

#visualize the cell-cell communication mediated by a single ligand-receendoor pair
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[c(1,2),] # show one ligand-receendoor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver,sources.use =c(1,10,11,19,21),targets.use = c(2:9,12:18,20))
#> [[1]]
# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
# Chord diagram
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

#Visualize cell-cell communication mediated by multiple ligand-receendoors or signaling pathways
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = c(6), signaling = pathways.show,remove.isolate = FALSE)
# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = c(1,15,16,34,51), targets.use = c(48), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = c(21,22,23,24,27,2,36), targets.use = c(48), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = c(4,5,6,7,8,46), targets.use = c(48), remove.isolate = FALSE)
# show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","TGFb","EGF"))
netVisual_bubble(cellchat, sources.use = c(1), targets.use = c(5:8), pairLR.use = pairLR.use, remove.isolate = TRUE)

#Plot the signaling gene expression distribution using violin/dot plot
plotGeneExpression(cellchat, signaling = "CXCL")


###########Systems analysis of cell-cell communication network
#Compute and visualize the network centrality scores
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 16, height = 5, font.size = 10)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("MIF", "CCL"))

######Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together
library(NMF)
library(ggalluvial)
#Identify and visualize outcoming communication pattern of target cells
selectK(cellchat, pattern = "outgoing")
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "outgoing")
# dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")

#Identify and visualize incoming communication pattern of target cells
selectK(cellchat, pattern = "incoming")
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "incoming")
# dot plot
netAnalysis_dot(cellchat, pattern = "incoming")

##Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)

#Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)


###################subset
cellchat@meta$subtype[cellchat@meta$subtype =='Macro0' ] <- 'TMA-HIF1A' 
cellchat@meta$subtype[cellchat@meta$subtype =='Macro1' ] <- 'TMA-CCL3'
cellchat@meta$subtype[cellchat@meta$subtype =='Macro2' ] <- 'TMA-CD163' 
cellchat@meta$subtype[cellchat@meta$subtype =='Macro3' ] <- 'TMA-CD74' 
cellchat@meta$subtype[cellchat@meta$subtype =='Macrophage(proliferating)' ] <- 'TMA-MKI67' 
cellchat <- setIdent(cellchat, ident.use = "subtype") # set "cell_type" as default cell identity
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

##macrophage endo
subcells<-subsetCellChat(cellchat,idents.use = c('AEC','GCEC','PCEC','VEC','ENDO(proliferating)','TMA-HIF1A','TMA-CCL3','TMA-CD163','TMA-CD74','TMA-MKI67'),group.by = 'subtype')
netVisual_aggregate(subcells, signaling = pathways.show,  vertex.receiver = vertex.receiver)
netVisual_aggregate(subcells, signaling = pathways.show, layout = "chord")
netVisual_heatmap(subcells, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_contribution(subcells, signaling = pathways.show)
pairLR.CXCL <- extractEnrichedLR(subcells, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receendoor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(subcells, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
# Chord diagram
netVisual_individual(subcells, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

netVisual_bubble(subcells,sources.use = c(5:9),targets.use = c(1:4,10), remove.isolate = FALSE)

###all subtype
subcells<-subsetCellChat(cellchat,idents.use = c("AEC", "B cell", "CD16+ monocyte",           
                                                 "CD4+ Tn","CD8+ T(proliferating)","CD8+ Temra",               
                                                 "CD8+ Tex","CD8+ Trm","CNT",                      
                                                 "DC","DCT",             
                                                 "ENDO-PT","ENDO(proliferating)",      
                                                 "GCEC","ICA","ICB",                      
                                                 "LOH", "LOH-PT",                  
                                                 'TMA-HIF1A','TMA-CCL3','TMA-CD163','TMA-CD74','TMA-MKI67',                   
                                                 "Mast","MES","NKbright",                 
                                                 "NKdim","PC",                       
                                                 "PCEC","PEC","Plasma",                   
                                                 "PODO","PROM1_PT-S1","PROM1_PT-S2",              
                                                 "PROM1_PT-S3","PT-S1","PT-S1/S2",                 
                                                 "PT-S2","PT-S3",                    
                                                 "Treg","type 1","type 2",                   
                                                 "type 3","type 4","VEC" ),group.by = 'subtype')

pdf('cellchat.pdf',width = 300,height = 100)
netVisual_bubble(cellchat,remove.isolate = FALSE)
dev.off()

#macro_others
netVisual_bubble(subcells,sources.use = c(35),remove.isolate = FALSE)
netVisual_bubble(subcells,sources.use = c(36),remove.isolate = FALSE)
netVisual_bubble(subcells,sources.use = c(37),remove.isolate = FALSE)
netVisual_bubble(subcells,sources.use = c(38),remove.isolate = FALSE)
netVisual_bubble(subcells,sources.use = c(39),remove.isolate = FALSE)
#aec macro
netVisual_bubble(subcells,sources.use = c(35:39),targets.use = c(1),remove.isolate = FALSE)
netVisual_bubble(subcells,sources.use = c(1),targets.use = c(35:39),remove.isolate = FALSE)

####
df.net<-subsetCommunication(cellchat)
df<-data.frame(matrix(0,138,5))
colnames(df)<-c('cancertype','subtype','celltype','ligand_count','receptor_count')
df$cancertype<-c(rep('type 2',46),rep('type 3',46),rep('type 4',46))
subtype<-unique
df$subtype<-rep(c(unique(df.net$source))[!c(unique(df.net$source))%in%c('type 1','type 2','type 3','type 4','other cancer cell')],3)

for(i in 1:138){
  ligand<-subset(df.net,source==df[i,1]& target==df[i,2])
  receptor<-subset(df.net,source==df[i,2]&target==df[i,1])
  df[i,4]<-nrow(ligand)
  df[i,5]<-nrow(receptor)
}
write.csv(df,'lrcount.csv',quote = F,row.names = F)

###plot
df<-read.csv('lrcount.csv')
plot_data<-subset(df,cancertype=='Tumor_FOS')
rownames(plot_data)<-plot_data$subtype
colnames(plot_data)[4:5]<-c('Ligand','Receptor')

library(pheatmap)
library(RColorBrewer)
anno_row<-data.frame(plot_data[,3])
rownames(anno_row)<-rownames(plot_data)
pheatmap(plot_data[,4:5],cluster_rows = F,cluster_cols = F,display_numbers = T,number_format ="%.0f",gaps_row = c(5,15,23),
         #annotation_row = anno_row,annotation_legend = F,annotation_names_row = F
         border_color = NA,color=colorRampPalette(c("white", "#CD0000"))(200)
         )



########endo
subcells<-subsetCellChat(cellchat,idents.use = c("AEC", "B cell", "CD16+ monocyte",           
                                                 "CD4+ Tn","CD8+ T(proliferating)","CD8+ Temra",               
                                                 "CD8+ Tex","CD8+ Trm",            
                                                 "ENDO(proliferating)",      
                                                 "GCEC",                
                                                 'TMA-HIF1A','TMA-CCL3','TMA-CD163','TMA-CD74','TMA-MKI67',                   
                                                 "Mast","MES","NKbright",                 
                                                 "NKdim",                       
                                                 "PCEC","Plasma",                   
                                                                    
                                                 "Treg","VEC" ),group.by = 'subtype')
df.net<-subsetCommunication(subcells)
df<-data.frame(matrix(0,23,23))
colnames(df)<-unique(subcells@idents)
rownames(df)<-unique(subcells@idents)
for(i in colnames(df)){
  for(j in rownames(df)){
    subdf1<-subset(df.net,source==i)
    subdf2<-subset(df.net,source==j)
    df[j,i]<-sum(subdf1$target==j)+sum(subdf2$target==i)
  }
}

#heatmap
library(pheatmap)
pheatmap(log(df),border_color = NA,color = colorRampPalette(colors = c("#2A5A9A",'#EAD0BA',"#BF3939"))(200),
         legend = TRUE)

#net plot
targets<-c("B cell", "CD16+ monocyte",           
           "CD4+ Tn","CD8+ T(proliferating)","CD8+ Temra",               
           "CD8+ Tex","CD8+ Trm",
           'TMA-HIF1A','TMA-CCL3','TMA-CD163','TMA-CD74','TMA-MKI67',                   
           "Mast","MES","NKbright",                 
           "NKdim", "Plasma",
           "Treg")
sources<-c('AEC','ENDO(proliferating)','GCEC','PCEC','VEC')
df<-data.frame(source=c(rep('AEC',18),targets),target=c(targets,rep('AEC',18)),num=rep(0,36))
for(i in 1:36){
  subdf<-subset(df.net,source==df$source[i] & target==df$target[i])
  df$num[i]<-nrow(subdf)
}

library(tidyverse)
library(igraph)
library(dplyr)
name<-data.frame(c(df$source,df$destination))
nodes<-name%>%
  distinct()%>%
  mutate(location=c("western"))
colnames(nodes)<-c("label","location")
edges<-df%>%
  rename(from=source,to=target,weight=num)

net_pc<-graph_from_data_frame(
  d=edges,vertices=nodes,
  directed=TRUE)
V(net_pc)
plot(net_pc)

layouts <- grep("^layout_", ls("package:igraph"), value=TRUE)[-1]
l <- do.call(layouts[2], list(net_pc))
plot(net_pc,layout=l)

deg<-c()
for(i in 1:18){
  deg<-c(deg,df$num[i]+df$num[i+18])
}
deg<-c(150,deg)
vcolor<-c('#8CBB80FF',"#9AD3EBFF","#2F7CA6FF","#9AD3EBFF","#9AD3EBFF","#9AD3EBFF","#9AD3EBFF","#9AD3EBFF",
          "#2F7CA6FF","#2F7CA6FF","#2F7CA6FF","#2F7CA6FF","#2F7CA6FF","#125A84FF","#8CBB80FF","#9AD3EBFF",
          "#9AD3EBFF","#2F7CA6FF","#9AD3EBFF")
V(net_pc)$color<-vcolor[factor(V(net_pc))]
E(net_pc)$width<-deg*0.08
plot(net_pc,vertex.size=deg*0.2,
     vertex.label.cex=.7,vertex.label.dist=0,
     edge.color=c(rep('#0e88ff',18),rep('#f25f30',18)),edge.arrow.size=deg*0.005, edge.curved=.1,layout=l)


##########dotplot
subdf<-subset(df.net,source=='AEC' & (target %in% targets))
shoupei<-table(subdf$interaction_name_2)[order(table(subdf$interaction_name_2),decreasing = T)]
shoupei<-names(shoupei)[1:10]
subdf<-subset(subdf,interaction_name_2%in% shoupei)

df<-data.frame(cell=rep(unique(subdf$target),each=10),Ligand_Receptor=rep(shoupei,14),pval=rep(0,140),Prob=rep(0,140))
for(i in 1:140){
  tmp<-subset(subdf,target==df[i,1] & interaction_name_2==df[i,2])
  if(nrow(tmp)==0){
    next
  }
  df[i,3]<-tmp$pval
  df[i,4]<-tmp$prob
}

ggplot(data = df, aes(x = cell, y =Ligand_Receptor ))+scale_color_gradientn(values = seq(0,1,0.2),colours = c('#003366','#FBFCB5','brown3'))+
  theme_bw()+
  geom_point(aes(size=Pval,color=Prob),size = 10)+
  theme(panel.grid = element_blank(),axis.text.x =element_text(angle =90,hjust =0.5,vjust = 0.5))+
  xlab(NULL) + ylab(NULL)+
  theme(panel.border = element_rect(fill=NA,color = 'black',size = 2,linetype = 'solid'))+
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 90, hjust = 1, vjust = 0.5, face = "bold"),
        axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0.5, face = "bold"))





