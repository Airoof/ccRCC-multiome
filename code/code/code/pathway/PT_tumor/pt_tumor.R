
setwd('E:\\work\\RCC')

library(Seurat)
library(ggplot2)
library(GSVA)
library(readxl)
genedata<- read_excel("E:/work/111_resources/13287_2020_1973_MOESM1_ESM/Supplementary Table 1.xlsx",col_names = FALSE)
rownames(genedata)<-genedata$...1
#genedata<-genedata[,-1]
genedata<-t(genedata)
genedata<-genedata[-1,]
genedata<-data.frame(genedata)
gs<-list()
for(i in 1:26){
  t<-na.omit(genedata[,i])
  gs[[colnames(genedata)[i]]]<-t[1:length(t)]
}

library(clusterProfiler)
genelist<-read.gmt('pathway/pt_tumor/CSC-related pathways-1.gmt')
genelist<-genelist[!(genelist$gene==''),]
gcSample = split(genelist$gene,
                 genelist$term)

#######subset cell
rcc.combined<-readRDS('data/combined/rcc_combined_4_73477.rds')
meta<-readRDS('data/combined/meta.rds')
meta<-meta[colnames(rcc.combined),]
rcc.combined$subtype<-meta$subtype_cancer

pttumor<-subset(rcc.combined,subtype %in% c('PROM1_PT-S1','PROM1_PT-S3','PROM1_PT-S2','PT-S1','PT-S2','PT-S3',
                                            'PT-S1/S2','type 1','type 2','type 3','type 4'))

######aucell
data<-pttumor@assays$RNA@data
library(AUCell)
cells_rankings <- AUCell_buildRankings(data, nCores = 1)
cells_AUC <- AUCell_calcAUC(gcSample, cells_rankings,nCores = 1)
saveRDS(cells_AUC,'pathway/pt_tumor/cells_AUC_ycy.rds')

x<- cells_AUC@assays@data@listData[["AUC"]]
x<-t(x)
meta<-pttumor@meta.data
meta<-cbind(meta,x)
pttumor@meta.data<-meta
DimPlot(pttumor,group.by = 'subtype')


######retu
x<-t(x)
library(pheatmap)
anno_col<-meta[,1:2]
colnames(anno_col)[2]<-'celltype'
anno_col$celltype<-pttumor$subtype
anno_col$celltype<-factor(anno_col$celltype,levels = c('PROM1_PT-S1','PROM1_PT-S2','PROM1_PT-S3','PT-S1','PT-S1/S2','PT-S2','PT-S3','type 1','type 2','type 3','type 4'),ordered = T)
anno_col<-anno_col[order(anno_col$celltype),]
x<-x[,rownames(anno_col)]
y<-t(scale(t(x)))
y<-ifelse(y>3,3,y)
y<-ifelse(y<(-3),-3,y)
pheatmap(y,cluster_rows = F,cluster_cols = F,labels_col = FALSE,annotation_col = anno_col)

######ridge plot
library(ggridges)
library(sciplot)
meta$subtype<-factor(meta$subtype,levels = c('PROM1_PT-S1','PROM1_PT-S2','PROM1_PT-S3','PT-S1','PT-S1/S2','PT-S2','PT-S3','type 1','type 2','type 3','type 4'),ordered = T)
ggplot(meta,aes(x = GOBP_NOTCH_SIGNALING_PATHWAY,y = subtype,fill = subtype))+
  geom_density_ridges(scale = 2,size = 0.75,alpha = 0.75)+
  theme_bw()

#####


#############################tumor 50
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSVA)
library(GSEABase)
genelist <- getGmt('pathway/pt_tumor/h.all.v2022.1.Hs.symbols.gmt')
# genelist<-genelist[!(genelist$gene==''),]
# gcSample = split(genelist$gene,
#                  genelist$term)
#######subset cell
rcc.combined<-readRDS('data/combined/rcc_combined_4_73477.rds')
meta<-readRDS('data/combined/meta.rds')
meta<-meta[colnames(rcc.combined),]
rcc.combined$subtype<-meta$subtype_cancer

pttumor<-subset(rcc.combined,subtype %in% c('PROM1_PT-S1','PROM1_PT-S3','PROM1_PT-S2','PT-S1','PT-S2','PT-S3',
                                            'PT-S1/S2','type 1','type 2','type 3','type 4'))

######aucell
data<-pttumor@assays$RNA@data
library(AUCell)
cells_rankings <- AUCell_buildRankings(data, nCores = 1)
cells_AUC <- AUCell_calcAUC(genelist, cells_rankings,nCores = 1)
saveRDS(cells_AUC,'pathway/pt_tumor/cells_AUC_tumor50.rds')

x<- cells_AUC@assays@data@listData[["AUC"]]
x<-t(x)
meta<-pttumor@meta.data
meta<-cbind(meta,x)
pttumor@meta.data<-meta
DimPlot(pttumor,group.by = 'subtype')

######ridge plot
library(ggridges)
library(sciplot)
meta$subtype<-factor(meta$subtype,levels = c('PROM1_PT-S1','PROM1_PT-S2','PROM1_PT-S3','PT-S1','PT-S1/S2','PT-S2','PT-S3','type 1','type 2','type 3','type 4'),ordered = T)
ggplot(meta,aes(x = HALLMARK_PANCREAS_BETA_CELLS,y = subtype,fill = subtype))+
  geom_density_ridges(scale = 2,size = 0.75,alpha = 0.75)+
  theme_bw()

# [1] "HALLMARK_TNFA_SIGNALING_VIA_NFKB"           "HALLMARK_HYPOXIA"                          
# [3] "HALLMARK_CHOLESTEROL_HOMEOSTASIS"           "HALLMARK_MITOTIC_SPINDLE"                  
# [5] "HALLMARK_WNT_BETA_CATENIN_SIGNALING"        "HALLMARK_TGF_BETA_SIGNALING"               
# [7] "HALLMARK_IL6_JAK_STAT3_SIGNALING"           "HALLMARK_DNA_REPAIR"                       
# [9] "HALLMARK_G2M_CHECKPOINT"                    "HALLMARK_APOPTOSIS"                        
# [11] "HALLMARK_NOTCH_SIGNALING"                   "HALLMARK_ADIPOGENESIS"                     
# [13] "HALLMARK_ESTROGEN_RESPONSE_EARLY"           "HALLMARK_ESTROGEN_RESPONSE_LATE"           
# [15] "HALLMARK_ANDROGEN_RESPONSE"                 "HALLMARK_MYOGENESIS"                       
# [17] "HALLMARK_PROTEIN_SECRETION"                 "HALLMARK_INTERFERON_ALPHA_RESPONSE"        
# [19] "HALLMARK_INTERFERON_GAMMA_RESPONSE"         "HALLMARK_APICAL_JUNCTION"                  
# [21] "HALLMARK_APICAL_SURFACE"                    "HALLMARK_HEDGEHOG_SIGNALING"               
# [23] "HALLMARK_COMPLEMENT"                        "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"        
# [25] "HALLMARK_PI3K_AKT_MTOR_SIGNALING"           "HALLMARK_MTORC1_SIGNALING"                 
# [27] "HALLMARK_E2F_TARGETS"                       "HALLMARK_MYC_TARGETS_V1"                   
# [29] "HALLMARK_MYC_TARGETS_V2"                    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
# [31] "HALLMARK_INFLAMMATORY_RESPONSE"             "HALLMARK_XENOBIOTIC_METABOLISM"            
# [33] "HALLMARK_FATTY_ACID_METABOLISM"             "HALLMARK_OXIDATIVE_PHOSPHORYLATION"        
# [35] "HALLMARK_GLYCOLYSIS"                        "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"  
# [37] "HALLMARK_P53_PATHWAY"                       "HALLMARK_UV_RESPONSE_UP"                   
# [39] "HALLMARK_UV_RESPONSE_DN"                    "HALLMARK_ANGIOGENESIS"                     
# [41] "HALLMARK_HEME_METABOLISM"                   "HALLMARK_COAGULATION"                      
# [43] "HALLMARK_IL2_STAT5_SIGNALING"               "HALLMARK_BILE_ACID_METABOLISM"             
# [45] "HALLMARK_PEROXISOME"                        "HALLMARK_ALLOGRAFT_REJECTION"              
# [47] "HALLMARK_SPERMATOGENESIS"                   "HALLMARK_KRAS_SIGNALING_UP"                
# [49] "HALLMARK_KRAS_SIGNALING_DN"                 "HALLMARK_PANCREAS_BETA_CELLS" 