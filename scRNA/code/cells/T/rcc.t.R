
setwd('E:\\work\\RCC\\immune\\T')


library(Seurat)
library(ggplot2)
rcc.combined<-readRDS('/public/home/xwtang/work/RCC/data/intergrated/rcc_combined_4.rds')
rcc.t<-subset(rcc.combined,celltype.refined=='Lymphocyte')
DefaultAssay(rcc.t)<-'integrated'

#plot
FeaturePlot(rcc.t,'PTPRC')
VlnPlot(rcc.t,c('PTPRC','CD247','EPCAM'))
DimPlot(rcc.t,label = T)

rcc.t<-FindVariableFeatures(rcc.t,nFeatures=2000)
rcc.t <- ScaleData(rcc.t, verbose = FALSE)
rcc.t <- RunPCA(rcc.t, npcs = 30, verbose = FALSE)
rcc.t <- RunUMAP(rcc.t, reduction = "pca", dims = 1:15)
rcc.t <- FindNeighbors(rcc.t, reduction = "pca", dims = 1:15)
rcc.t <- FindClusters(rcc.t, resolution = 1)

#plot
DimPlot(rcc.t,label = T)
DimPlot(rcc.t,group.by = 'orig.ident')
FeaturePlot(rcc.t,'CD4')

###
cg=c('CA9','SLC31A2','CD68','CD163','CD247','CD8A','CD4','GNLY','GZMB','KIT','SLC18A2','PTPRC','PECAM1','PDGFRB','NPHS1','NPHS2','EMCN','CLEC4C','JCHAIN','IRF8','UBE2C','TOP2A',
     'CTSS','CD19','MS4A1','FCGR3A','NCAM1','KDR','FLT1','VCAM1','CFH','NPHS1','NPHS2')
DotPlot(rcc.t, features = unique(cg),
        assay='RNA'  )  + coord_flip()+ RotatedAxis() 


ms=c('CD4','ANXA1','GNLY','TCF7','CXCR6','CXCR5','GZMK','IL23R','CXCL13','FOXP3',
     'CTLA4','CD8A','LEF1','GPR183','CX3CR1','KLRG1','FGFBP2','PRF1','CD6','CD160','LAYN','SLC4A10','FCGR3A','NCAM1','KLRB1')
DotPlot(rcc.t, features = unique(ms),
        assay='RNA'  )  + coord_flip()+ RotatedAxis() 
FeaturePlot(rcc.t,features = ms)

#plot cd8 gzmk
cd8gzmk<-c('GZMK','CXCR4','CXCR3','CD44')
DotPlot(rcc.t, features = unique(cd8gzmk),
        assay='RNA'  )  + coord_flip()+ RotatedAxis() 
FeaturePlot(rcc.t,features = cd8gzmk)


#cd8
cd8<-c('CD3E','CD8A','CCR7','SELL','FASLG','CD44','FAS','PDCD1','UBE2C','TOP2A')
pdf('/public/home/xwtang/work/RCC/immune_cell/T/pic/cd8_D_marker.pdf',width = 8,height = 6)
DotPlot(rcc.t, features = unique(cd8),
        assay='RNA'  )  + coord_flip()+ RotatedAxis() 
FeaturePlot(rcc.t,features = cd8)
dev.off()

#cd8 temra
cd8emra<-c('S1PR1','S1PR5','ITGB7','CX3CR1','KLRG1','FGFBP2')
DotPlot(rcc.t, features = unique(cd8emra),
        assay='RNA'  )  + coord_flip()+ RotatedAxis() 
FeaturePlot(rcc.t,features = cd8emra)

#cd8 layn
DefaultAssay(rcc.t)<-'RNA'
cd8lyan<-c('HAVCR2','CXCL13','PDCD1','LAYN','TOX','IFNG','GZMB','GZMA','CTLA4','TIGIT','LAG3')
DotPlot(rcc.t, features = unique(cd8lyan))  + coord_flip()+ RotatedAxis() 
FeaturePlot(rcc.t,features = cd8lyan)

#cd8 tem
cd8cd6<-c('EOMES','SUB1','CD74','GZMK','CST7','DUSP2','CMC1','DKK3')
FeaturePlot(rcc.t,features = cd8cd6)


#####
pdf('/public/home/xwtang/work/RCC/immune_cell/T/pic/test.pdf',width = 8,height = 6)
FeaturePlot(rcc.t,'TRA')
VlnPlot(rcc.t,c('PTPRC','TRA'))
DimPlot(rcc.t,label = T)
dev.off()

########
#0 CD8+ Tex
#1 CD8+ Tem cd44/gzmk
#2 CD4+ Tn lef1/sell/tcf7/cd28
#3 CD8+ Temra  cx3cr1/klrg1/fgfbp2/prf1  S1PR1,S1PR5,ITGB7  Temra
#4 T-cancer cell
#5 T-cancer cell
#6 CD4+ Tm cd69/il7r/
#7 NKdim cd56dim cd16high(fcgr3a)/cd56(ncam1)/gnly
#8 NKbright cd56bright nkg2a(KLRC1) cd16(fcgr3a)/cd56high(ncam1)/gnly   
#9  cd8+Tem cd44/gzmk/il7r
#10 NK-cancer cell
#11 Treg ctla4/foxp3
#12 CD8+ Tex proliferative top2a

#0 CD8+ Trm
#1 CD8+ Tex
#2 CD4+ Tn
#3 CD8+ Temra  cx3cr1/klrg1/fgfbp2/prf1  S1PR1,S1PR5,ITGB7  Temra
#4 NKdim cd56dim cd16high(fcgr3a)/cd56(ncam1)/gnly
#5 T-EPC
#6 T-EPC
#7 T-EPC
#8 Treg
#9 NKbright cd56bright nkg2a(KLRC1) cd16(fcgr3a)/cd56high(ncam1)/gnly  
#10 CD8+ Tex proliferative top2a
#11 B cell




type<-c('0'='CD8+ Trm','1'='CD8+ Tex','2'='CD4+ Tn','3'='CD8+ Temra','4'='NKdim',
        '5'='T-EPC','6'='T-EPC','7'='T-EPC','8'='Treg','9'=
          'NKbright','10'='CD8+ T(proliferating)','11'='B cell')
rcc.t[['celltype.refined']] = unname(type[rcc.t@meta.data$seurat_clusters])

saveRDS(rcc.t,"E:\\work\\RCC\\immune\\T\\rcc.t.rds")


####
DefaultAssay(rcc.t)<-'RNA'
FeaturePlot(rcc.t,c('LRP2','CUBN','CA9','CD247','CD4','CD8A','PTPRC','FCGR3A','NCAM1'))

rcc.t$celltype.refined<-factor(rcc.t$celltype.refined,levels = c('NKbright','NKdim','CD4+ Tn','CD4+ Tm','Treg','CD8+ Temra','CD8+ Tem','CD8+ Tex','CD8+ Tex(proliferative)','NK-EPI','T-EPI'),ordered = T)
cg=c('FCGR3A','NCAM1','CD4','LEF1','IL7R','FOXP3','CTLA4','CX3CR1','KLRG1','GZMK','HAVCR2','TOP2A',
     'CA9','LRP2','CUBN')
DotPlot(rcc.t, features = unique(cg),
        assay='RNA',group.by = 'celltype.refined'  )  + coord_flip()+ RotatedAxis() 

FeaturePlot(rcc.t,c('FCGR3A','NCAM1'))

cd4<-subset(rcc.t, celltype.refined %in% c('CD4+ Tm','CD4+ Tn','Treg'))
VlnPlot(cd4, c('LEF1','TCF7','CD28'),group.by = 'celltype.refined',pt.size = 0)
VlnPlot(cd4, c('S100A4','KLRB1','GZMK'),group.by = 'celltype.refined',pt.size = 0)
VlnPlot(cd4, c('TNFRSF18','CTLA4','LAYN'),group.by = 'celltype.refined',pt.size = 0)

cd8<-subset(rcc.t, celltype.refined %in% c('CD8+ Tem','CD8+ Temra','CD8+ Tex','CD8+ Tex(proliferative)'))
VlnPlot(cd8, c('CX3CR1','FGFBP2','GZMH'),group.by = 'celltype.refined',pt.size = 0)
VlnPlot(cd8, c('GZMK','CXCR4','CD44'),group.by = 'celltype.refined',pt.size = 0)
VlnPlot(cd8, c('HAVCR2','PDCD1','TNFRSF9'),group.by = 'celltype.refined',pt.size = 0)
VlnPlot(cd8, c('UBE2C','TOP2A','MKI67'),group.by = 'celltype.refined',pt.size = 0,ncol = 2)
VlnPlot(cd8, c('CX3CR1','FGFBP2','GZMH','GZMK','CXCR4','CD44','HAVCR2','PDCD1','TNFRSF9','UBE2C','TOP2A','MKI67'),group.by = 'celltype.refined',pt.size = 0,ncol = 2)

nk<-subset(rcc.t,celltype.refined %in% c('NKdim','NKbright'))
VlnPlot(nk, c('FCGR3A','NCAM1'),group.by = 'celltype.refined',pt.size = 0)


####样本信息
type<-c('normala'='T1a','normald'='T3a','hsma'='T3a','hsmb'='T1a','RCCa'='T1a','RCCd'='T3a')
rcc.t[['stage']] = unname(type[rcc.t@meta.data$orig.ident])

meta<-rcc.t@meta.data
pie(as.vector(table(meta$stage)),labels = names(table(meta$stage)),main = '总体分布') 
meta1<-subset(meta,stage=='T1a')
pie(as.vector(table(meta1$celltype.refined)),labels = names(table(meta1$celltype.refined)),main = 'T1a细胞分布')
meta2<-subset(meta,stage=='T3a')
pie(as.vector(table(meta2$celltype.refined)),labels = names(table(meta2$celltype.refined)),main = 'T3a细胞分布')
metaex<-subset(meta,celltype.refined=='CD8+ Tex')
pie(as.vector(table(metaex$stage)),labels = names(table(metaex$stage)),main = 'Tex分布')
metatreg<-subset(meta,celltype.refined=='Treg')
pie(as.vector(table(metatreg$stage)),labels = names(table(metatreg$stage)),main = 'Treg分布')
metaepc<-subset(meta,celltype.refined=='T-EPC')
pie(as.vector(table(metaepc$stage)),labels = names(table(metaepc$stage)),main = 'T-EPC分布')
metarm<-subset(meta,celltype.refined=='CD8+ Trm')
pie(as.vector(table(metarm$stage)),labels = names(table(metarm$stage)),main = 'Trm分布')
metaemra<-subset(meta,celltype.refined=='CD8+ Temra')
pie(as.vector(table(metaemra$stage)),labels = names(table(metaemra$stage)),main = 'Temra分布')


####pheatmap
rcc.t<-readRDS('monocle/T/rcc.t_1666.rds')
rcc.t$celltype.refined<-ifelse(rcc.t$celltype.refined=='CD4+ Tn','CD4+ Tconv',rcc.t$celltype.refined)
rcc.t$celltype.refined<-ifelse(rcc.t$celltype.refined=='Treg','CD4+ Treg',rcc.t$celltype.refined)
rcc.t$celltype.refined<-ifelse(rcc.t$celltype.refined=='CD8+ Temra','CD8+ Teff',rcc.t$celltype.refined)
rcc.t$celltype.refined<-ifelse(rcc.t$celltype.refined=='CD8+ Trm','CD8+ Tmem',rcc.t$celltype.refined)
rcc.t$celltype.refined<-ifelse(rcc.t$celltype.refined=='CD8+ T(proliferating)','CD8+ Tprolif',rcc.t$celltype.refined)

rcc.t<-subset(rcc.t,celltype.refined %in% c('CD4+ Tconv','CD4+ Treg','CD8+ Teff','CD8+ Tmem','CD8+ Tex','CD8+ Tprolif'))
rcc.t$celltype.refined<-factor(rcc.t$celltype.refined,levels = c('CD4+ Tconv','CD4+ Treg','CD8+ Teff','CD8+ Tmem','CD8+ Tex','CD8+ Tprolif'),ordered = T)
Idents(rcc.t)<-rcc.t$celltype.refined
markers<-FindAllMarkers(rcc.t,only.pos = T)
markers<-subset(markers,p_val_adj<0.05)
markers<-subset(markers,avg_log2FC>0.7)

data<-rcc.t@assays$RNA@data
data<-data.frame(data[markers$gene,])
df<-data.frame(matrix(0,nrow = nrow(data),ncol = 6))
meta<-rcc.t@meta.data
colnames(df)<-c('CD4+ Tconv','CD4+ Treg','CD8+ Teff','CD8+ Tmem','CD8+ Tex','CD8+ Tprolif')
rownames(df)<-rownames(data)
for(i in c('CD4+ Tconv','CD4+ Treg','CD8+ Teff','CD8+ Tmem','CD8+ Tex','CD8+ Tprolif')){
  tmp<-subset(meta,celltype.refined==i)
  df[,i]<-rowMeans(data[,gsub('-','.',rownames(tmp))])
}


library(ComplexHeatmap)
df<-t(scale(t(df)))

library("BuenColors")
col <- jdb_color_maps[1:25]    #选取了25个颜色
celltype_info<-factor(c('CD4+ Tconv','CD4+ Treg','CD8+ Teff','CD8+ Tmem','CD8+ Tex','CD8+ Tprolif'),
                      levels = c('CD4+ Tconv','CD4+ Treg','CD8+ Teff','CD8+ Tmem','CD8+ Tex','CD8+ Tprolif'),ordered = T)
names(col) <- levels(celltype_info)
top_anno <- HeatmapAnnotation(
        cluster = anno_block(gp = gpar(fill = col), # 设置填充色
        labels = levels(celltype_info), 
        labels_gp = gpar(cex = 0.8, col = "white"))) # 设置字体
mark_gene<-c('TGFBR', 'IL7R', 'THEMIS', 'DOCK9', 'BTBD11', 'PDE4B', 'KALRN', 'APBA2', 'GRAP2',  
         'IKZF2', 'TOX', 'PLCL1', 'IL2RA', 'PELI1','TNFRSF4','ARHGAP31','SKAP1','NAMPT','CADM1',
         'TOX', 'TOX2','ETV1', 'CCDC141', 'CXCL13', 'TRPS1', 'CRTAM', 'BIRC3',  'TNFRSF1B', 'GZMK', 'TNFRSF9', 'CRIM1','HAVCR2',
         'CX3CR1', 'GNLY', 'PRF1', 'GZMH', 'GZMB', 'A2M', 'KLRC1','KLRC3', 'KLRB1', 'KLRG1','KLRF1','KLRD1','FGFBP2','MYBL1', 
          'CAMK4', 'CD44', 'RASGRF2', 'CD69', 'ITGA1', 
          'TOP2A',  'ASPM' , 'MKI67',  'CENPF', 'CENPP', 'CENPE', 'HMGB2', 'HMGB1','KIF20B', 'TUBA1B' )

gene_pos <- which(rownames(df) %in% unique(mark_gene))

row_anno <- rowAnnotation(mark_gene = anno_mark(at = gene_pos, 
                                                labels = rownames(df)[gene_pos],
                                                labels_gp=gpar(fontsize = 8)))
#修改颜色
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("#2A5A9A", "white",'#BF3939'))
p<-Heatmap(df,
        col = col_fun,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        column_split = celltype_info,
        top_annotation = top_anno,
        right_annotation = row_anno,
        column_title = NULL,
        heatmap_legend_param = list(
                title = "Z-Score"
        ))
p 

