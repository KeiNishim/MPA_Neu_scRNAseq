library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(patchwork)
library(BiocManager)
library(data.table)
library(tidyr)
library(monocle3)
library(SeuratWrappers)

#Read BD Rhapsody output files (sample_info includes file_path, disease name and sample name)
#################################################
load_and_prep_ADT <- function(x){
  counts <- read.table( x['path'], skip=6, sep = ",", header = TRUE, row.names = 1)
  x.ab <- counts[,1:38,]
  x.adt <- CreateSeuratObject(counts = t(x.ab), assay ="ADT")
  x.adt <- AddMetaData(x.adt, metadata = x['disease'], col.name="disease")
  x.adt <- AddMetaData(x.adt, metadata = x['sample'], col.name="sample")
  x.adt <- NormalizeData(x.adt, normalization.method="CLR", margin = 2) ## normalization across cells 
  return(x.adt)}
samples_adt <- apply(sample_info, 1, load_and_prep_ADT)


load_and_prep_RNA <- function(x){
  counts <- read.table( x['path'], skip=6, sep = ",", header = TRUE, row.names = 1)
  x.rn <- counts[,40:length(counts),]
  x.rna <- CreateSeuratObject(counts = t(x.rn), min.features=100)
  x.rna <- AddMetaData(x.rna, metadata = x['disease'], col.name="disease")
  x.rna <- AddMetaData(x.rna, metadata = x['sample'], col.name="sample")
  ###############################################
  return(x.rna)}
samples_rna <- apply(sample_info, 1, load_and_prep_RNA)


#Integration for ADT
features_ADT <- rownames(samples_adt[[1]])
for(i in 1:14){
  features_ADT <- intersect(features_ADT, rownames(samples_adt[[i]]))
}
samples_adt <- lapply(X=samples_adt, FUN=function(x){
  x <- ScaleData(x, features=features_ADT, verbose=FALSE)
  x <- RunPCA(x, features=features_ADT, verbose=FALSE)
})
adt.anchors <- FindIntegrationAnchors(object.list = samples_adt, 
                                      anchor.features=features_ADT,reduction="rpca",ref=c(1,3,5,7,11),dims=1:12)
s.int <- IntegrateData(anchorset = adt.anchors, dims=1:12)

#Integration for RNA (SCTransform)
samples_rna <- lapply(X=samples_rna, FUN=NormalizeData)
samples_rna <- lapply(X=samples_rna, FUN=SCTransform)
samples_rna <- lapply(X=samples_rna, FUN=RunPCA)
features_RNA <- SelectIntegrationFeatures(object.list = samples_rna, assay=rep('SCT', length(samples_rna)), nfeatures = 3000)
samples_rna <- PrepSCTIntegration(object.list = samples_rna, anchor.features = features_RNA)
rna.anchors <- FindIntegrationAnchors(object.list = samples_rna, normalization.method='SCT', 
                                      anchor.features = features_RNA, reduction="cca", dims=1:20)
s.int <- IntegrateData(anchorset = rna.anchors, normalization.method="SCT")
#Then integrate RNA and ADT file as s.int


#Run PCA analysis
#######################################
DefaultAssay(s.int) <- "integrated"
s.int <- RunPCA(s.int, verbose=TRUE)
s.int <- FindNeighbors(s.int, reduction="pca")
s.int <- RunUMAP(s.int, reduction="pca", dims=1:15)

#Then manually remove doublet and dead cells

pp <-DimPlot(subset(s.int, idents=c("HD","MPA"),downsample=74406), group.by="seurat_clusters",split.by="disease")
pp
ggsave('Fig.1b.png',pp, width=8.5, height=4.8)
write.csv(table(s.int@meta.data$sample, s.int@meta.data$seurat_clusters),"MPAall_population.csv")

all <- FindAllMarkers(subset(s.int,downsample=50000),only.pos = TRUE,
                      logfc.threshold = 0.25, min.pct=0.25, max.cells.per.ident=5000)

top10 <- all %>% group_by(cluster) %>% top_n(n = -4, wt = p_val_adj)
top10 <- top10 %>% group_by(cluster) %>% top_n(n =4, wt = avg_log2FC)
X <- unique(top10['gene'])
Y <- apply(X,2,rev)
DefaultAssay(s.int) <- "SCT"
pp <- DotPlot(s.int,features=Y,scale.by="size", scale=T,col.max=5, col.min=-5)+
  coord_flip()+ theme(axis.text.x = element_text(angle=90,size = 8,vjust=0.5,hjust=1))+ 
  theme(axis.text.y = element_text(size = 7.5))+scale_size(range=c(0.1,3))+
  scale_colour_gradient2(low="blue", mid="white", high="red", midpoint=0)
pp
ggsave('Fig.1c.png',pp, width=7.5, height=7.8)

DefaultAssay(s.int) <- "ADT"
pp <- DotPlot(s.int,dot.min=0,features=rev(c(
  "CD11b.ICRF44.ITGAM.AHS0184.pAbO",
  "CD15.FUT4.AHS0196.pAbO",
  "CD183.CXCR3.AHS0031.pAbO",
  "CD11c.B.LY6.ITGAX.AHS0056.pAbO",
  "CD14.MPHIP9.CD14.AHS0037.pAbO",
  "CD16.3G8.FCGR3A.AHS0053.pAbO",
  "CD62L.DREG.56.SELL.AHS0049.pAbO",
  "CD193.CCR3.AHS0159.pAbO",
  "CD134.ACT35.TNFRSF4.AHS0013.pAbO",
  "FCER1A.FCER1A.AHS0129.pAbO",
  "HLA.DR.CD74.AHS0035.pAbO",
  "CD3.UCHT1.CD3E.AHS0231.pAbO",
                               "CD4.SK3.CD4.AHS0032.pAbO",
                               "CD28.L293.CD28.AHS0138.pAbO",
                               "CD25.2A3.IL2RA.AHS0026.pAbO",
                               "CD27.M.T271.CD27.AHS0025.pAbO",
                               "CD127.IL7R.AHS0028.pAbO",
                               "CD45RA.HI100.PTPRC.AHS0009.pAbO",
                               "CD278.ICOS.AHS0012.pAbO",
                               "CD8.SK1.CD8A.AHS0228.pAbO",
                               "CD279.EH12.1.PDCD1.AHS0014.pAbO",
                               "CD161.HP.3G10.KLRB1.AHS0205.pAbO",
                               "CD56.NCAM16.2.NCAM1.AHS0019.pAbO",
                               "GITR.TNFRSF18.AHS0104.pAbO",
                               "CD19.SJ25C1.CD19.AHS0030.pAbO",
                               "CD22.CD22.AHS0195.pAbO",
                               "CXCR5.CXCR5.AHS0039.pAbO",
                               "CD196.CCR6.AHS0034.pAbO",
                               "CD137.TNFRSF9.AHS0003.pAbO",
                               "CD272.BTLA.AHS0052.pAbO"
                               
                               
)),
              scale.by="size", scale=T,col.max=5, col.min=-5)+
  coord_flip()+ theme(axis.text.x = element_text(angle=90,size = 8,vjust=0.5,hjust=1))+ 
  theme(axis.text.y = element_text(size = 7.5))+scale_size(range=c(0.1,3))+
  scale_colour_gradient2(low="blue", mid="white", high="red", midpoint=0)
pp
ggsave('Fig.1d.png',pp, width=10, height=7.5)

pp <-DimPlot(subset(s.int, idents=rownames(table(s.int@meta.data[["sample"]])),downsample=7395), group.by="seurat_clusters",split.by="sample")
pp
ggsave('FigS1b.png',pp, width=40, height=4.8)

pp <-DimPlot(s.int,group.by="seurat_clusters", reduction="umap")
pp
ggsave('FigS1a.png',pp, width=6, height=5)

s.markerssct <- FindMarkers(s.int,ident.1=c("MPA"), ident.2=c("HD"),slot="data", 
                            only.pos = T, min.pct=0.1,logfc.threshold = 0.25, max.cells.per.ident=50000)
write.csv(s.markerssct,"MPADEgenes.csv")

#Analysis of Neutrophil
pp <-DimPlot(subset(s.int,idents=c("HD","MPA"),downsample=50000), group.by="seurat_clusters",split.by="disease", reduction="umap")
pp
ggsave('Fig2A.png',pp, width=8.5, height=4.8)

pp <-DimPlot(s.int, group.by="seurat_clusters")
pp
ggsave('FigS4-1.png',pp, width=5.5, height=4.8)

pp <-DimPlot(subset(s.int, idents=rownames(table(s.int@meta.data[["sample"]])),downsample=2727), group.by="seurat_clusters",split.by="sample", reduction="umap")
pp
ggsave('FigS4-2.png',pp, width=40, height=4.8)

s.int <- subset(s.int, idents=c("HD","MPA"), downsample=50000)
pp <- FeaturePlot(s.int, reduction="umap",features = c("MPO","MMP9"), max.cutoff="q99",ncol=3, split.by="disease")
pp
ggsave('Figs5b-1.png',pp, width=8, height=9)

pp <- FeaturePlot(s.int, reduction="umap",features = c("PADI4","ITGAM"), max.cutoff="q99",ncol=3, split.by="disease")
pp
ggsave('Figs5b-2.png',pp, width=8, height=9)

pp <- FeaturePlot(s.int, reduction="umap",features = c("GBP1","GBP5"),max.cutoff="q99",ncol=3, split.by="disease")
pp
ggsave('Figs5b-3.png',pp, width=8, height=9)

pp <- FeaturePlot(s.int, reduction="umap",features = c("FCGR1A","FCGR2A"),max.cutoff="q99.5",ncol=3, split.by="disease")
pp
ggsave('Figs5b-4.png',pp, width=8, height=9)

pp <- FeaturePlot(s.int, reduction="umap",features = c("FCGR3A","FCGR3B"),max.cutoff="q99",ncol=3, split.by="disease")
pp
ggsave('Figs5b5.png',pp, width=8, height=9)

#Analysis of monocyte
pp <-DimPlot(subset(s.int,idents=c("HD","MPA"),downsample=4695),reduction="umap", group.by="seurat_clusters",split.by="disease")
pp
ggsave('Figs7-1.png',pp, width=8.5, height=4.8)

pp <-DimPlot(s.int,reduction="umap", group.by="seurat_clusters",split.by="sample")
pp
ggsave('Figs7-4.png',pp, width=40, height=4.8)

all <- FindAllMarkers(s.int,only.pos = TRUE,
                      logfc.threshold = 0.25, min.pct=0.25, max.cells.per.ident=5000)
top10 <- all %>% group_by(cluster) %>% top_n(n = -4, wt = p_val_adj)
top10 <- top10 %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC)
X <- unique(top10['gene'])
Y <- apply(X,2,rev)
pp <- DotPlot(s.int,features=Y,scale.by="size", scale=T,col.max=5, col.min=-5)+
  coord_flip()+ theme(axis.text.x = element_text(angle=90,size = 8,vjust=0.5,hjust=1))+ 
  theme(axis.text.y = element_text(size = 7.5))+scale_size(range=c(0.1,3))+
  scale_colour_gradient2(low="blue", mid="white", high="red", midpoint=0)
pp
ggsave('FigS7-1.png',pp, width=5, height=5)

pp <- DotPlot(s.int,features=rev(c( "CD11b.ICRF44.ITGAM.AHS0184.pAbO","CD14.MPHIP9.CD14.AHS0037.pAbO","CD62L.DREG.56.SELL.AHS0049.pAbO",
                                   "CD11c.B.LY6.ITGAX.AHS0056.pAbO",
                                   "CD16.3G8.FCGR3A.AHS0053.pAbO","HLA.DR.CD74.AHS0035.pAbO",
                                   "FCER1A.FCER1A.AHS0129.pAbO"
)),
scale.by="size", scale=T,col.max=5, col.min=-5)+
  coord_flip()+ theme(axis.text.x = element_text(angle=90,size = 8,vjust=0.5,hjust=1))+ 
  theme(axis.text.y = element_text(size = 7.5))+scale_size(range=c(0.1,3))+
  scale_colour_gradient2(low="blue", mid="white", high="red", midpoint=0)
pp
ggsave('FigS7-3.png',pp, width=6.5, height=4)

pp <- DotPlot(s.int,features=rev(c( "CD15.FUT4.AHS0196.pAbO",
                                    "CD62L.DREG.56.SELL.AHS0049.pAbO",
                                    "CD16.3G8.FCGR3A.AHS0053.pAbO",
                                    "CD183.CXCR3.AHS0031.pAbO",
                                   "CD11b.ICRF44.ITGAM.AHS0184.pAbO",
                                   "CD11c.B.LY6.ITGAX.AHS0056.pAbO",
                                   "HLA.DR.CD74.AHS0035.pAbO",
                                   "CD278.ICOS.AHS0012.pAbO",
                                   "CXCR5.CXCR5.AHS0039.pAbO"
)),
scale.by="size", scale=T,col.max=5, col.min=-5)+
  coord_flip()+ theme(axis.text.x = element_text(angle=90,size = 8,vjust=0.5,hjust=1))+ 
  theme(axis.text.y = element_text(size = 7.5))+scale_size(range=c(0.1,2))+
  scale_colour_gradient2(low="blue", mid="white", high="red", midpoint=0)
pp
ggsave('FigS4-3.png',pp, width=6.5, height=2.38)
