library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(Seurat)

#Use seurat object from MPA clustering.R
pbmc_small_sce <- as.SingleCellExperiment(s.int, dimreducs = c("pca","umap"), graphs=c("integrated_snn", "integrated_nn"), assay="integrated")
traj_milo <- Milo(pbmc_small_sce)
plotUMAP(traj_milo, colour_by="disease", point_size=0.1)
traj_milo <- buildGraph(traj_milo, k = 10, d = 15, reduced.dim="PCA")
traj_milo <- makeNhoods(traj_milo, prop = 0.02, k = 10, d=15,reduced_dim="PCA", refined = TRUE)#0.1, 10, 30
traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), sample="sample")
traj_milo <- calcNhoodDistance(traj_milo, d=15, reduced.dim="PCA")
saveRDS(traj_milo, "Neu_milo3.rds")

traj_design <- data.frame(colData(traj_milo))[,c("sample", "disease")]
traj_design <- distinct(traj_design)
age <- data.frame(c(79, 80, 69, 70, 82, 79, 64, 62, 70, 62, 70, 59, 84))

names(age) <- "age"
traj_design <- cbind(traj_design, age)
rownames(traj_design) <- traj_design$sample
traj_design$class <- c("2_MPA","2_MPA","2_MPA","2_MPA","2_MPA","2_MPA",  "1_HD", "1_HD", "1_HD", "1_HD", "1_HD", "1_HD", "1_HD")
da_results <- testNhoods(traj_milo,reduced.dim="PCA", design = ~ age+class, design.df = traj_design, )
da_results %>%
  arrange(- SpatialFDR) %>%
  head()
traj_milo <- buildNhoodGraph(traj_milo)
ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)
ggplot(da_results, aes(logFC, -log10(SpatialFDR)))+ 
  geom_point() +
  geom_hline(yintercept = 1)

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(traj_milo, da_results,size_range=c(0.5,3), layout="UMAP",alpha=0.05)+#3 and 0.05 for pbmc 
  scale_fill_gradient2(low='blue', mid='white', high="red")
nh_graph_pl

ggsave('Fig.1f.png', nh_graph_pl +plot_layout(guides="collect"), width=7.5, height=7.5)
ggsave('Fig.S9-1.png', nh_graph_pl +plot_layout(guides="collect"), width=7.5, height=7.5)
ggsave('Fig.2e.png', nh_graph_pl +plot_layout(guides="collect"), width=7.5, height=6.5)

da_results <- annotateNhoods(traj_milo, da_results, coldata_col = "seurat_clusters")
da_results$seurat_clusters <- factor(da_results$seurat_clusters,
                                     levels= (c("PI","T2IFN","T1IFN","Aged","Mature","Immature","Myelocyte")))#Neutrophil
da_results$seurat_clusters <- factor(da_results$seurat_clusters,
                                     levels= rev(c("CD4 Naive","CD4 Memory","CD8 Naive/TCM","CD8 TEM/CTL","NK",
                                                "Bcell","PBPC",
                                                "Neutrophil_1","Neutrophil_2",
                                                "CD14 Mono","CD16 Mono","Eosino/Baso","Myelocyte","cDC","pDC")))#WBC
da_results$seurat_clusters <- factor(da_results$seurat_clusters,
                                     levels= rev(c("CD14Mono_Activated","CD14Mono_VCAN","CD14Mono_HLA","CD16Mono","cDC")))
pDa <- plotDAbeeswarm(da_results, group.by = "seurat_clusters", alpha=0.05)+
  scale_color_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" )+
  geom_hline(yintercept=0, linetype="dashed")+geom_point(size=0.5)+geom_boxplot(outlier.shape = NA)
pDa
ggsave('Fig.1g.png', pDa, width=9, height=15)
ggsave('Fig.S9-2.png', pDa, width=12, height=15)
ggsave('Fig.2f.png', pDa, width=8.5, height=10)

median(da_results[da_results$seurat_clusters=="T2IFN", "logFC"])
