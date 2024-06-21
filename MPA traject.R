library(velocyto.R)


#SeuratWrappers####
cds <- as.cell_data_set(s.neu)
cds <- cluster_cells(cds, reduction_method="UMAP")
cds <- learn_graph(cds, use_partition=T, close_loop=F, verbose=T, learn_graph_control=list(ncenter=1200))
cds <- order_cells(cds)
pp <- plot_cells(cds=cds, 
           color_cells_by="seurat_clusters",
           label_cell_group=F,
           label_groups_by_cluster=F,
           show_trajectory_graph=T,
           label_leaves=F,
           label_branch_points=F,
           label_roots =F,
           cell_size=0.1)
pp
ggsave('Fig3B.png',pp, width=5.8, height=4.8)

pp <- plot_cells(cds=cds, 
           color_cells_by="pseudotime",
           label_groups_by_cluster=F,
           show_trajectory_graph=T,
           label_leaves=F,
           label_branch_points=F,
           label_roots = F,
           label_principal_points=F,
           cell_size=0.1)
pp
ggsave('Fig3C.png',pp, width=5.8, height=4.8)
saveRDS(s.neu, "MPAcds.rds")

pt.matrix <- exprs(cds)[match(Y,rownames(rowData(cds))),order(pseudotime(cds))]
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes;
#K means with 6 groups
htkm <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 6,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)
htkm