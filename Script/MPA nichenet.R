library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(patchwork)
library(BiocManager)
library(metap)
library(data.table)
library(BiocGenerics)
library("org.Hs.eg.db")
library("hgu133plus2.db")
library(jetset)
library(icellnet)
library(gridExtra)
library(devtools)
library(circlize)
library(nichenetr)
library(tidyverse)

##Use seurat object from MPA clustering.R
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
s.int <- RenameIdents(s.int, "Neutrophil_1"="Neutrophil", "Neutrophil_2"="Neutrophil","Myelocyte"="Neutrophil", "CD4 Naive"="CD4Naive","CD4 Memory"="CD4Memory",
                      "CD8 Naive/TCM"="CD8NaiveTCM", "CD8 TEM/CTL"="CD8TEMCTL","CD14 Mono"="CD14Mono","CD16 Mono"="CD16Mono"
                      )
DimPlot(s.neu)
Ligand =c("Neutrophil", "CD4Naive","CD4Memory","CD8NaiveTCM","CD8TEMCTL","NK","Bcell","PBPC","CD14Mono","CD16Mono","cDC","pDC")
Receptor = c("T2IFN")
expressed_genes_sender = rownames(FindMarkers(s.int, ident.1=Ligand, ident.2=c("Neutrophil"),
                                              only.pos = FALSE, min.pct=0.0050,logfc.threshold = -Inf, max.cells.per.ident=5000))
expressed_genes_receiver = rownames(FindMarkers(s.neu, ident.1=Receptor,only.pos = FALSE, min.pct=0.1,logfc.threshold = -Inf, max.cells.per.ident=2000))
interested_genes_receiver <- rownames(FindMarkers(s.neu, ident.1=Receptor,logfc.threshold=0.25, min.pct=0.1, only.pos=T, max.cells.per.ident=2000))
geneset_oi = interested_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)] 
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
ligands= lr_network %>% pull(from) %>% unique()
ligands <-intersect(colnames(ligand_target_matrix),ligands)#
expressed_ligands = intersect(ligands,expressed_genes_sender)
receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_receiver)
lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors)
potential_ligands = lr_network_expressed %>% pull(from) %>% unique()

ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities %>% arrange(-pearson) 
write.csv(ligand_activities,"T2IFNligands.csv")
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
best_upstream_ligands

# show histogram of ligand activity scores
p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(20, pearson) %>% pull(pearson))), color="red", linetype="dashed", linewidth=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
active_ligand_target_links_df <- na.omit(active_ligand_target_links_df)
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix)
nrow(active_ligand_target_links_df)
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
order_targets <- intersect(order_targets, rownames(active_ligand_target_links))
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

#Circus output
avg_expression_ligands = AverageExpression(s.int, features = best_upstream_ligands, assay="RNA")
sender_ligand_assignment = avg_expression_ligands$RNA %>% apply(1, function(ligand_expression){
  ligand_expression > (ligand_expression %>% mean() + ligand_expression %>% sd())
}) %>% t()
sender_ligand_assignment = sender_ligand_assignment %>% apply(2, function(x){x[x == TRUE]}) %>% purrr::keep(function(x){length(x) > 0})
names(sender_ligand_assignment)
all_assigned_ligands = sender_ligand_assignment %>% lapply(function(x){names(x)}) %>% unlist()
unique_ligands = all_assigned_ligands %>% table() %>% .[. == 1] %>% names()
general_ligands = best_upstream_ligands %>% setdiff(unique_ligands)

PBPC_specific_ligands = sender_ligand_assignment$PBPC %>% names() %>% setdiff(general_ligands)
NK_specific_ligands = sender_ligand_assignment$NK %>% names() %>% setdiff(general_ligands)
Neutrophil_specific_ligands = sender_ligand_assignment$Neutrophil %>% names() %>% setdiff(general_ligands)
CD4Naive_specific_ligands = sender_ligand_assignment$`CD4Naive` %>% names() %>% setdiff(general_ligands)
CD4Memory_specific_ligands = sender_ligand_assignment$`CD4Memory` %>% names() %>% setdiff(general_ligands)
CD8NaiveTCM_specific_ligands = sender_ligand_assignment$`CD8NaiveTCM` %>% names() %>% setdiff(general_ligands)
CD8TEMCTL_specific_ligands = sender_ligand_assignment$`CD8TEMCTL` %>% names() %>% setdiff(general_ligands)
Bcell_specific_ligands = sender_ligand_assignment$Bcell %>% names() %>% setdiff(general_ligands)
CD14Mono_specific_ligands = sender_ligand_assignment$`CD14Mono` %>% names() %>% setdiff(general_ligands)
CD16Mono_specific_ligands = sender_ligand_assignment$`CD16Mono` %>% names() %>% setdiff(general_ligands)
cDC_specific_ligands = sender_ligand_assignment$`cDC` %>% names() %>% setdiff(general_ligands)
pDC_specific_ligands = sender_ligand_assignment$`pDC` %>% names() %>% setdiff(general_ligands)

ligand_type_indication_df = tibble(
  ligand_type = c(rep("PBPC-specific", times = PBPC_specific_ligands %>% length()),
                  rep("NK-specific", times = NK_specific_ligands %>% length()),
                  rep("Neutrophil-specific", times = Neutrophil_specific_ligands %>% length()),
                  rep("CD4Naive-specific", times = CD4Naive_specific_ligands %>% length()),
                  rep("CD4Memory-specific", times = CD4Memory_specific_ligands %>% length()),
                  rep("CD8NaiveTCM-specific", times = CD8NaiveTCM_specific_ligands %>% length()),
                  rep("CD8TEMCTL-specific", times = CD8TEMCTL_specific_ligands %>% length()),
                  rep("Bcell-specific", times = Bcell_specific_ligands %>% length()),
                  rep("CD14Mono-specific", times = CD14Mono_specific_ligands %>% length()),
                  rep("CD16Mono-specific", times = CD16Mono_specific_ligands %>% length()),
                  rep("cDC-specific", times = cDC_specific_ligands %>% length()),
                  rep("pDC-specific", times = pDC_specific_ligands %>% length()),
                  rep("General", times = general_ligands %>% length())),
  
  ligand = c(PBPC_specific_ligands, NK_specific_ligands, Neutrophil_specific_ligands,
             CD4Naive_specific_ligands,CD4Memory_specific_ligands,CD8NaiveTCM_specific_ligands,CD8TEMCTL_specific_ligands,
             Bcell_specific_ligands,CD14Mono_specific_ligands,CD16Mono_specific_ligands,cDC_specific_ligands,
             pDC_specific_ligands,general_ligands))


active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
active_ligand_target_links_df = active_ligand_target_links_df %>% mutate(target_type = "T2IFN") %>% inner_join(ligand_type_indication_df) # if you want ot make circos plots for multiple gene sets, combine the different data frames and differentiate which target belongs to which gene set via the target type
active_ligand_target_links_df <-active_ligand_target_links_df[-150,]

cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(0.67)
active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)
ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())
circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)

grid_col_ligand =c("General" = "lawngreen",
                   "NK-specific" = "royalblue",
                   "Neutrophil-specific" = "darkgreen",
                   "CD4Naive-specific" = "steelblue",
                   "CD4Memory-specific" = "steelblue4",
                   "CD8NaiveTCM-specific"="orange",
                   "CD8TEMCTL-specific"="orange4",
                   "Bcell-specific"="violetred",
                   "PBPC-specific"="violetred4",
                   "CD14Mono-specific"="lightpink",
                   "CD16Mono-specific"="lightpink4",
                   "cDC-specific"="khaki",
                   "pDC-specific"="khaki4")
grid_col_target =c(
  "T2IFN" = "tomato")

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target = tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
library(dplyr)
links_circle = circos_links %>%select(ligand,target, weight)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(target,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$target)
grid_col =c(grid_ligand_color,grid_target_color)

# give the option that links in the circos plot will be transparant ~ ligand-target potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 
target_order = circos_links$target %>% unique()
ligand_order = c(general_ligands,
                 CD4Naive_specific_ligands,CD4Memory_specific_ligands,CD8NaiveTCM_specific_ligands,CD8TEMCTL_specific_ligands,
                 NK_specific_ligands, Bcell_specific_ligands,PBPC_specific_ligands, 
                 CD14Mono_specific_ligands,CD16Mono_specific_ligands,cDC_specific_ligands,pDC_specific_ligands,
                 Neutrophil_specific_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)

order = c(ligand_order,target_order)
width_same_cell_same_ligand_type = 0.5
width_different_cell = 1
width_ligand_target = 1
width_same_cell_same_target_type = 0.5
table(circos_links[,"ligand_type"])

gaps = c(
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "General") %>% distinct(ligand) %>% nrow() -1)),width_different_cell,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "CD4Naive-specific") %>% distinct(ligand) %>% nrow() -1)),width_different_cell,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "CD4Memory-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "CD8NaiveTCM-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "CD8TEMCTL-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "NK-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Bcell-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "PBPC-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "CD14Mono-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "CD16Mono-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "cDC-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "pDC-specific") %>% distinct(ligand) %>% nrow() -1)),width_ligand_target,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Neutrophil-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_type == "T2IFN") %>% distinct(target) %>% nrow() -1)),width_ligand_target
)
circos.clear()
circos.par(gap.degree = gaps)

chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
}, bg.border = NA) #
#Save plot image

#Immature CCI analysis
Ligand =c("Neutrophil", "CD4Naive","CD4Memory","CD8NaiveTCM","CD8TEMCTL","NK","Bcell","PBPC","CD14Mono","CD16Mono","cDC","pDC")
Receptor = c("Immature")
expressed_genes_sender = rownames(FindMarkers(s.int, ident.1=Ligand, ident.2=c("Neutrophil"),
                                              only.pos = FALSE, min.pct=0.005,logfc.threshold = -Inf, max.cells.per.ident=5000))
expressed_genes_receiver = rownames(FindMarkers(s.neu, ident.1=Receptor,only.pos = FALSE, min.pct=0.1,logfc.threshold = -Inf, max.cells.per.ident=2000))
interested_genes_receiver <- rownames(FindMarkers(s.neu, ident.1=Receptor,logfc.threshold=0.25, min.pct=0.1, only.pos=T, max.cells.per.ident=2000))
geneset_oi = interested_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)] 
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
ligands = lr_network %>% pull(from) %>% unique()
ligands <-intersect(colnames(ligand_target_matrix),ligands)#
expressed_ligands = intersect(ligands,expressed_genes_sender)
receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_receiver)
lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors)
potential_ligands = lr_network_expressed %>% pull(from) %>% unique()

ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities %>% arrange(-pearson)
write.csv(ligand_activities,"Immatureligands.csv")
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
best_upstream_ligands

# show histogram of ligand activity scores
p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(20, pearson) %>% pull(pearson))), color="red", linetype="dashed", linewidth=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
active_ligand_target_links_df <- na.omit(active_ligand_target_links_df)
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix)
nrow(active_ligand_target_links_df)
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
order_targets <- intersect(order_targets, rownames(active_ligand_target_links))
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

#Circus output
avg_expression_ligands = AverageExpression(s.int, features = best_upstream_ligands, assay="RNA")
sender_ligand_assignment = avg_expression_ligands$RNA %>% apply(1, function(ligand_expression){
  ligand_expression > (ligand_expression %>% mean() + ligand_expression %>% sd())
}) %>% t()
sender_ligand_assignment = sender_ligand_assignment %>% apply(2, function(x){x[x == TRUE]}) %>% purrr::keep(function(x){length(x) > 0})
names(sender_ligand_assignment)
all_assigned_ligands = sender_ligand_assignment %>% lapply(function(x){names(x)}) %>% unlist()
unique_ligands = all_assigned_ligands %>% table() %>% .[. == 1] %>% names()
general_ligands = best_upstream_ligands %>% setdiff(unique_ligands)

PBPC_specific_ligands = sender_ligand_assignment$PBPC %>% names() %>% setdiff(general_ligands)
NK_specific_ligands = sender_ligand_assignment$NK %>% names() %>% setdiff(general_ligands)
Neutrophil_specific_ligands = sender_ligand_assignment$Neutrophil %>% names() %>% setdiff(general_ligands)
CD4Naive_specific_ligands = sender_ligand_assignment$`CD4Naive` %>% names() %>% setdiff(general_ligands)
CD4Memory_specific_ligands = sender_ligand_assignment$`CD4Memory` %>% names() %>% setdiff(general_ligands)
CD8NaiveTCM_specific_ligands = sender_ligand_assignment$`CD8NaiveTCM` %>% names() %>% setdiff(general_ligands)
CD8TEMCTL_specific_ligands = sender_ligand_assignment$`CD8TEMCTL` %>% names() %>% setdiff(general_ligands)
Bcell_specific_ligands = sender_ligand_assignment$Bcell %>% names() %>% setdiff(general_ligands)
CD14Mono_specific_ligands = sender_ligand_assignment$`CD14Mono` %>% names() %>% setdiff(general_ligands)
CD16Mono_specific_ligands = sender_ligand_assignment$`CD16Mono` %>% names() %>% setdiff(general_ligands)
cDC_specific_ligands = sender_ligand_assignment$`cDC` %>% names() %>% setdiff(general_ligands)
pDC_specific_ligands = sender_ligand_assignment$`pDC` %>% names() %>% setdiff(general_ligands)

ligand_type_indication_df = tibble(
  ligand_type = c(rep("PBPC-specific", times = PBPC_specific_ligands %>% length()),
                  rep("NK-specific", times = NK_specific_ligands %>% length()),
                  rep("Neutrophil-specific", times = Neutrophil_specific_ligands %>% length()),
                  rep("CD4Naive-specific", times = CD4Naive_specific_ligands %>% length()),
                  rep("CD4Memory-specific", times = CD4Memory_specific_ligands %>% length()),
                  rep("CD8NaiveTCM-specific", times = CD8NaiveTCM_specific_ligands %>% length()),
                  rep("CD8TEMCTL-specific", times = CD8TEMCTL_specific_ligands %>% length()),
                  rep("Bcell-specific", times = Bcell_specific_ligands %>% length()),
                  rep("CD14Mono-specific", times = CD14Mono_specific_ligands %>% length()),
                  rep("CD16Mono-specific", times = CD16Mono_specific_ligands %>% length()),
                  rep("cDC-specific", times = cDC_specific_ligands %>% length()),
                  rep("pDC-specific", times = pDC_specific_ligands %>% length()),
                  rep("General", times = general_ligands %>% length())),
  
  ligand = c(PBPC_specific_ligands, NK_specific_ligands, Neutrophil_specific_ligands,
             CD4Naive_specific_ligands,CD4Memory_specific_ligands,CD8NaiveTCM_specific_ligands,CD8TEMCTL_specific_ligands,
             Bcell_specific_ligands,CD14Mono_specific_ligands,CD16Mono_specific_ligands,cDC_specific_ligands,
             pDC_specific_ligands,general_ligands))


active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n =250) %>% bind_rows()
active_ligand_target_links_df = active_ligand_target_links_df %>% mutate(target_type = "Immature") %>% inner_join(ligand_type_indication_df) # if you want ot make circos plots for multiple gene sets, combine the different data frames and differentiate which target belongs to which gene set via the target type
active_ligand_target_links_df
active_ligand_target_links_df <-active_ligand_target_links_df[-170,]

cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(0.75)
active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)
ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())
circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)

grid_col_ligand =c("General" = "lawngreen",
                   "NK-specific" = "royalblue",
                   "Neutrophil-specific" = "darkgreen",
                   "CD4Naive-specific" = "steelblue",
                   "CD4Memory-specific" = "steelblue4",
                   "CD8NaiveTCM-specific"="orange",
                   "CD8TEMCTL-specific"="orange4",
                   "Bcell-specific"="violetred",
                   "PBPC-specific"="violetred4",
                   "CD14Mono-specific"="lightpink",
                   "CD16Mono-specific"="lightpink4",
                   "cDC-specific"="khaki",
                   "pDC-specific"="khaki4")
grid_col_target =c(
  "Immature" = "tomato")

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target = tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
library(dplyr)
links_circle = circos_links %>% select(ligand,target, weight)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(target,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$target)
grid_col =c(grid_ligand_color,grid_target_color)

# give the option that links in the circos plot will be transparant ~ ligand-target potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 
target_order = circos_links$target %>% unique()
ligand_order = c(general_ligands,
                 CD4Naive_specific_ligands,CD4Memory_specific_ligands,CD8NaiveTCM_specific_ligands,CD8TEMCTL_specific_ligands,
                 NK_specific_ligands, Bcell_specific_ligands,PBPC_specific_ligands, 
                 CD14Mono_specific_ligands,CD16Mono_specific_ligands,cDC_specific_ligands,pDC_specific_ligands,
                 Neutrophil_specific_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,target_order)
width_same_cell_same_ligand_type = 0.5
width_different_cell = 1
width_ligand_target = 1
width_same_cell_same_target_type = 0.5
table(circos_links[,"ligand_type"])

gaps = c(
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "General") %>% distinct(ligand) %>% nrow() -1)),width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "CD4Naive-specific") %>% distinct(ligand) %>% nrow() -1)),width_different_cell,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "CD4Memory-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "CD8NaiveTCM-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "CD8TEMCTL-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "NK-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Bcell-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "PBPC-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "CD14Mono-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "CD16Mono-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "cDC-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "pDC-specific") %>% distinct(ligand) %>% nrow() -1)),width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Neutrophil-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_type == "Immature") %>% distinct(target) %>% nrow() -1)),width_ligand_target
)
circos.clear()
circos.par(gap.degree = gaps)

chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
}, bg.border = NA) #

#HD v.s. MPA
Ligand =c("MPA")
Receptor = c("MPA")
expressed_genes_sender = rownames(FindMarkers(s.int, ident.1=Ligand,
                                              only.pos = FALSE, min.pct=0.005,logfc.threshold = -Inf, max.cells.per.ident=5000))
expressed_genes_receiver = rownames(FindMarkers(s.neu, ident.1=Receptor,only.pos = FALSE, min.pct=0.1,logfc.threshold = -Inf, max.cells.per.ident=2000))
interested_genes_receiver <- rownames(FindMarkers(s.neu, ident.1=Receptor,logfc.threshold=0.25, min.pct=0.1, only.pos=T, max.cells.per.ident=2000))
geneset_oi = interested_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)] 
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
ligands = lr_network %>% pull(from) %>% unique()
ligands <-intersect(colnames(ligand_target_matrix),ligands)#
expressed_ligands = intersect(ligands,expressed_genes_sender)
receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_receiver)
lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors)
potential_ligands = lr_network_expressed %>% pull(from) %>% unique()

ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities %>% arrange(-pearson) 
write.csv(ligand_activities,"T2IFNligands.csv")
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)

# show histogram of ligand activity scores
p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(20, pearson) %>% pull(pearson))), color="red", linetype="dashed", linewidth=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
active_ligand_target_links_df <- na.omit(active_ligand_target_links_df)
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix)
nrow(active_ligand_target_links_df)
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
order_targets <- intersect(order_targets, rownames(active_ligand_target_links))
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

#Circus output
Idents(s.int) <- "seurat_clusters"
avg_expression_ligands = AverageExpression(s.int, features = best_upstream_ligands, assay="RNA")
sender_ligand_assignment = avg_expression_ligands$RNA %>% apply(1, function(ligand_expression){
  ligand_expression > (ligand_expression %>% mean() + ligand_expression %>% sd())
}) %>% t()
sender_ligand_assignment = sender_ligand_assignment %>% apply(2, function(x){x[x == TRUE]}) %>% purrr::keep(function(x){length(x) > 0})
names(sender_ligand_assignment)
all_assigned_ligands = sender_ligand_assignment %>% lapply(function(x){names(x)}) %>% unlist()
unique_ligands = all_assigned_ligands %>% table() %>% .[. == 1] %>% names()
general_ligands = best_upstream_ligands %>% setdiff(unique_ligands)

PBPC_specific_ligands = sender_ligand_assignment$PBPC %>% names() %>% setdiff(general_ligands)
NK_specific_ligands = sender_ligand_assignment$NK %>% names() %>% setdiff(general_ligands)
Neutrophil_specific_ligands = sender_ligand_assignment$Neutrophil %>% names() %>% setdiff(general_ligands)
CD4Naive_specific_ligands = sender_ligand_assignment$`CD4Naive` %>% names() %>% setdiff(general_ligands)
CD4Memory_specific_ligands = sender_ligand_assignment$`CD4Memory` %>% names() %>% setdiff(general_ligands)
CD8NaiveTCM_specific_ligands = sender_ligand_assignment$`CD8NaiveTCM` %>% names() %>% setdiff(general_ligands)
CD8TEMCTL_specific_ligands = sender_ligand_assignment$`CD8TEMCTL` %>% names() %>% setdiff(general_ligands)
Bcell_specific_ligands = sender_ligand_assignment$Bcell %>% names() %>% setdiff(general_ligands)
CD14Mono_specific_ligands = sender_ligand_assignment$`CD14Mono` %>% names() %>% setdiff(general_ligands)
CD16Mono_specific_ligands = sender_ligand_assignment$`CD16Mono` %>% names() %>% setdiff(general_ligands)
cDC_specific_ligands = sender_ligand_assignment$`cDC` %>% names() %>% setdiff(general_ligands)
pDC_specific_ligands = sender_ligand_assignment$`pDC` %>% names() %>% setdiff(general_ligands)

ligand_type_indication_df = tibble(
  ligand_type = c(rep("PBPC-specific", times = PBPC_specific_ligands %>% length()),
                  rep("NK-specific", times = NK_specific_ligands %>% length()),
                  rep("Neutrophil-specific", times = Neutrophil_specific_ligands %>% length()),
                  rep("CD4Naive-specific", times = CD4Naive_specific_ligands %>% length()),
                  rep("CD4Memory-specific", times = CD4Memory_specific_ligands %>% length()),
                  rep("CD8NaiveTCM-specific", times = CD8NaiveTCM_specific_ligands %>% length()),
                  rep("CD8TEMCTL-specific", times = CD8TEMCTL_specific_ligands %>% length()),
                  rep("Bcell-specific", times = Bcell_specific_ligands %>% length()),
                  rep("CD14Mono-specific", times = CD14Mono_specific_ligands %>% length()),
                  rep("CD16Mono-specific", times = CD16Mono_specific_ligands %>% length()),
                  rep("cDC-specific", times = cDC_specific_ligands %>% length()),
                  rep("pDC-specific", times = pDC_specific_ligands %>% length()),
                  rep("General", times = general_ligands %>% length())),
  
  ligand = c(PBPC_specific_ligands, NK_specific_ligands, Neutrophil_specific_ligands,
             CD4Naive_specific_ligands,CD4Memory_specific_ligands,CD8NaiveTCM_specific_ligands,CD8TEMCTL_specific_ligands,
             Bcell_specific_ligands,CD14Mono_specific_ligands,CD16Mono_specific_ligands,cDC_specific_ligands,
             pDC_specific_ligands,general_ligands))


active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
active_ligand_target_links_df = active_ligand_target_links_df %>% mutate(target_type = "T2IFN") %>% inner_join(ligand_type_indication_df) # if you want ot make circos plots for multiple gene sets, combine the different data frames and differentiate which target belongs to which gene set via the target type
#active_ligand_target_links_df <-active_ligand_target_links_df[-115,]

cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(0.67)
active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)
ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())
circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)

grid_col_ligand =c("General" = "lawngreen",
                   "NK-specific" = "royalblue",
                   "Neutrophil-specific" = "darkgreen",
                   "CD4Naive-specific" = "steelblue",
                   "CD4Memory-specific" = "steelblue4",
                   "CD8NaiveTCM-specific"="orange",
                   "CD8TEMCTL-specific"="orange4",
                   "Bcell-specific"="violetred",
                   "PBPC-specific"="violetred4",
                   "CD14Mono-specific"="lightpink",
                   "CD16Mono-specific"="lightpink4",
                   "cDC-specific"="khaki",
                   "pDC-specific"="khaki4")
grid_col_target =c(
  "T2IFN" = "tomato")

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target = tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
library(dplyr)
links_circle = circos_links %>% select(ligand,target, weight)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(target,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$target)
grid_col =c(grid_ligand_color,grid_target_color)

# give the option that links in the circos plot will be transparant ~ ligand-target potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 
target_order = circos_links$target %>% unique()
ligand_order = c(general_ligands,
                 CD4Naive_specific_ligands,CD4Memory_specific_ligands,CD8NaiveTCM_specific_ligands,CD8TEMCTL_specific_ligands,
                 NK_specific_ligands, Bcell_specific_ligands,PBPC_specific_ligands, 
                 CD14Mono_specific_ligands,CD16Mono_specific_ligands,cDC_specific_ligands,pDC_specific_ligands,
                 Neutrophil_specific_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,target_order)
width_same_cell_same_ligand_type = 0.5
width_different_cell = 1
width_ligand_target = 1
width_same_cell_same_target_type = 0.5
table(circos_links[,"ligand_type"])

gaps = c(
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "General") %>% distinct(ligand) %>% nrow() -1)),width_different_cell,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "CD4Naive-specific") %>% distinct(ligand) %>% nrow() -1)),width_different_cell,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "CD4Memory-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "CD8NaiveTCM-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "CD8TEMCTL-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "NK-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Bcell-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "PBPC-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "CD14Mono-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "CD16Mono-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "cDC-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "pDC-specific") %>% distinct(ligand) %>% nrow() -1)),width_ligand_target,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Neutrophil-specific") %>% distinct(ligand) %>% nrow() -1)), width_ligand_target,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_type == "T2IFN") %>% distinct(target) %>% nrow() -1)),width_ligand_target
)
circos.clear()
circos.par(gap.degree = gaps)

chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
}, bg.border = NA) #
#Save plot image
