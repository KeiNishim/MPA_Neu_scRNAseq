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
library(dittoSeq)

#For Heatmap analysis###
#Use seurat object from MPA clustering.R
all <- FindAllMarkers(subset(s.int,downsample=20000),only.pos = TRUE,
                          logfc.threshold = 0.25, min.pct=0.25, max.cells.per.ident=5000)

write.csv(all, "Genemarkers.csv")

top10 <- all %>% group_by(cluster) %>% top_n(n = -4, wt = p_val_adj)
top10 <- top10 %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC)
X <- unique(top10['gene'])
Y <- apply(X,2,rev)
pp <- DotPlot(s.int,features=Y,scale.by="size", scale=T,col.max=5, col.min=-5)+
  coord_flip()+ theme(axis.text.x = element_text(angle=90,size = 8,vjust=0.5,hjust=1))+ 
  theme(axis.text.y = element_text(size = 7.5))+scale_size(range=c(0.1,3))+
  scale_colour_gradient2(low="blue", mid="white", high="red", midpoint=0)
ggsave('Fig2B.png',pp, width=5.5, height=5.5)

pp <- dittoHeatmap(subset(s.int, idents=c("AI143","AI193","AI54","AI55","AI62","AI75","HD"),
                    downsample=1000),
             genes=rownames(s.markerssct),
             annot.by=c("disease","sample2","seurat_clusters"),
             order.by=c("disease","sample2","seurat_clusters"),
             scale.to.max=T,
             slot="data",
             border_color="black")
pp
ggsave('Fig.4A.png', pp, width=15, height=11.2)

pp <- dittoHeatmap(subset(s.int, idents=c("AI143","AI193","AI54","AI55","AI62","AI75","HD"),
                          downsample=1000),
                   genes=c("MMP9","PADI4","RFLNB","ALOX5AP","MME","HMGB2","FCN1","VIM","APMAP","CDA","QPCT","PLP2","CYP4F3","CKAP4"),
                   annot.by=c("disease","sample2","seurat_clusters"),
                   order.by=c("disease","sample2","seurat_clusters"),
                   
                   scale.to.max=T,
                   slot="data",
                   border_color="black")
pp
ggsave('Fig.S9b.png', pp, width=10, height=3)

pp <- dittoHeatmap(subset(s.int, idents=c("AI143","AI193","AI54","AI55","AI62","AI75","HD"),
                          downsample=1000),
                   genes=c("APOL6",
                           "CARD16",
                           "EPSTI1",
                           "FCGR1A",
                           "GBP1",
                           "GBP5",
                           "IFI16",
                           "IFIT3",
                           "IFITM1",
                           "IFITM3",
                           "PARP9",
                           "PLSCR1",
                           "RNF213",
                           "SAMD9L",
                           "SERPING1",
                           "TNFAIP6",
                           "TNFSF13B",
                           "TRIM22",
                           "XAF1"),
                   annot.by=c("disease","sample2","seurat_clusters"),
                   order.by=c("disease","sample2","seurat_clusters"),
                   
                   scale.to.max=T,
                   slot="data",
                   border_color="black")
pp
ggsave('Fig.s9a.png', pp, width=10, height=3.8)

#Enrichment score analysis
features <- list(c("MPO","AZU1","ELANE","LTF","TOP2A","CEACAM8","CEBPE"))
s.int <- AddModuleScore(s.int, features=features, ctrl=1,name="IFNA.score", slot="data")
p <- FeaturePlot(s.int, "IFNA.score1",pt.size=0.5)+
  scale_colour_gradient2(low="blue", mid="white", high="red", midpoint = 0.5)
p
ggsave('Fig.s4-1.png', p, width=4, height=4)

features <- list(c("P2RY14","STEAP4","GBP4","GBP5","CD274","GBP1","GK","ANKRD22","CSF2RB","GBP2","LCP2",
                   "LRRK2","WARS","PRR5L","GCLM","LPCAT2","CD69","FCGR1A"))
s.int <- AddModuleScore(s.int, features=features, ctrl=1,name="IFNA.score", slot="data")
p <- FeaturePlot(s.int, "IFNA.score1",pt.size=0.5)+
  scale_colour_gradient2(low="blue", mid="white", high="red", midpoint = 0.25)
p
ggsave('Fig.s4-2.png', p, width=4, height=4)

features <- list(c("HERC5","IFIT1","RSAD2","ISG15","OASL","MX1","TRIM22"))
s.int <- AddModuleScore(s.int, features=features, ctrl=1,name="IFNA.score", slot="data")
p <- FeaturePlot(s.int, "IFNA.score1",pt.size=0.5)+
  scale_colour_gradient2(low="blue", mid="white", high="red", midpoint = 0.25)
p
ggsave('Fig.s4-3.png', p, width=4, height=4)

features <- list(c("PI3","SOD2","SLPI","SLC7A11","TNFAIP6","IL1A","CD22","SERPINB9","CCL2","IL1B"))
s.int <- AddModuleScore(s.int, features=features, ctrl=1,name="IFNA.score", slot="data")
p <- FeaturePlot(s.int, "IFNA.score1",pt.size=0.5)+
  scale_colour_gradient2(low="blue", mid="white", high="red", midpoint = 0.5)
p
ggsave('Fig.s4-4.png', p, width=4, height=4)



features <- list(c("A1BG","ABCA13","ACAA1","ACLY","ACP3","ACTR10","ACTR1B","ACTR2","ADA2","ADAM10","ADAM8","ADGRE3","ADGRE5","ADGRG3","AGA","AGL","AGPAT2","AHSG","ALAD","ALDH3B1","ALDOA","ALDOC","ALOX5","AMPD3","ANO6","ANPEP","ANXA2","AOC1","AP1M1","AP2A2","APAF1","APEH","APRT","ARG1","ARHGAP45","ARHGAP9","ARL8A","ARMC8","ARPC5","ARSA","ARSB","ASAH1","ATAD3B","ATG7","ATP11A","ATP11B","ATP6AP2","ATP6V0A1","ATP6V0C","ATP6V1D","ATP8A1","ATP8B4","AZU1","B2M","B4GALT1","BIN2","BPI","BRI3","BST1","BST2","C1orf35","C3","C3AR1","C5AR1","C6orf120","CAB39","CALML5","CAMP","CAND1","CANT1","CAP1","CAPN1","CAT","CCT2","CCT8","CD14","CD177","CD300A","CD33","CD36","CD44","CD47","CD53","CD55","CD58","CD59","CD63","CD68","CD93","CDA","CDK13","CEACAM1","CEACAM3","CEACAM6","CEACAM8","CEP290","CFD","CFP","CHI3L1","CHIT1","CHRNB4","CKAP4","CLEC12A","CLEC4C","CLEC4D","CLEC5A","CMTM6","CNN2","COMMD3","COMMD9","COPB1","COTL1","CPNE1","CPNE3","CPPED1","CR1","CRACR2A","CREG1","CRISP3","CRISPLD2","CSNK2B","CST3","CSTB","CTSA","CTSB","CTSC","CTSD","CTSG","CTSH","CTSS","CTSZ","CXCL1","CXCR1","CXCR2","CYB5R3","CYBA","CYBB","CYFIP1","CYSTM1","DBNL","DDOST","DDX3X","DEFA1","DEFA1B","DEFA4","DEGS1","DERA","DGAT1","DIAPH1","DNAJC13","DNAJC3","DNAJC5","DNASE1L1","DOCK2","DOK3","DPP7","DSC1","DSG1","DSN1","DSP","DYNC1H1","DYNC1LI1","DYNLL1","DYNLT1","EEF1A1","EEF2","ELANE","ENPP4","EPX","ERP44","FABP5","FAF2","FCAR","FCER1G","FCGR2A","FCGR3B","FCN1","FGL2","FGR","FLG2","FOLR3","FPR1","FPR2","FRK","FRMPD3","FTH1","FTL","FUCA1","FUCA2","GAA","GALNS","GCA","GDI2","GGH","GHDC","GLA","GLB1","GLIPR1","GM2A","GMFG","GNS","GOLGA7","GPI","GPR84","GRN","GSDMD","GSN","GSTP1","GUSB","GYG1","HBB","HEBP2","HEXB","HGSNAT","HK3","HLA-B","HLA-C","HMGB1","HMOX2","HP","HPSE","HRNR","HSP90AA1","HSP90AB1","HSPA1A","HSPA1B","HSPA6","HSPA8","HUWE1","HVCN1","IDH1","IGF2R","ILF2","IMPDH1","IMPDH2","IQGAP1","IQGAP2","IRAG2","IST1","ITGAL","ITGAM","ITGAV","ITGAX","ITGB2","JUP","KCMF1","KCNAB2","KPNB1","KRT1","LAIR1","LAMP1","LAMP2","LAMTOR1","LAMTOR2","LAMTOR3","LCN2","LGALS3","LILRA3","LILRB2","LILRB3","LPCAT1","LRG1","LRRC7","LTA4H","LTF","LYZ","MAGT1","MAN2B1","MANBA","MAPK1","MAPK14","MCEMP1","METTL7A","MGAM","MGST1","MIF","MLEC","MME","MMP25","MMP8","MMP9","MNDA","MOSPD2","MPO","MS4A3","MVP","NAPRT","NBEAL2","NCKAP1L","NCSTN","NDUFC2","NEU1","NFAM1","NFASC","NFKB1","NHLRC3","NIT2","NME2","NPC2","NRAS","OLFM4","OLR1","ORM1","ORM2","ORMDL3","OSCAR","OSTF1","P2RX1","PA2G4","PADI2","PAFAH1B2","PDAP1","PDXK","PECAM1","PFKL","PGAM1","PGLYRP1","PGM1","PGM2","PGRMC1","PIGR","PKM","PKP1","PLAC8","PLAU","PLAUR","PLD1","PLEKHO2","PNP","PPBP","PPIA","PPIE","PRCP","PRDX4","PRDX6","PRG2","PRG3","PRKCD","PRSS2","PRSS3","PRTN3","PSAP","PSEN1","PSMA2","PSMA5","PSMB1","PSMB7","PSMC2","PSMC3","PSMD1","PSMD11","PSMD12","PSMD13","PSMD14","PSMD2","PSMD3","PSMD6","PSMD7","PTAFR","PTGES2","PTPN6","PTPRB","PTPRC","PTPRJ","PTPRN2","PTX3","PYCARD","PYGB","PYGL","QPCT","QSOX1","RAB10","RAB14","RAB18","RAB24","RAB27A","RAB31","RAB37","RAB3A","RAB3D","RAB44","RAB4B","RAB5B","RAB5C","RAB6A","RAB7A","RAB9B","RAC1","RAP1A","RAP1B","RAP2B","RAP2C","RETN","RHOA","RHOF","RHOG","RNASE2","RNASE3","RNASET2","ROCK1","S100A11","S100A12","S100A7","S100A8","S100A9","S100P","SCAMP1","SDCBP","SELL","SERPINA1","SERPINA3","SERPINB1","SERPINB10","SERPINB12","SERPINB3","SERPINB6","SIGLEC14","SIGLEC5","SIGLEC9","SIRPA","SIRPB1","SLC11A1","SLC15A4","SLC27A2","SLC2A3","SLC2A5","SLC44A2","SLCO4C1","SLPI","SNAP23","SNAP25","SNAP29","SPTAN1","SRP14","STBD1","STING1","STK10","STK11IP","STOM","SURF4","SVIP","SYNGR1","TARM1","TBC1D10C","TCIRG1","TCN1","TICAM2","TIMP2","TLR2","TMBIM1","TMC6","TMEM179B","TMEM30A","TMEM63A","TNFAIP6","TNFRSF1B","TOLLIP","TOM1","TRAPPC1","TRPM2","TSPAN14","TTR","TUBB","TUBB4B","TXNDC5","TYROBP","UBR4","UNC13D"
                   ,"VAMP8","VAPA","VAT1","VCL","VCP","VNN1","VPS35L","XRCC5","XRCC6","YPEL5"))         
s.int <- AddModuleScore(s.int, features=features, ctrl=1,name="DEG.score", slot="data")

features <- list(c("B2M","CAMK2A","CAMK2B","CAMK2D","CAMK2G","CD44","CIITA","FCGR1A","FCGR1BP","GBP1","GBP2","GBP3","GBP4","GBP5","GBP6","GBP7","HLA-A","HLA-B","HLA-C","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQA2","HLA-DQB1","HLA-DQB2","HLA-DRA","HLA-DRB1","HLA-DRB3","HLA-DRB4","HLA-DRB5","HLA-E","HLA-F","HLA-G","HLA-H","ICAM1","IFI30","IFNG","IFNGR1","IFNGR2","IRF1","IRF2","IRF3","IRF4","IRF5","IRF6","IRF7","IRF8","IRF9","JAK1","JAK2","MAPK1","MAPK3","MID1","MT2A","NCAM1","OAS1","OAS2","OAS3","OASL","PIAS1","PML","PRKCD","PTAFR","PTPN1","PTPN11","PTPN2","PTPN6","RAF1","SMAD7","SOCS1","SOCS3","SP100","STAT1","SUMO1","TRIM10","TRIM14","TRIM17","TRIM2","TRIM21","TRIM22","TRIM25","TRIM26","TRIM29","TRIM3","TRIM31","TRIM34","TRIM35","TRIM38","TRIM45","TRIM46","TRIM48","TRIM5","TRIM6","TRIM62","TRIM68","TRIM8","VCAM1","YBX1"))         
s.int <- AddModuleScore(s.int, features=features, ctrl=1,name="IFNGscore", slot="data")

features=list(c("ALPL",                "BST1",                "CD93",                "CEACAM3",
                "CREB5",                "CRISPLD2",                "CSF3R",                "CYP4F3",
                "DYSF",                "FCAR",                "FCGR3B",                "CPPED1",
                "FPR1",                "FPR2",                "G0S2",                "HIST1H2BC",
                "HPSE",                "CXCR1",                "CXCR2",                "KCNJ15",
                "LILRB2",                "MGAM",                "MME",                "PDE4B",
                "S100A12",              "SIGLEC5",                "SLC22A4",                "SLC25A37",
                "TECPR2",                "TNFRSF10C",                "VNN3",                "AKT1",
                "AKT2",                "ATG7",                "CLEC6A",                "CSF3",
                "CTSG",                "CYBB",                "DNASE1",                "ELANE",
                "ENTPD4",                "F3",                "HMGB1",                "IL17A",
                "IL1B",                "IL6",                "IL8",                "IRAK4",
                "ITGAM",                "ITGB2",                "KCNN3",                "MAPK1",
                "MAPK3",                "MMP9",                "MPO",                "MTOR",
                "PADI4",                "PTAFR",                "PIK3CA",                "RIPK1",
                "RIPK3",                "SELP",                "SELPLG",                "SIGLEC14",
                "TLR2",                "TLR4",                "TLR7",                "TLR8",
                "TNF"
))

features=list(c("FCGR3A","FCGR3b","SYK","MAP3k7","RAF1","MAP2K1","MAP2k2","MAPK1","MAPK3","SIGLEC9","CYBB","NCF1","CYBA","NCF2",
                "NCF4","RAC1","RAC2","TLR7","TLR8","ELANE","MPO","ACTG1","ACTB","VDAC1","VDAC2","VDAC3","SLC25A4",
                "SLC25A5","SLC25A6","SLC25A31","PPIF","PADI4","IGH","FCGR1A","FCGR2A","ITGAM","ITGB2",
                "ITGAL","CLEC7A","SRC","PLCB1","PLCB2","PLCB3","PLCB4","PLCG1","PLCG2","PRKCA","PRKCB","PRKCG","ATG7","FPR1",
              "FPR2","PIK3CA","PIK3CD","PIK3CB","PIK3R1","PIK3R2","PIK3R3","AKT1","AKT2","AKT3","MTOR","NFKB1",
              "RELA","C3","CR1","CR1L"," C5","C5AR","HMGB1","TLR2","TLR4","MAPK11","MAPK12","MAPK13","MAPK14",
              "ITGA2B","ITGB3","FGA","FGB","FGG","GP1BA","VWF","SELP","SELPG","AGER",
              "CASP4","CASP1","GSDMD","HAT1","AZU1","CTSG","CAMP","CLCN3","CLCN4","CLCN5","AQP9"))
s.int <- AddModuleScore(s.int, features=features, ctrl=1,name="NETscore", slot="data")

features=list(c("FCGR1A","FCGR2A","FCGR2B","FCGR2C","FCGR3A","FCGR3B","FCGRT","TRIM21","FCRL5"))
s.int <- AddModuleScore(s.int, features=features, ctrl=1,name="FCG.score", slot="data")

FeaturePlot(s.int, c("NET.score1"), pt.size=0.5)+
  scale_colour_gradient2(low="blue", mid="white", high="red", midpoint=0)#REACTOME NEUTROPHIL and 69 NET
pp <- DotPlot(s.int,features=c("NET.score1"),group.by="seurat_clusters",split.by="disease",
              scale.by="size", scale=T,col.max=5, col.min=-5,
              cols=c("lightgrey", "blue","red","red","red","red","red","red","red","red","red","red","red","red","red","red"))+
  coord_flip()+ theme(axis.text.x = element_text(angle=90,size = 8,vjust=0.5,hjust=1))+ 
  theme(axis.text.y = element_text(size = 7.5))+scale_size(range=c(0.1,3))+
  scale_colour_gradient2(low="blue", mid="white", high="red", midpoint=0)
write.csv(pp[["data"]],"SampleScoreNET.csv")
pp
ggsave('Fig.3C.png', p, width=7, height=3.2)
