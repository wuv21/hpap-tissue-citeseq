# note that this code is written to be run from the project base directory
# renv::load("/data/hpap-citeseq/hpap-citeseq-analysis")

source("figures/genericFigureSettings.R")
source("scripts/dimPlots.R")
source("scripts/deg.R")
library(ComplexHeatmap)
library(Seurat)
library(cowplot)
library(presto)

set.seed(42)

tsaCatalog <- readRDS("rds/tsa_catalog.rds")

seu <- tryCatch(
  {
    print(head(seu$TissueCondensed))
    message("Seurat object already exists")
    return(seu)
  },
  error = function(cond) {
    message("Seurat object doesn't exist. Loading now.")
    tmp <- readRDS("outs/rds/seuMergedPostHSP_forFigures_2023-09-17_09-03-10.rds")
    
    return(tmp)
  })

manualClusterOrder <- unique(seu$manualAnnot)
manualClusterOrder <- factor(manualClusterOrder,
  levels = customSortAnnotation(manualClusterOrder),
  labels = stringr::str_trim(customSortAnnotation(manualClusterOrder)))

seu$Disease_Status <- factor(seu$Disease_Status, levels = c("ND", "AAb+", "T1D"))


################################################################################
# S4A - dot plots!
################################################################################
baseMarkers <- c(
  "A0034", #CD3
  "A0072", #CD4
  "A0046", #CD8
  "A0081", #CD14
  "A0083", #CD16
  "A0047", #CD56
  "A0050", #CD19
  "A0087", #CD45RO
  "A0063", #CD45RA
  "A0154", #CD27
  "A0386", #CD28
  "A0156", #CD95
  "A0390", #CD127
  "A0146", #CD69
  "A0085", #CD25
  "A0089", #TIGIT
  "A0088", #PD1
  "A0141", #CCR5
  "A0144", #CXCR5
  "A0149"  #CD161
)

labellerAdt <- tsaCatalog$cleanName
names(labellerAdt) <- paste0("adt_", tsaCatalog$DNA_ID)

# base adt dot style
figABase <- DotPlot(seu, features = paste0("adt_", baseMarkers), group.by = "manualAnnot")

figA <- figABase +
  scale_y_discrete(limits = levels(manualClusterOrder)) +
  textSizeOnlyTheme +
  theme(
    panel.grid.major.y = element_line(color = "#eeeeee"),
    axis.text.x = element_text(angle = 45, size = BASEPTFONTSIZE, hjust = 1)) +
  scale_x_discrete(labels = labellerAdt)


basePanelRNA <- c(
  "IL7R", # cd127
  "CCR7",
  "CD14",
  "CSF3R", #cd114
  "LYZ",
  "S100A4",
  "MS4A1",
  "CD8A",
  "FCGR3A",
  "MS4A7",
  "GNLY",
  "NKG7",
  "FCER1A",
  "CST3",
  "PPBP", #cxcl7 (platelet basic protein)
  "CXCR5",
  "TOX",
  "TCF7",
  "TBX21",
  "CD69",
  "CD3D",
  "CD3G",
  "CD40LG",
  "ICOS",
  "FOXP3",
  "GATA3",
  "IKZF2", # helios
  "RORC",
  "BCL6",
  "PRDM1", #blimp
  "CD38",
  "TNFRSF17" #BCMA
)

# base rna dot style
figBBase <- DotPlot(seu, features = basePanelRNA, group.by = "manualAnnot")

figB <- figBBase +
  scale_y_discrete(limits = levels(manualClusterOrder)) +
  textSizeOnlyTheme +
  theme(
    panel.grid.major.y = element_line(color = "#eeeeee"),
    axis.text.x = element_text(angle = 45, size = BASEPTFONTSIZE, hjust = 1))


################################################################################
# Final layout and plot all
################################################################################
layout <- c(
  patchwork::area(1, 1, 4, 12), # a
  patchwork::area(5, 1, 8, 12) # b
)

p <- wrap_elements(plot = figA) +
  wrap_elements(plot = figB) +
  plot_annotation(tag_levels = list(LETTERS[1:2])) +
  plot_layout(design = layout) &
  plotTagTheme

saveFinalFigure(
  plot = p,
  prefixDir = "figures/outs",
  fn = "sfig4_final",
  devices = c("pdf", "png"),
  addTimestamp = TRUE,
  gwidth = 11,
  gheight = 13)

