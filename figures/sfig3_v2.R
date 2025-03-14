# %% note that this code is written to be run from the project base directory
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
    tmp <- readRDS("outs/rds/seuMergedPostHSP_forFigures_2025-01-12_04-07-24.rds")
    
    return(tmp)
  })

sfigClusterOrder <- read.csv("figures/sfig4_ordering/cluster_order.csv")

sfigClusterOrder <- factor(sfigClusterOrder$clusterOrder,
  levels = sfigClusterOrder$clusterOrder,
  labels = stringr::str_trim(sfigClusterOrder$clusterOrder))

seu$Disease_Status <- factor(seu$Disease_Status, levels = c("ND", "AAb+", "T1D"))


# one of the clusters has an extra space so just removing it here
seu$manualAnnot <- stringr::str_trim(seu$manualAnnot)

# fix cluster naming
seu$manualAnnot <- gsub("CD8 Tem/Temra", "CD8 Tcm/Tem/Temra", seu$manualAnnot)


################################################################################
# specific annot genes
################################################################################
annotGenes <- read.csv("figures/sfig4_ordering/annot_genes_suppfig.csv")


################################################################################
# %% S4A/B - dot plots!
################################################################################
labellerAdt <- tsaCatalog$cleanName
names(labellerAdt) <- paste0("adt_", tsaCatalog$DNA_ID)

# base adt dot style
figABase <- DotPlot(seu,
  features = paste0("adt_", annotGenes[annotGenes$modality == "ADT", "DNA_ID"]),
  group.by = "manualAnnot",
  cols = "RdYlBu")

figA <- figABase +
  scale_y_discrete(limits = rev(levels(sfigClusterOrder))) +
  textSizeOnlyTheme +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_line(color = "#eeeeee"),
    axis.text.x = element_text(angle = 45, size = BASEPTFONTSIZE, hjust = 1)) +
  scale_x_discrete(labels = labellerAdt)


# base rna dot style
figBBase <- DotPlot(seu,
  features = toupper(annotGenes[annotGenes$modality == "RNA", "geneName"]),
  group.by = "manualAnnot",
  cols = "RdYlBu")

figB <- figBBase +
  scale_y_discrete(limits = rev(levels(sfigClusterOrder))) +
  textSizeOnlyTheme +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_line(color = "#eeeeee"),
    axis.text.x = element_text(angle = 45, size = BASEPTFONTSIZE, hjust = 1))


################################################################################
# %% Final layout and plot all
################################################################################
layout <- c(
  patchwork::area(1, 1, 4, 12), # a
  patchwork::area(5, 1, 8, 12) # b
)

p <- wrap_elements(plot = figA) +
  wrap_elements(plot = figB) +
  plot_annotation(tag_levels = list(letters[1:2])) +
  plot_layout(design = layout) &
  plotTagTheme

pdf("/srv/http/betts/hpap/final_figures/amsesk/pdf/sfig3_v3_final.pdf", width=11, height=13, family="sans")
print(p)
dev.off()

# %%
# saveFinalFigure(
#   plot = p,
#   prefixDir = "figures/outs",
#   fn = "sfig3_v3_final",
#   devices = c("pdf", "png"),
#   addTimestamp = TRUE,
#   gwidth = 11,
#   gheight = 13)
#
