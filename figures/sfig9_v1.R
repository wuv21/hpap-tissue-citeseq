# note that this code is written to be run from the project base directory
# renv::load("/data/hpap-citeseq/hpap-citeseq-analysis")

source("figures/genericFigureSettings.R")
source("scripts/dimPlots.R")
source("scripts/deg.R")
source("scripts/generic.R")
library(ComplexHeatmap)
library(dplyr)
library(Seurat)
library(cowplot)
library(patchwork)
library(stringr)

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
# S9A ridgeline of cd27 for b cells
# include naive t cell as pos control
################################################################################
labellerAdt <- tsaCatalog$cleanName
names(labellerAdt) <- paste0("adt_", tsaCatalog$DNA_ID)

specClusters <- levels(manualClusterOrder)[grepl("^(B|CD4 na|CD8 na)", levels(manualClusterOrder))]
seu_spec_pln <- subset(seu, subset = TissueCondensed == "pLN" & manualAnnot %in% specClusters)

featuresOfInterest <- c(
  "A0154"
)

tmp <- RidgePlot(seu_spec_pln,
  features = featuresOfInterest, group.by = "manualAnnot", combine = FALSE)

tmp <- lapply(tmp, function(x) {
  x$labels$title <- labellerAdt[x$labels$title]
  
  x <- x +
    stat_summary(aes(group = ident, color = ident),
      fun = median, fun.min = "median", fun.max= "median", geom = "crossbar", linetype = "dotted", show.legend = FALSE) +
    subplotTheme +
    theme(
      plot.margin = margin(3,3,3,3),
      axis.text = element_text(size = BASEPTFONTSIZE),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "blank") 
  
  x$layers[[1]]$geom$default_aes$alpha <- 0.3
  
  return(x)
})


################################################################################
# Final layout and plot all
################################################################################
layout <- c(
  area(1, 1, 6, 12)
)

p <- wrap_elements(plot = wrap_plots(tmp, nrow = 1)) +
  plot_annotation(tag_levels = list(LETTERS[1:2])) +
  plot_layout(design = layout) &
  plotTagTheme

saveFinalFigure(
  plot = p,
  prefixDir = "figures/outs",
  fn = "sfig9_final",
  devices = c("pdf", "png"),
  addTimestamp = TRUE,
  gwidth = 6,
  gheight = 5)

