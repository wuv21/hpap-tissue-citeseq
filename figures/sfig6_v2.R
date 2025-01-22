# note that this code is written to be run from the project base directory
renv::load("/data/hpap-citeseq/hpap-citeseq-analysis")

source("figures/genericFigureSettings.R")
source("scripts/dimPlots.R")
source("scripts/deg.R")
source("scripts/generic.R")
library(ComplexHeatmap)
library(dplyr)
library(Seurat)
library(cowplot)
library(presto)
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
    tmp <- readRDS("outs/rds/seuMergedPostHSP_forFigures_2025-01-12_04-07-24.rds")
    
    return(tmp)
  })


customSortAnnotation <- function(x) {
  priority <- c("B", "CD4", "CD8", "NK", "DC", "Monocyte")
  
  x <- sort(unique(x))
  newX <- c()
  
  for (p in priority) {
    pIndices <- grepl(p, x)
    newX <- append(newX, x[pIndices])
    x <- x[!pIndices]
  }
  
  if (length(newX) > 0) {
    newX <- append(newX, x)
  }
  
  
  return(newX)
}


# one of the clusters has an extra space so just removing it here
seu$manualAnnot <- stringr::str_trim(seu$manualAnnot)

# fix cluster naming
seu$manualAnnot <- gsub("CD8 Tem/Temra", "CD8 Tcm/Tem/Temra", seu$manualAnnot)

# cluster order
manualClusterOrder <- unique(seu$manualAnnot)
manualClusterOrder <- factor(manualClusterOrder,
  levels = customSortAnnotation(manualClusterOrder),
  labels = stringr::str_trim(customSortAnnotation(manualClusterOrder)))

# disease status order
seu$Disease_Status <- factor(seu$Disease_Status, levels = c("ND", "AAb+", "T1D"))


################################################################################
# S7A - naive t cell adt.
################################################################################
labellerAdt <- tsaCatalog$cleanName
names(labellerAdt) <- paste0("adt_", tsaCatalog$DNA_ID)

clustersForFigure <- c(
  "CD4 naive #1",
  "CD4 naive #2",
  "CD8 naive #1",
  "CD8 naive #2",
  "CD4 Tcm/Tem",
  "CD8 Tcm/Tem/Temra"
)

seu_pln <- subset(seu, subset = TissueCondensed == "pLN" & manualAnnot %in% clustersForFigure)

featuresOfInterest <- c(
  "A0154",
  "A0390",
  "A0156",
  "A0386",
  "A0063",
  "A0087"
)

tmp <- RidgePlot(
  seu_pln,
  features = featuresOfInterest,
  group.by = "manualAnnot",
  combine = FALSE)

tmp <- lapply(tmp, function(x) {
  x$labels$title <- labellerAdt[x$labels$title]
  
  # plotting up to the 95th percentile for ease of visualization only because
  # of outliers.
  newXMax <- quantile(x$data[, 1], probs = 0.995)
  newXMax <- ceiling(newXMax)
  
  x <- x +
    stat_summary(aes(group = ident, color = ident),
      fun = median, fun.min = "median", fun.max= "median", geom = "crossbar", linetype = "dotted", show.legend = FALSE) +
    scale_y_discrete(limits = rev(clustersForFigure)) +
    scale_x_continuous(limits = c(-0.25, as.numeric(newXMax)), expand = expansion(mult = c(0.05, 0.3))) +
    subplotTheme +
    theme(
      plot.margin = margin(3,3,3,3),
      plot.title = element_text(size = BASEPTFONTSIZE + 1),
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
p <- wrap_elements(plot = wrap_plots(tmp, nrow = 3))

saveFinalFigure(
  plot = p,
  prefixDir = "figures/outs",
  fn = "sfig6_v3_final",
  devices = c("pdf", "png"),
  addTimestamp = TRUE,
  gwidth = 8.5,
  gheight = 7)

