# note that this code is written to be run from the project base directory
# renv::load("/data/hpap-citeseq/hpap-citeseq-analysis")

source("figures/genericFigureSettings.R")
source("scripts/deg.R")
library(ComplexHeatmap)
library(Seurat)
library(cowplot)
library(presto)
library(ggtext)

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


seuSpleen <- subset(seu, subset = TissueCondensed == "Spleen")
commonGenesFn <- "outs/csv/pLN_findAllMarkers_byManualAnnot_inMoreThan12Clusters.csv"
commonGenes <- read.csv(commonGenesFn)
commonGenesByDisease <- split(commonGenes, commonGenes$cluster)

dcDeg <- findMarkersCombinatorial(
  subset(seuSpleen, subset = manualAnnot == "DC"),
  combVar = "Disease_Status"
) %>%
  select(-p_val_adj)


dcPlot <- plotCombinatorialDEGLollipop(dcDeg %>%
  dplyr::filter(!(gene %in% commonGenesByDisease$`ND`$gene) &
    !(gene %in% commonGenesByDisease$`AAb+`$gene) &
    !(gene %in% commonGenesByDisease$`T1D`$gene)), title = "test") +
  coord_cartesian(clip = "off") +
  theme_classic() +
  subplotTheme +
  theme(
    axis.line = element_line(color = "#000000"),
    panel.grid.major.y = element_line(color = "#dddddd"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(color = "#000000", size = 6),
    axis.text.y = element_text(color = "#000000", size = 5),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.position = "bottom",
    legend.box.margin = margin(t = -10, b = 0),
    plot.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_markdown(hjust = 0.5, size = 6)) +
  facet_wrap(~ matchup,
    scales = "free",
    labeller = as_labeller(c(
      "ND_vs_AAb+" = glue("<span style='color: {COLORS$disease['AAb+']};'>up in AAb+</span> | <span style='color: {COLORS$disease['ND']};'>up in ND</span>"),
      "T1D_vs_ND" = glue("<span style='color: {COLORS$disease['ND']};'>up in ND</span> | <span style='color: {COLORS$disease['T1D']};'>up in T1D</span>"),
      "T1D_vs_AAb+" = glue("<span style='color: {COLORS$disease['AAb+']};'>up in AAb+</span> | <span style='color: {COLORS$disease['T1D']};'>up in T1D</span>"))))


monoOneDeg <- findMarkersCombinatorial(
  subset(seuSpleen, subset = manualAnnot == "Monocyte #1"),
  combVar = "Disease_Status"
) %>%
  select(-p_val_adj)


monoPlot <- plotCombinatorialDEGLollipop(monoOneDeg %>%
    dplyr::filter(!(gene %in% commonGenesByDisease$`ND`$gene) &
        !(gene %in% commonGenesByDisease$`AAb+`$gene) &
        !(gene %in% commonGenesByDisease$`T1D`$gene)), title = "test") +
  coord_cartesian(clip = "off") +
  theme_classic() +
  subplotTheme +
  theme(
    axis.line = element_line(color = "#000000"),
    panel.grid.major.y = element_line(color = "#dddddd"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(color = "#000000", size = 6),
    axis.text.y = element_text(color = "#000000", size = 5),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.position = "bottom",
    legend.box.margin = margin(t = -10, b = 0),
    plot.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_markdown(hjust = 0.5, size = 6)) +
  facet_wrap(~ matchup,
    scales = "free",
    labeller = as_labeller(c(
      "ND_vs_AAb+" = glue("<span style='color: {COLORS$disease['AAb+']};'>up in AAb+</span> | <span style='color: {COLORS$disease['ND']};'>up in ND</span>"),
      "T1D_vs_ND" = glue("<span style='color: {COLORS$disease['ND']};'>up in ND</span> | <span style='color: {COLORS$disease['T1D']};'>up in T1D</span>"),
      "T1D_vs_AAb+" = glue("<span style='color: {COLORS$disease['AAb+']};'>up in AAb+</span> | <span style='color: {COLORS$disease['T1D']};'>up in T1D</span>"))))


################################################################################
# Final layout and plot all
################################################################################
# layout <- c(
#   patchwork::area(1, 1, 4, 6), # a flow 1
#   patchwork::area(5, 1, 6, 6), # c heatmap
# )

p <- dcPlot /
  monoPlot +
  plot_annotation(tag_levels = list(LETTERS[1:2])) &
  plotTagTheme

saveFinalFigure(
  plot = p,
  prefixDir = "figures/outs",
  fn = "sfigZ_v3_final",
  devices = c("pdf", "png"),
  addTimestamp = TRUE,
  gwidth = 7,
  gheight = 7)

