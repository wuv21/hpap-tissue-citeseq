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
    tmp <- readRDS("outs/rds/seuMergedPostHSP_forFigures_2023-09-17_09-03-10.rds")
    
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


manualClusterOrder <- unique(seu$manualAnnot)
manualClusterOrder <- factor(manualClusterOrder,
  levels = customSortAnnotation(manualClusterOrder),
  labels = stringr::str_trim(customSortAnnotation(manualClusterOrder)))

seu$Disease_Status <- factor(seu$Disease_Status, levels = c("ND", "AAb+", "T1D"))




################################################################################
# S5A - heatmap of deg
################################################################################
clusterDEG <- read.csv("outs/tsv/DE_RNA_cluster.tsv", sep = "\t")
top10DegByCluster <- clusterDEG %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.05) %>%
  slice_max(avg_log2FC, n = 10)

top10DegAvgExpression <- AverageExpression(
  object = seu,
  assays = "RNA",
  return.seurat = FALSE,
  features = top10DegByCluster$gene,
  group.by = "manualAnnot",
  slot = "data")

top10DegAvgExpression <- top10DegAvgExpression$RNA
figA <- Heatmap(
  matrix = t(scale(t(top10DegAvgExpression))),
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_names = FALSE,
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 5),
  column_names_rot = 45,
  row_km = length(colnames(top10DegAvgExpression)),
  right_annotation = rowAnnotation(
    foo = anno_block(gp = gpar(fill = "#000000"), width = unit(0.7, "mm")),
    bar = anno_block(
      graphics = function(index, levels) {
        lbls <- rownames(top10DegAvgExpression)[index]
        lbls <- str_wrap(paste(lbls, collapse = " "), width = 30)
        
        grid.rect(gp = gpar(fill = NA, col = NA))
        txt = paste(lbls, collapse = ",")
        grid.text(txt, 0.01, 0.5, rot = 0, hjust = 0, gp = gpar(fontsize = 4))
      },
      width = unit(3, "cm"))
  ),
  name = "Scaled Average Expression",
  heatmap_legend_param = list(
    direction = "horizontal",
    title_position = "topcenter",
    title_gp = gpar(fontsize = 6, fontface = "plain"),
    labels_gp = gpar(fontsize = 6),
    legend_height = unit(0.08, "cm"),
    legend_width = unit(2, "cm")),
  row_title_gp = gpar(fontsize = 5)
)


clusterDEA <- read.csv("outs/tsv/DE_ADT_cluster.tsv", sep = "\t")
top10DeaByCluster <- clusterDEA %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.05) %>%
  slice_max(avg_log2FC, n = 10)

top10DeaAvgExpression <- AverageExpression(
  object = seu,
  assays = "adt",
  features = top10DeaByCluster$gene,
  group.by = "manualAnnot",
  slot = "data")

tsaCatalog <- readRDS("rds/tsa_catalog.rds")
labellerAdt <- tsaCatalog$cleanName
names(labellerAdt) <- paste0("adt_", tsaCatalog$DNA_ID)

top10DeaAvgExpression <- top10DeaAvgExpression$adt
rownames(top10DeaAvgExpression) <- labellerAdt[paste0("adt_", rownames(top10DeaAvgExpression))]

################################################################################
# S5B - heatmap of dea
################################################################################
figB <- Heatmap(
  matrix = t(scale(t(top10DeaAvgExpression))),
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_names = FALSE,
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 5),
  column_names_rot = 45,
  row_km = length(colnames(top10DeaAvgExpression)),
  right_annotation = rowAnnotation(
    foo = anno_block(gp = gpar(fill = "#000000"), width = unit(0.7, "mm")),
    bar = anno_block(
      graphics = function(index, levels) {
        lbls <- rownames(top10DeaAvgExpression)[index]
        lbls <- str_wrap(paste(lbls, collapse = " | "), width = 30)
        
        grid.rect(gp = gpar(fill = NA, col = NA))
        txt = paste(lbls, collapse = ",")
        grid.text(txt, 0.01, 0.5, rot = 0, hjust = 0, gp = gpar(fontsize = 4))
      },
      width = unit(3, "cm"))
  ),
  name = "Scaled Average Expression",
  heatmap_legend_param = list(
    direction = "horizontal",
    title_position = "topcenter",
    title_gp = gpar(fontsize = 6, fontface = "plain"),
    labels_gp = gpar(fontsize = 6),
    legend_height = unit(0.08, "cm"),
    legend_width = unit(2, "cm")),
  row_title_gp = gpar(fontsize = 5)
)

################################################################################
# Final layout and plot all
################################################################################
layout <- c(
  area(1, 1, 6, 6), #a
  area(1, 7, 6, 12) #b
)

p <- wrap_elements(full = 
    grid.grabExpr(
      draw(figA,
        heatmap_legend_side = "bottom",
        annotation_legend_side = "bottom",
        background = "transparent",
        padding = unit(c(0, 1.5, 0.5, -1), "lines"),
        merge_legend = TRUE)), clip = FALSE) +
  wrap_elements(full = 
      grid.grabExpr(
        draw(figB,
          heatmap_legend_side = "bottom",
          annotation_legend_side = "bottom",
          background = "transparent",
          padding = unit(c(0, 1.5, 0.5, -1), "lines"),
          merge_legend = TRUE)), clip = FALSE) +
  plot_annotation(tag_levels = list(LETTERS[1:2])) +
  plot_layout(design = layout) &
  plotTagTheme

saveFinalFigure(
  plot = p,
  prefixDir = "figures/outs",
  fn = "sfig5_final",
  devices = c("pdf", "png"),
  addTimestamp = TRUE,
  gwidth = 8.5,
  gheight = 9)

