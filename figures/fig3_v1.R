# note that this code is written to be run from the project base directory

source("figures/genericFigureSettings.R")
source("scripts/dimPlots.R")
library(Seurat)
library(cowplot)

set.seed(42)

seu <- tryCatch(
  {
    print(head(seu$TissueCondensed))
    message("Seurat object already exists")
    return(seu)
  },
  error = function(cond) {
    message("Seurat object doesn't exist. Loading now.")
    tmp <- readRDS("outs/rds/seuMergedPostHSP_forFigures_2023-09-12_16-27-41.rds")
    
    return(tmp)
})


manualClusterOrder <- unique(seu$manualAnnot)
manualClusterOrder <- factor(manualClusterOrder,
  levels = customSortAnnotation(manualClusterOrder),
  labels = stringr::str_trim(customSortAnnotation(manualClusterOrder)))

################################################################################
# Figures A-C
################################################################################
figA <- DimPlot(seu, reduction = "rna.umap", group.by = "TissueCondensed") +
  scale_color_manual(values = COLORS[["tissue"]]) +
  labs(title = "Tissue")
figB <- DimPlot(seu, reduction = "rna.umap", group.by = "Disease_Status", shuffle = TRUE) +
  scale_color_manual(values = COLORS[["disease"]], limits = c("ND", "AAb+", "T1D")) +
  labs(title = "Disease Status")
figC <- DimPlot(seu, reduction = "rna.umap", group.by = "manualAnnot") +
  scale_color_discrete(limits = levels(manualClusterOrder)) +
  labs(title = "Cell Type")

legendTheme <- theme(
  legend.margin = margin(0,0,0,0),
  legend.box.margin = margin(-10,-10,-10,-10),
  legend.justification = "left",
  legend.key.size = unit(0.2, 'lines'),
  legend.spacing.x = unit(2, "points"),
  legend.text = element_text(margin = margin(r = 20, unit = "pt"), size = 5),
  legend.title = element_text(size = 7, color = "#000000", hjust = 0))

a_legend <- ggpubr::as_ggplot(get_legend(figA + 
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1), nrow = 1,
      title = "Legend for (A)",
      title.position = "top",
      title.hjust = 0)) +
    legendTheme +
    theme(legend.box.margin = margin(-10,0,-10,-10))))

b_legend <- ggpubr::as_ggplot(get_legend(figB + 
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1), nrow = 1,
      title = "Legend for (B)",
      title.position = "top",
      title.hjust = 0)) +
    legendTheme +
    theme(legend.box.margin = margin(-10,-10,-10,0))))

c_legend <- ggpubr::as_ggplot(get_legend(figC + 
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1), nrow = 15,
      title = "Legend for (C)",
      title.position = "top",
      title.hjust = 0))  +
    legendTheme))

figA <- figA +
  subplotTheme +
  umapPlotThemeNoLeg 

figB <- figB + 
  subplotTheme +
  umapPlotThemeNoLeg + 
  theme(axis.title.y = element_text(color = "#FFFFFF00"),
    plot.margin = unit(c(0, 0, 0, 0), "pt"))

figC <- figC + 
  subplotTheme +
  umapPlotThemeNoLeg +
  theme(axis.title.y = element_text(color = "#FFFFFF00"),
    plot.margin = unit(c(0, 0, 0, 0), "pt"))

figABC <- (figA + figB + figC & labs(x = "UMAP 1", y = "UMAP 2"))

figABCLegend <- (a_legend + b_legend + plot_layout(widths = c(1, 1))) / c_legend + plot_layout(heights = c(1, 15))

figABCLegend[[1]] <- figABCLegend[[1]] +
  theme(plot.margin = margin(t = 0, b = 0, l = 10))

figABCLegend[[2]] <- figABCLegend[[2]] +
  theme(plot.margin = margin(t = 40, b = 30, l = 10))


################################################################################
# Figure D
################################################################################
frequencyDf <- data.frame(
  cluster = seu$manualAnnot,
  tissue = seu$TissueCondensed,
  donor = seu$DonorID,
  diseaseStatus = seu$Disease_Status) %>%
  mutate(diseaseStatus = factor(diseaseStatus, levels = c("ND", "AAb+", "T1D"))) %>%
  mutate(cluster = factor(stringr::str_trim(cluster), levels = levels(manualClusterOrder)))

frequencyPlotTheme <- list(
  theme_bw(),
  coord_cartesian(clip = "off"),
  theme(
    axis.text.x = element_text(color = "#000000", angle = 45, hjust = 1, vjust = 1, size = 4, margin = margin(b = -5, unit = "lines")),
    axis.title.y = element_text(color = "#000000", size = 6),
    axis.text.y = element_text(color = "#000000", size = 4),
    axis.title.x = element_blank(),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_rect(fill = NULL, colour = "#000000", size = 0.25)
    ),
  guides(fill = "none",
    shape = guide_legend(
      override.aes = list(size = 2.5, alpha = 1), nrow = 1,
      title = "Tissue",
      title.position = "top",
      title.hjust = 0.5))
)

# number of cells in each cluster
figD <- frequencyDf %>%
  group_by(cluster, tissue, diseaseStatus) %>%
  summarize(nCells = n()) %>%
  ggplot(aes(x = cluster, y = nCells, shape = tissue, fill = diseaseStatus, group = tissue)) +
  geom_point(color = "#000000", alpha = 0.8, size = 0.75, position = position_dodge(width = 0.75)) +
  # facet_wrap(~ tissue, ncol = 1, strip.position = "right") +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = COLORS[["disease"]]) +
  scale_shape_manual(values = c(21:23)) +
  labs(x = "Total number of cells", y = "Cluster") +
  subplotTheme +
  frequencyPlotTheme +
  legendTheme +
  theme(
    legend.justification = "center", 
    legend.box.margin = margin(-20,0,-20,0),
    legend.title = element_blank())
    # legend.title = element_text(size = 6))
  # coord_flip()

################################################################################
# Fig E heatmap
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
figE <- Heatmap(
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

################################################################################
# Fig F heatmap
################################################################################
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

figF <- Heatmap(
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
  area(1, 1, 2, 7), #abc
  area(1, 8, 2, 12), #legend
  area(3, 1, 4, 12), #d
  area(5, 1, 11, 6), #e
  area(5, 7, 11, 12) #f
)

figABC_final <- figABC +
  plot_annotation(tag_levels = list(LETTERS[1:3])) &
  plotTagTheme

p <- wrap_elements(full = figABC_final, ignore_tag = TRUE) +
  wrap_elements(full = figABCLegend, ignore_tag = TRUE, clip = FALSE) +
  wrap_elements(full = figD + theme(plot.margin = margin(l = 10, t = 10, b = 30)), clip = FALSE) +
  wrap_elements(full = 
    grid.grabExpr(
      draw(figE,
        heatmap_legend_side = "bottom",
        annotation_legend_side = "bottom",
        background = "transparent",
        padding = unit(c(0, 1.5, 0.5, -1), "lines"),
        merge_legend = TRUE)), clip = FALSE) +
  wrap_elements(full = 
      grid.grabExpr(
        draw(figF,
          heatmap_legend_side = "bottom",
          annotation_legend_side = "bottom",
          background = "transparent",
          padding = unit(c(0, 1.5, 0.5, -1), "lines"),
          merge_legend = TRUE)), clip = FALSE) +
  plot_annotation(tag_levels = list(c("D", "E", "F"))) +
  plot_layout(design = layout) &
  plotTagTheme

saveFinalFigure(plot = p, prefixDir = "outs", fn = "fig3_citeseq_intro", devices = c("png", "pdf"), gwidth = 8.5, gheight = 11)
