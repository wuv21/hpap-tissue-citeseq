# renv::load("/data/hpap-citeseq/hpap-citeseq-analysis")

source("figures/genericFigureSettings.R")
source("scripts/dimPlots.R")
library(Seurat)
library(ggpubr)
library(ggrastr)
library(ComplexHeatmap)
library(circlize)
library(forcats)
library(cluster)
library(ggrepel)

set.seed(42)

parentDir <- "figures/greg_flow_data"

################################################################################
# A and B - flow umaps from omiq
# code is modified from @Greg
################################################################################
omiq_umap_coordinate_lnmc <- read.csv(paste0(parentDir, "/OMIQ_HPAP099_Panc-Tail_CD45+.csv"))
omiq_umap_coordinate_snmc <- read.csv(paste0(parentDir, "/OMIQ_HPAP099_Spleen_CD45+.csv"))

legendTheme <- theme(
  legend.margin = margin(0,0,0,0),
  legend.box.margin = margin(-10,-10,-10,-10),
  legend.justification = "left",
  legend.key.size = unit(0.1, 'lines'),
  legend.spacing.x = unit(3, "points"),
  legend.text = element_text(margin = margin(r = 12, unit = "pt"), size = 4.5),
  legend.title = element_text(size = 8, color = "#000000", hjust = 0))

figA <- omiq_umap_coordinate_lnmc %>%
  ggplot(aes(x = umap_1, y = umap_2, color = OmiqFilter)) +
  geom_point(size = 0.2) +
  labs(
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_classic() +
  scale_color_manual(values = COLORS$expanded_palette)

figB <- omiq_umap_coordinate_snmc %>%
  ggplot(aes(x = umap_1, y = umap_2, color = OmiqFilter)) +
  geom_point(size = 0.2) +
  labs(
    x = "UMAP 1",
    y = "UMAP 2",
  ) +
  theme_classic() +
  scale_color_manual(values = COLORS$expanded_palette)

ab_legend <- ggpubr::as_ggplot(get_legend(figA + 
    guides(colour = guide_legend(
      override.aes = list(size = 2.5, alpha = 1),
      nrow = 10,
      title = "Legend for (A) and (B)",
      title.position = "top",
      title.hjust = 0)) +
    legendTheme +
    theme(legend.box.margin = margin(0,0,0,0))))

figA <- figA +
  subplotTheme +
  umapPlotThemeNoLeg  +
  theme(plot.margin = margin(0,0,0,0))

figB <- figB + 
  subplotTheme +
  umapPlotThemeNoLeg + 
  theme(
    plot.margin = margin(0,0,0,0),
    axis.title.y = element_text(color = "#FFFFFF00"))


################################################################################
# C - lineage data heatmap
# code is modified from @Greg
################################################################################
dfLineageFilter <- readRDS(paste0(parentDir, "/rds/dfLineageFilter.rds"))
dfHeatmapNorm <- dfLineageFilter %>%
  group_by(metric) %>%
  mutate(z_score = scale(value)) %>%
  dplyr::select("HPAP Donor":metric, z_score)

LineageList <- c("Eosinophils", "Neutrophils", "B cells","ILCs",
  "Imm. Granulocytes", "HSCs", "DCs", "Monocytes", "NK cells",
  "T cells", "CD8+ T cells", "CD8 Tcm", "CD8 Tem", "CD8 Temra",
  "CD8 Tn", "CD8 Tnl", "CD4+ T cells", "CD4 Tcm", "CD4 Tem", "CD4 Temra",
  "CD4 Tn", "CD4 Tnl", "CD4 Mem Tregs", "CD4 Mem gcTfh", "CD4 Mem mTfh")
LinTissueOrder <- c("Spleen", "Head pLN", "Body pLN", "Tail pLN", "SMA mLN", "MES mLN")
SampleVariablesLin <- c("HPAP Donor", "Tissue")

dfTotalLineageHeatmap <- dfHeatmapNorm %>%
  pivot_wider(names_from = metric, values_from = z_score) %>%
  dplyr::select(all_of(c(SampleVariablesLin, LineageList))) %>%
  mutate(Sample_ID = paste("S", as.character(row_number()), sep = "")) %>%
  relocate(Sample_ID, .before = `HPAP Donor`) %>%
  mutate(Tissue = str_replace_all(Tissue,"_", " "))

dfTotalLineageHeatmap <- as.data.frame(dfTotalLineageHeatmap)
rownames(dfTotalLineageHeatmap) <- dfTotalLineageHeatmap$Sample_ID
dfTotalLineageHeatmap <- dplyr::select(dfTotalLineageHeatmap, `HPAP Donor`:ncol(dfTotalLineageHeatmap))

dfLinHeatmapMetadata <- dplyr::select(dfTotalLineageHeatmap, Tissue)
dfTotalLineageHeatmap <- dplyr::select(dfTotalLineageHeatmap, "Eosinophils":ncol(dfTotalLineageHeatmap))
dfTotalLineageHeatmap <- t(dfTotalLineageHeatmap)

col_fun_LineageHM <- circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")) #set maxima at 5th and 95th percentiles
LinAnn <- HeatmapAnnotation(
  df = dfLinHeatmapMetadata,
  which = "col",
  col = list(`Tissue` = COLORS$tissue_detailed),
  simple_anno_size = unit(5, "points"),
  height = unit(5, "points"),
  annotation_name_gp = gpar(fontsize = 8),
  annotation_legend_param = list(
    title = "Tissue Legend for (C) and (D)",
    at = LinTissueOrder,
    labels = LinTissueOrder,
    nrow = 1,
    legend_height = unit(1, "mm"),
    legend_width = unit(0.5, "mm"),
    labels_gp = gpar(fontsize = 5),
    grid_width = unit(1, "mm"),
    grid_height = unit(1, "mm"),
    title_gp = gpar(fontsize = 8, fontface = "plain")
  ))

figC <- Heatmap(
  dfTotalLineageHeatmap,
  name = "Z-score",
  col = col_fun_LineageHM,
  show_column_names = FALSE,
  top_annotation = LinAnn,
  cluster_columns = diana,
  row_km = 3,
  row_title_gp = gpar(fontsize = 8),
  row_title_rot = 0,
  row_names_gp = gpar(fontsize = 6),
  row_dend_width = unit(3, "points"),
  column_dend_height = unit(12, "points"),
  heatmap_legend_param = list(
    direction = "horizontal",
    title_position = "topleft",
    labels_gp = gpar(fontsize = 5),
    title_gp = gpar(fontsize = 8, fontface = "plain"),
    grid_height = unit(1.5, "mm"),
    legend_height = unit(2, "mm"),
    legend_width = unit(10, "mm")),
  )

################################################################################
# D - Tissue Distribution -
# Generate PCA-compatible dataframe and metadata for entire dataset
# code is modified from @Greg
################################################################################

dfPCA <- dfLineageFilter %>%
  dplyr::select(1:Tissue, metric, value) %>%
  pivot_wider(id_cols = `HPAP Donor`:`Tissue`, names_from = metric, values_from = value) %>%
  mutate(Sample_ID = paste("S", as.character(row_number()), sep = "")) %>%
  relocate(Sample_ID, .before = `HPAP Donor`) %>%
  as.data.frame()
rownames(dfPCA) <- dfPCA$Sample_ID

dfLineagePCA <- dplyr::select(dfPCA, LineageList)
TotalLinPCA <- prcomp(dfLineagePCA, scale. = TRUE)

dfTotalLinPCA <- TotalLinPCA$x %>%
  as.data.frame() %>%
  mutate(Sample_ID = rownames(dfLineagePCA)) %>%
  full_join(dfPCA, by = "Sample_ID")

Lin_var_explained <- TotalLinPCA$sdev^2/sum(TotalLinPCA$sdev^2)

figD <- dfTotalLinPCA %>%
  mutate(Tissue = str_replace(Tissue, "_", " ")) %>%
  mutate(Tissue = fct_relevel(Tissue, "Spleen", "SMA mLN", "MES mLN", "Head pLN", "Body pLN", "Tail pLN")) %>%
  ggplot(aes(PC1, PC2, color = Tissue)) +
  geom_point(size = 0.5) +
  labs(
    x = paste0("PC1: ", round(Lin_var_explained[1]*100, 1), "%"),
    y = paste0("PC2: ", round(Lin_var_explained[2]*100, 1), "%")) +
  scale_color_manual(values = COLORS$tissue_detailed) +
  theme_classic() +
  theme(aspect.ratio = 1) +
  subplotTheme +
  umapPlotThemeNoLeg + 
  theme(
    plot.margin = margin(0,0,0,0))

dfTotalLinBiplot <- TotalLinPCA$rotation %>%
  as.data.frame() %>%
  mutate(Population = rownames(TotalLinPCA$rotation))

figE <- dfTotalLinBiplot %>%
  ggplot(aes(PC1, PC2)) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
    arrow = arrow(angle = 20, 
      length = unit(0.02, "inches"),
      ends = "last",
      type = "closed"),
    color = "grey",
    size = 0.2) +
  labs(
    x = paste0("PC1: ", round(Lin_var_explained[1]*100, 1), "%"),
    y = paste0("PC2: ", round(Lin_var_explained[2]*100, 1), "%")) +
  geom_text_repel(label = dfTotalLinBiplot$Population, 
    size = 1, 
    min.segment.length = 0.2,
    segment.color = "darkblue",
    segment.alpha = 0.5,
    segment.size = 0.2,
    seed = 42) +
  theme_classic() +
  theme(aspect.ratio = 1) +
  subplotTheme +
  umapPlotThemeNoLeg + 
  theme(
    plot.margin = margin(0,0,0,0),
    axis.title.y = element_text(color = "#FFFFFF00"))

################################################################################
# Figures F-H
################################################################################
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


# one of the clusters has an extra space so just removing it here
seu$manualAnnot <- stringr::str_trim(seu$manualAnnot)

# fix cluster naming
seu$manualAnnot <- gsub("CD8 Tem/Temra", "CD8 Tcm/Tem/Temra", seu$manualAnnot)

# cluster order
manualClusterOrder <- unique(seu$manualAnnot)
sortedClust <- customSortAnnotation(manualClusterOrder)
manualClusterOrder <- factor(manualClusterOrder,
  levels = sortedClust,
  labels = sortedClust)


figF <- DimPlot(seu, reduction = "rna.umap", group.by = "TissueCondensed") +
  scale_color_manual(values = COLORS[["tissue"]]) +
  labs(title = "Tissue")
figG <- DimPlot(seu, reduction = "rna.umap", group.by = "Disease_Status", shuffle = TRUE) +
  scale_color_manual(values = COLORS[["disease"]], limits = c("ND", "AAb+", "T1D")) +
  labs(title = "Disease Status")
figH <- DimPlot(seu, reduction = "rna.umap", group.by = "manualAnnot") +
  scale_color_discrete(limits = levels(manualClusterOrder)) +
  labs(title = "Cell Type")

f_legend <- ggpubr::as_ggplot(get_legend(figF + 
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1), nrow = 1,
      title = "Legend for (F)",
      title.position = "top",
      title.hjust = 0)) +
    legendTheme +
    theme(legend.box.margin = margin(0,0,-10,-10))))

g_legend <- ggpubr::as_ggplot(get_legend(figG + 
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1), nrow = 1,
      title = "Legend for (G)",
      title.position = "top",
      title.hjust = 0)) +
    legendTheme +
    theme(legend.box.margin = margin(0,0,-10,0))))

h_legend <- ggpubr::as_ggplot(get_legend(figH + 
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1), nrow = 15,
      title = "Legend for (H)",
      title.position = "top",
      title.hjust = 0))  +
    legendTheme))

figF <- figF +
  subplotTheme +
  umapPlotThemeNoLeg

figG <- figG + 
  subplotTheme +
  umapPlotThemeNoLeg + 
  theme(axis.title.y = element_text(color = "#FFFFFF00"),
    plot.margin = unit(c(0, 0, 0, 0), "pt"))

figH <- figH + 
  subplotTheme +
  umapPlotThemeNoLeg +
  theme(axis.title.y = element_text(color = "#FFFFFF00"),
    plot.margin = unit(c(0, 0, 0, 0), "pt"))

figFGH <- (figF + figG + figH & labs(x = "UMAP 1", y = "UMAP 2"))

figFGHLegend <- (f_legend + g_legend + plot_layout(widths = c(1, 1))) / h_legend + plot_layout(heights = c(1, 15))

figFGHLegend[[1]] <- figFGHLegend[[1]] +
  theme(plot.margin = margin(t = -30, b = 0, l = 10))

figFGHLegend[[2]] <- figFGHLegend[[2]] +
  theme(plot.margin = margin(t = 30, b = -60, l = 10))


################################################################################
# Figure I
################################################################################
frequencyDf <- data.frame(
  cluster = seu$manualAnnot,
  tissue = seu$TissueCondensed,
  donor = seu$DonorID,
  diseaseStatus = seu$Disease_Status) %>%
  mutate(diseaseStatus = factor(diseaseStatus, levels = c("ND", "AAb+", "T1D"))) %>%
  mutate(cluster = factor(stringr::str_trim(cluster), levels = levels(manualClusterOrder))) %>%
  mutate(tissue = factor(tissue, levels = c("pLN", "mesLN", "Spleen"), labels = c("pLN", "mLN", "Spleen")))

frequencyPlotTheme <- list(
  theme_bw(),
  coord_cartesian(clip = "off"),
  theme(
    axis.text.x = element_text(color = "#000000", angle = 45, hjust = 1, vjust = 1, size = 4, margin = margin(b = -5, unit = "lines")),
    axis.title.y = element_text(color = "#000000", size = 8),
    axis.text.y = element_text(color = "#000000", size = 6),
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
figI <- frequencyDf %>%
  group_by(cluster, tissue, diseaseStatus) %>%
  summarize(nCells = n()) %>%
  ggplot(aes(x = cluster, y = nCells, shape = tissue, fill = diseaseStatus, group = tissue)) +
  geom_point(color = "#000000", alpha = 0.95, size = 0.75, position = position_dodge(width = 0.75)) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = COLORS[["disease"]]) +
  scale_shape_manual(values = c(21:23)) +
  labs(y = "Cell count")+
  subplotTheme +
  frequencyPlotTheme +
  legendTheme +
  theme(
    plot.margin = margin(b = 15),
    plot.background = element_rect(fill = "#ffffff00", color = "#ffffff00"),
    legend.justification = "center", 
    legend.box.margin = margin(-20,0,-20,0),
    legend.title = element_blank())

################################################################################
# Final layout and plot all
################################################################################
layout <- c(
  area(t = 1, l = 1, b = 2, r = 8), #AB
  area(t = 1, l = 9, b = 2, r = 11), #ab legend
  area(t = 1, l = 12, b = 4, r = 20), #C
  area(t = 3, l = 1, b = 4, r = 11), #DE
  area(t = 5, l = 1, b = 6, r = 12), #fgh
  area(t = 5, l = 13, b = 6, r = 20), #fgh legend
  area(t = 7, l = 1, b = 8, r = 20) #i
)

figAB_final <- (figA + figB) +
  plot_annotation(tag_levels = list(LETTERS[1:2])) &
  plotTagTheme

figDE_final <- (figD + figE) +
  plot_annotation(tag_levels = list(LETTERS[4:5])) &
  plotTagTheme

figFGH_final <- figFGH +
  plot_annotation(tag_levels = list(LETTERS[6:8])) &
  plotTagTheme

p <- wrap_elements(full = figAB_final, ignore_tag = TRUE) +
  wrap_elements(full = ab_legend, clip = FALSE, ignore_tag = TRUE) +
  wrap_elements(full = grid.grabExpr(
    draw(figC,
      heatmap_legend_side = "bottom",
      annotation_legend_side = "bottom",
      background = "transparent",
      padding = unit(c(1, 3.5, 1.5, 1.5), "lines"),
      merge_legend = TRUE)), clip = FALSE) +
  wrap_elements(full = figDE_final, ignore_tag = TRUE) +
  wrap_elements(full = figFGH_final, ignore_tag = TRUE) +
  wrap_elements(full = figFGHLegend, ignore_tag = TRUE, clip = FALSE) +
  wrap_elements(full = figI + theme(plot.margin = margin(l = 7.5, t = 20, b = 40)), clip = FALSE, ) +
  plot_annotation(tag_levels = list(c("C", "I"))) +
  plot_layout(design = layout) &
  plotTagTheme

saveFinalFigure(
  plot = p,
  prefixDir = "figures/outs",
  fn = "fig1_v2_final",
  devices = c("pdf", "png"),
  addTimestamp = TRUE,
  gwidth = 8.5,
  gheight = 8)

