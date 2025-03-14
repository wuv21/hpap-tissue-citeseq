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

################################################################################
# %% A and B - flow umaps from omiq
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
  legend.title = element_text(size = 6, color = "#000000", hjust = 0))

figA <- omiq_umap_coordinate_lnmc %>%
  ggplot(aes(x = umap_1, y = umap_2, color = OmiqFilter)) +
  geom_point(size = 0.1) +
  labs(
    x = "UMAP 1",
    y = "UMAP 2",
    title = "Lymph Node: HPAP-099"
  ) +
  theme_classic() +
  scale_color_manual(values = COLORS$expanded_palette)

figB <- omiq_umap_coordinate_snmc %>%
  ggplot(aes(x = umap_1, y = umap_2, color = OmiqFilter)) +
  geom_point(size = 0.1) +
  labs(
    x = "UMAP 1",
    y = "UMAP 2",
    title = "Spleen: HPAP-099"
  ) +
  theme_classic() +
  scale_color_manual(values = COLORS$expanded_palette)

ab_legend <- ggpubr::as_ggplot(get_legend(figA + 
    guides(colour = guide_legend(
      override.aes = list(size = 2.5, alpha = 1),
      nrow = 10,
      title = "Cell Type",
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
# %% C - lineage data heatmap
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
DisStatusOrder = c("ND", "AAb+", "T1D")
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
tissueLegend = ComplexHeatmap::Legend(
    title = "Tissue",
    at = LinTissueOrder,
    labels = LinTissueOrder,
    nrow = 6,
    legend_gp = gpar(fill= COLORS$tissue_detailed),
    legend_height = unit(6, "mm"),
    legend_width = unit(0.5, "mm"),
    labels_gp = gpar(fontsize = 5),
    grid_width = unit(1, "mm"),
    grid_height = unit(1, "mm"),
    title_gp = gpar(fontsize = 6, fontface = "plain")
)
LinAnn <- HeatmapAnnotation(
  df = dfLinHeatmapMetadata,
  which = "col",
  col = list(`Tissue` = COLORS$tissue_detailed),
  simple_anno_size = unit(5, "points"),
  height = unit(5, "points"),
  annotation_name_gp = gpar(fontsize = 8),
  show_legend = TRUE,
  annotation_legend_param = list(
    title = "Tissue",
    at = LinTissueOrder,
    labels = LinTissueOrder,
    nrow = 1,
    legend_height = unit(1, "mm"),
    legend_width = unit(0.5, "mm"),
    labels_gp = gpar(fontsize = 5),
    grid_width = unit(1, "mm"),
    grid_height = unit(1, "mm"),
    title_gp = gpar(fontsize = 6, fontface = "plain")
  )
)

# figCheatmapLegend = ComplexHeatmap::Legend(
#     title = "Z-score",
#     col_fun = col_fun_LineageHM,
#     direction = "horizontal",
#     title_position = "topcenter",
#     labels_gp = gpar(fontsize = 5),
#     title_gp = gpar(fontsize = 6, fontface = "plain"),
#     grid_height = unit(1.5, "mm"),
#     legend_height = unit(2, "mm"),
#     legend_width = unit(10, "mm")
# )
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
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    direction = "horizontal",
    title_position = "topcenter",
    labels_gp = gpar(fontsize = 5),
    title_gp = gpar(fontsize = 6, fontface = "plain"),
    grid_height = unit(1.5, "mm"),
    legend_height = unit(2, "mm"),
    legend_width = unit(10, "mm")
    )
  )


################################################################################
# %% D - Tissue Distribution -
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
  theme(
    plot.margin = margin(0,0,0,0)
  )

d_legend = ggpubr::as_ggplot(get_legend(figD + 
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1), ncol=1,
      title = "Tissue",
      title.position = "top",
      title.hjust = 0)) +
    legendTheme +
    theme(legend.box.margin = margin(0,0,-10,-10))))

figD = figD + 
  umapPlotThemeNoLeg

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
# %% Figures F-H
################################################################################

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

# rename mesLN as mLN
seu$TissueCondensed <- ifelse(seu$TissueCondensed == "mesLN", "mLN", seu$TissueCondensed)

figF <- DimPlot(seu, reduction = "rna.umap", shuffle = TRUE, group.by = "TissueCondensed") +
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
      # title = "Legend for (F)",
      title = "Tissue (e)",
      title.position = "top",
      title.hjust = 0)) +
    legendTheme +
    theme(legend.box.margin = margin(0,0,-10,-10))))

g_legend <- ggpubr::as_ggplot(get_legend(figG + 
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1), nrow = 1,
      title = "Disease Status (f)",
      title.position = "top",
      title.hjust = 0)) +
    legendTheme +
    theme(legend.box.margin = margin(0,0,-10,0))))

h_legend <- ggpubr::as_ggplot(get_legend(figH + 
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1), nrow = 15,
      title = "Cell Type (g)",
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
# %% Figure I
################################################################################
frequencyDf <- data.frame(
  cluster = seu$manualAnnot,
  tissue = seu$TissueCondensed,
  donor = seu$DonorID,
  diseaseStatus = seu$Disease_Status) %>%
  mutate(diseaseStatus = factor(diseaseStatus, levels = c("ND", "AAb+", "T1D"))) %>%
  mutate(cluster = factor(stringr::str_trim(cluster), levels = levels(manualClusterOrder))) %>%
  mutate(tissue = factor(tissue, levels = c("pLN", "mLN", "Spleen"), labels = c("pLN", "mLN", "Spleen")))

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
  )
  # guides(
  #   fill = guide_legend(
  #     override.aes = list(size = 2.5, alpha = 1, color="#FFFFFF", shape=21), nrow = 1,
  #     title = "Disease Status",
  #     title.position = "top",
  #     title.hjust = 0.5),
  #   shape = guide_legend(
  #     override.aes = list(size = 2.5, alpha = 1), nrow = 1,
  #     title = "Tissue",
  #     title.position = "top",
  #     title.hjust = 0.5)
  # )
)

# number of cells in each cluster
SHAPES = c(21:23)
names(SHAPES)=c("pLN", "mLN", "Spleen")
SHAPES
figI <- frequencyDf %>%
  group_by(cluster, tissue, diseaseStatus) %>%
  summarize(nCells = n()) %>%
  ggplot(aes(x = cluster, y = nCells, shape = tissue, fill = diseaseStatus, group = tissue)) +
  geom_point(color = "#000000", alpha = 0.95, size = 0.75, position = position_dodge(width = 0.75)) +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = COLORS[["disease"]]) +
  # scale_shape_manual(values = SHAPES) +
  scale_shape_manual(values = c(21:23)) +
  labs(y = "Cell count")+
  subplotTheme +
  frequencyPlotTheme +
  legendTheme +
  guides(shape="none", fill="none", color="none")+
  theme(
    plot.margin = margin(b = 15),
    plot.background = element_rect(fill = "#ffffff00", color = "#ffffff00"),
    legend.justification = "center", 
    legend.box.margin = margin(-10,0,-10,0),
    legend.title = element_blank())

################################################################################
# %% Final layout and plot all
################################################################################
layout <- c(
  patchwork::area(t = 1, l = 1, b = 2, r = 8), #AB
  patchwork::area(t = 1, l = 9, b = 2, r = 11), #ab legend
  patchwork::area(t = 1, l = 12, b = 4, r = 20), #C

  # patchwork::area(t = 3, l = 1, b = 4, r = 11), #DE
  patchwork::area(t = 3, l = 1, b = 4, r = 5), #NewC
  patchwork::area(t = 3, l = 6, b = 4, r = 7), #NewC legend
  patchwork::area(t = 3, l = 7, b = 4, r = 11), #NewD
  patchwork::area(t = 5, l = 1, b = 6, r = 12), #fgh
  patchwork::area(t = 5, l = 13, b = 6, r = 20), #fgh legend
  patchwork::area(t = 7, l = 1, b = 8, r = 20) #i
)

figAB_final <- (figA + figB) +
  plot_annotation(tag_levels = list(c("a", ""))) &
  plotTagTheme

figD_final = figD +
  plot_annotation(tag_levels = list(c("c"))) &
  plotTagTheme

figE_final = figE +
  plot_annotation(tag_levels = list(c("d"))) &
  plotTagTheme

# figDE_final <- (figD + figE) +
#   plot_annotation(tag_levels = list(c("c","d"))) &
#   plotTagTheme
# figDE_final

# figDE_legends = grid.grabExpr(draw(ComplexHeatmap::packLegend(tissueLegend, figCheatmapLegend, direction="vertical")))
# figD_legend = grid.grabExpr(draw(ComplexHeatmap::packLegend(tissueLegend, direction="vertical"), just="right"))

figH_leg = ComplexHeatmap::Legend(
    title = NULL,
    at = c("pLN", "mLN", "Spleen", "ND", "AAb+", "T1D"),
    labels = c("pLN", "mLN", "Spleen", "ND", "AAb+", "T1D"),
    nrow = 1,
    legend_gp = gpar(col = COLORS$tissue_detailed),
    legend_height = unit(1, "mm"),
    legend_width = unit(0.5, "mm"),
    labels_gp = gpar(fontsize = 5),
    grid_width = unit(2, "mm"),
    grid_height = unit(1, "mm"),
    title_gp = gpar(fontsize = 6, fontface = "plain"),
    graphics = list(
                    function(x,y,w,h) grid.points(x, y, gp=gpar(fill="#ffffff", col="#000000"), pch=21, size=unit(3,"mm")),
                    function(x,y,w,h) grid.points(x, y, gp=gpar(fill="#ffffff", col="#000000"), pch=22, size=unit(3,"mm")),
                    function(x,y,w,h) grid.points(x, y, gp=gpar(fill="#ffffff", col="#000000"), pch=23, size=unit(3,"mm")),
                    function(x,y,w,h) grid.rect(x, y, w*0.33, h, gp=gpar(fill=COLORS$disease["ND"], col=COLORS$disease["ND"])),
                    function(x,y,w,h) grid.rect(x, y, w*0.33, h, gp=gpar(fill=COLORS$disease["AAb+"], col=COLORS$disease["AAb+"])),
                    function(x,y,w,h) grid.rect(x, y, w*0.33, h, gp=gpar(fill=COLORS$disease["T1D"], col=COLORS$disease["T1D"]))
    )
)

figFGH_final <- figFGH +
  plot_annotation(tag_levels = list(c("e","f","g"))) &
  plotTagTheme

figH_final = figI + 
  grid.grabExpr(draw(figH_leg)) + 
  plot_layout(ncol=1, heights=c(45,1))

p <- wrap_elements(full = figAB_final, ignore_tag = TRUE) +
  wrap_elements(full = ab_legend, clip = FALSE, ignore_tag = TRUE) +
  wrap_elements(full = grid.grabExpr(
    draw(figC,
      heatmap_legend_side = "bottom",
      annotation_legend_side = "bottom",
      background = "transparent",
      padding = unit(c(1, 3.5, 1.5, 1.5), "lines"),
      merge_legend = TRUE)), clip = FALSE) +
  wrap_elements(full = figD_final, ignore_tag = TRUE) +
  wrap_elements(full = d_legend, ignore_tag=TRUE) +
  wrap_elements(full = figE_final, ignore_tag = TRUE) +
  wrap_elements(full = figFGH_final, ignore_tag = TRUE) +
  wrap_elements(full = figFGHLegend, ignore_tag = TRUE, clip = FALSE) +
  wrap_elements(full = figH_final, clip=FALSE)+
  plot_annotation(tag_levels = list(c("b", "h"))) +
  plot_layout(design = layout) &
  plotTagTheme

pdf("/srv/http/betts/hpap/final_figures/amsesk/pdf/fig1_v3_final.pdf", width=8.5, height=8, family="sans")
print(p)
dev.off()

# %%
# saveFinalFigure(
#   plot = p,
#   prefixDir = "/srv/http/betts/hpap/final_figures",
#   fn = "fig1_v3_final",
#   # devices = c("Cairo", "png", "pdf"),
#   devices = c("pdf"),
#   addTimestamp = FALSE,
#   gwidth = 8.5,
#   gheight = 8)
#
# %%

pdf("/srv/http/betts/hpap/final_figures/leg.pdf", width=5, height=5)
print(f_legend)
# figH_leg = ComplexHeatmap::Legend(
#     title = NULL,
#     at = c("pLN", "mLN", "Spleen", "ND", "AAb+", "T1D"),
#     labels = c("pLN", "mLN", "Spleen", "ND", "AAb+", "T1D"),
#     nrow = 1,
#     legend_gp = gpar(col = COLORS$tissue_detailed),
#     legend_height = unit(1, "mm"),
#     legend_width = unit(0.5, "mm"),
#     labels_gp = gpar(fontsize = 5),
#     grid_width = unit(2, "mm"),
#     grid_height = unit(1, "mm"),
#     title_gp = gpar(fontsize = 6, fontface = "plain"),
#     graphics = list(
#                     function(x,y,w,h) grid.points(x, y, gp=gpar(fill="#ffffff", col="#000000"), pch=21, size=unit(3,"mm")),
#                     function(x,y,w,h) grid.points(x, y, gp=gpar(fill="#ffffff", col="#000000"), pch=22, size=unit(3,"mm")),
#                     function(x,y,w,h) grid.points(x, y, gp=gpar(fill="#ffffff", col="#000000"), pch=23, size=unit(3,"mm")),
#                     function(x,y,w,h) grid.rect(x, y, w*0.33, h, gp=gpar(fill=COLORS$disease["ND"], col=COLORS$disease["ND"])),
#                     function(x,y,w,h) grid.rect(x, y, w*0.33, h, gp=gpar(fill=COLORS$disease["AAb+"], col=COLORS$disease["AAb+"])),
#                     function(x,y,w,h) grid.rect(x, y, w*0.33, h, gp=gpar(fill=COLORS$disease["T1D"], col=COLORS$disease["T1D"]))
#                     )
# )
# draw(leg)
# COLORS$disease

dev.off()

# %%
    graphics = lapply(1:length(LinTissueOrder), function(label) function(x,y,w,h) {
                        grid.points(x,y,gp=gpar(col=unname(COLORS$tissue_detailed[label])),pch=16)
    })
label = "Spleen"
unname(COLORS$tissue_detailed[label])

# %%

