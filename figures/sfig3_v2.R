# note that this code is written to be run from the project base directory
renv::load("/data/hpap-citeseq/hpap-citeseq-analysis")

source("figures/genericFigureSettings.R")
source("scripts/dimPlots.R")
source("scripts/deg.R")
library(dplyr)
library(circlize)
library(ComplexHeatmap)
library(cluster)
library(forcats)

set.seed(42)


################################################################################
# A-F
# modified code from @Greg
################################################################################
dfLineageFilter <- readRDS("figures/greg_flow_data/rds/dfLineageFilter.rds")

LineageList <- c("Eosinophils", "Neutrophils", "B cells","ILCs", "Imm. Granulocytes",
  "HSCs", "DCs", "Monocytes", "NK cells", "T cells", "CD8+ T cells", "CD8 Tcm", "CD8 Tem",
  "CD8 Temra", "CD8 Tn", "CD8 Tnl", "CD4+ T cells", "CD4 Tcm", "CD4 Tem", "CD4 Temra",
  "CD4 Tn", "CD4 Tnl", "CD4 Mem Tregs", "CD4 Mem gcTfh", "CD4 Mem mTfh")

LinColorsOnly <- c("#BC4749", "#228B22", "#73A942", "#023E8A", "#0777B6", "#00B4D8")

# generate df for plotting lineage cell populations
dfTotalLineageFreq <- filter(dfLineageFilter, metric %in% LineageList)

# these populations were selected based upon similarity across subsets
# (B cells), differences between spleen and all LNs (T cells, CD8 Tn), 
# and subtle differences among LNs (NK cells, ILCs, and CD8 Tcm)
LinPops <- c("B cells", "T cells", "NK cells", "ILCs", "CD8 Tcm", "CD8 Tn")

LinPopsGraphs <- lapply(LinPops, function(l) {
  message(paste("Generating graph for", l, "across tissues"), sep = "")
  lplot <- dfTotalLineageFreq %>%
    filter(metric == as.character(l)) %>%
    mutate(Tissue = str_replace(Tissue, "_", " ")) %>%
    mutate(Tissue = fct_relevel(Tissue, "Spleen", "SMA mLN", "MES mLN", "Head pLN", "Body pLN", "Tail pLN")) %>%
    ggplot(aes(Tissue, value, color = Tissue)) +
    geom_boxplot(width = 0.8,
      outlier.shape = NA,
      show.legend = FALSE) +
    geom_point(
      size = 1,
      stroke = 0.2,
      alpha = 0.4,
      show.legend = FALSE,
      position = position_jitterdodge(jitter.width = 1, 
        dodge.width = 0.8)) +
    scale_color_manual(values = LinColorsOnly) +
    labs(y = "% of CD45+ cells", title = l) +
    theme_classic() +
    subplotTheme +
    theme(
      plot.title = element_text(size = BASEPTFONTSIZE, hjust = 0.5, margin = margin(b = 18)),
      plot.title.position = "panel",
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = BASEPTFONTSIZE, color = "#000000"),
      axis.title.y = element_text(size = BASEPTFONTSIZE)
    )
  
  return(lplot)
})

################################################################################
# G
# modified code from @Greg
################################################################################
# PCA compatible dataframe for the entire dataset
dfPCA <- dfLineageFilter %>%
  dplyr::select(1:Tissue, metric, value) %>%
  pivot_wider(id_cols = `HPAP Donor`:`Tissue`, names_from = metric, values_from = value) %>%
  mutate(Sample_ID = paste("S", as.character(row_number()), sep = "")) %>%
  relocate(Sample_ID, .before = `HPAP Donor`) %>%
  as.data.frame()

rownames(dfPCA) <- dfPCA$Sample_ID


# PCA for total lineage populatinos in LNs only
dfLineagePCAln <- dfPCA %>%
  filter(Tissue != "Spleen") %>%
  dplyr::select(LineageList)


TotalLinPCAln <- prcomp(dfLineagePCAln, scale. = TRUE)
dfTotalLinPCAln <- TotalLinPCAln$x %>%
  as.data.frame() %>%
  mutate(Sample_ID = rownames(dfLineagePCAln)) %>%
  left_join(dfPCA, by = "Sample_ID")

Lin_var_explainedln <- TotalLinPCAln$sdev^2/sum(TotalLinPCAln$sdev^2)

LinColorsOnlyNoSpleen <- c("#228B22", "#73A942", "#023E8A", "#0777B6", "#00B4D8")

pcaPlotTheme <- list(
  theme_classic(),
  subplotTheme,
  guides(color = guide_legend(byrow = TRUE)),
  coord_cartesian(clip = "off"),
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank(),
    legend.margin = margin(l = -5),
    legend.spacing.x = unit(0, "pt"),
    legend.text = element_text(size = BASEPTFONTSIZE - 4),
    legend.position = "right",
    legend.box.margin = margin(0,0,0,0),
    legend.direction = "vertical")
)


TotalLinPCAplotln <- dfTotalLinPCAln %>%
  mutate(Tissue = str_replace(Tissue, "_", " ")) %>%
  mutate(Tissue = fct_relevel(Tissue, "SMA mLN", "MES mLN", "Head pLN", "Body pLN", "Tail pLN")) %>%
  ggplot(aes(PC1, PC2, color = Tissue)) +
  geom_point(size = 0.5) +
  labs(x=paste0("PC1: ",round(Lin_var_explainedln[1]*100,1),"%"),
    y=paste0("PC2: ",round(Lin_var_explainedln[2]*100,1),"%")) +
  scale_color_manual(values = LinColorsOnlyNoSpleen) +
  pcaPlotTheme


################################################################################
# H
# modified code from @Greg
################################################################################
dfPCAfilter <- dfPCA %>%
  filter(!is.na(`CD4 Tcm CD25+`)) 

pcaCD4Pops <- c()
pcaCD4Pops <- append(pcaCD4Pops, str_match(colnames(dfPCA), "^CD4 .*$"))
pcaCD4Pops <- pcaCD4Pops[!is.na(pcaCD4Pops)]
pcaCD4Pops <- pcaCD4Pops %>%
  str_remove("CD4 Tn CD127.CD27.") %>%
  str_remove("CD4 Tn CD27.") %>%
  str_remove("CD4 Tn CD127.") %>%
  str_remove("CD4 Tnl CD127\\+CD27\\+") %>%
  str_remove("^CD4 Mem .*$")
pcaCD4Pops <- pcaCD4Pops[pcaCD4Pops != ""]

dfCD4forPCA <- dplyr::select(dfPCAfilter, all_of(pcaCD4Pops))

CD4PCA <- prcomp(dfCD4forPCA, scale. = TRUE)

dfPCAMetadata <- dplyr::select(dfPCA, Tissue, `HPAP Donor`) 
dfPCAMetadata$Sample_ID <- rownames(dfPCAMetadata)

dfCD4PCA <- CD4PCA$x %>%
  as.data.frame() %>%
  mutate(Sample_ID = rownames(dfCD4forPCA)) %>%
  left_join(dfPCAMetadata, by = "Sample_ID") %>%
  filter(Tissue != "Spleen")

CD4_var_explained <- CD4PCA$sdev^2/sum(CD4PCA$sdev^2)

CD4PCAplot <- dfCD4PCA %>%
  mutate(Tissue = str_replace(Tissue, "_", " ")) %>%
  mutate(Tissue = fct_relevel(Tissue, "SMA mLN", "MES mLN", "Head pLN", "Body pLN", "Tail pLN")) %>%
  ggplot(aes(PC1, PC2, color = Tissue)) +
  geom_point(size = 0.5) +
  labs(x=paste0("PC1: ",round(CD4_var_explained[1]*100,1),"%"),
    y=paste0("PC2: ",round(CD4_var_explained[2]*100,1),"%")) +
  scale_color_manual(values = LinColorsOnlyNoSpleen) +
  pcaPlotTheme

################################################################################
# I
# modified code from @Greg
################################################################################
LinTissueColorsHeatmap <- list(
  `Tissue` = c(
    "Spleen" = "#BC4749", 
    "SMA mLN" = "#228B22",
    "MES mLN" = "#73A942",
    "Head pLN" = "#023E8A",
    "Body pLN" = "#0777B6",
    "Tail pLN" = "#00B4D8"))

# generate a df that is compatible with heatmap generation
dfHeatmapNorm <- dfLineageFilter %>%
  group_by(metric) %>%
  mutate(z_score = scale(value)) %>%
  dplyr::select("HPAP Donor":metric,z_score)

# heatmap of major immune cell population distribution across HPAP tissues
LinTissueOrder <- c("Spleen", "Head pLN", "Body pLN", "Tail pLN", "SMA mLN", "MES mLN")
SampleVariablesLin <- c("HPAP Donor", "Tissue")

dfTotalLineageHeatmap <- dfHeatmapNorm %>%
  pivot_wider(names_from = metric, values_from = z_score) %>%
  dplyr::select(SampleVariablesLin, LineageList) %>%
  mutate(Sample_ID = paste("S", as.character(row_number()), sep = "")) %>%
  relocate(Sample_ID, .before = `HPAP Donor`) %>%
  mutate(Tissue = str_replace_all(Tissue,"_", " "))

dfTotalLineageHeatmap <- as.data.frame(dfTotalLineageHeatmap)
rownames(dfTotalLineageHeatmap) <- dfTotalLineageHeatmap$Sample_ID
dfTotalLineageHeatmap <- dplyr::select(dfTotalLineageHeatmap, `HPAP Donor`:ncol(dfTotalLineageHeatmap))

dfLinHeatmapMetadata <- dplyr::select(dfTotalLineageHeatmap, Tissue)

dfTotalLineageHeatmap <- dplyr::select(dfTotalLineageHeatmap, "Eosinophils":ncol(dfTotalLineageHeatmap))
dfTotalLineageHeatmap <- t(dfTotalLineageHeatmap)

col_fun_LineageHM <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red")) #set maxima at 5th and 95th percentiles
LinAnn <- HeatmapAnnotation(df = dfLinHeatmapMetadata,
  which = "col",
  col = LinTissueColorsHeatmap,
  simple_anno_size = unit(0.1, "in"),
  height = unit(0.1, "in"),
  annotation_name_gp = gpar(fontsize = 8),
  annotation_legend_param = list(
    title = "Tissue",
    at = LinTissueOrder,
    direction = "horizontal",
    nrow = 1,
    labels = LinTissueOrder,
    legend_height = unit(0.3, "in"),
    grid_width = unit(0.1, "in"),
    grid_height = unit(0.1, "in"),
    labels_gp = gpar(fontsize = 4),
    title_gp = gpar(fontsize = 6)
  ))

LinHeatmap <- Heatmap(dfTotalLineageHeatmap,
  name = "Z-score",
  col = col_fun_LineageHM,
  show_column_names = FALSE,
  top_annotation = LinAnn,
  cluster_columns = diana,
  row_km = 3,
  row_title_gp = gpar(fontsize = 10),
  row_title_rot = 0,
  row_names_gp = gpar(fontsize = 4),
  row_dend_width = unit(0.1, "in"),
  column_dend_height = unit(0.2, "in"),
  heatmap_legend_param = list(
    legend_height = unit(0.3, "in"),
    direction = "horizontal",
    grid_width = unit(0.1, "in"),
    grid_height = unit(0.1, "in"),
    labels_gp = gpar(fontsize = 4),
    title_gp = gpar(fontsize = 6)
  ))

################################################################################
# # Final layout and plot all
################################################################################
layout <- c(
  patchwork::area(1, 1, 2, 3),
  patchwork::area(1, 4, 2, 6),
  patchwork::area(3, 1, 4, 3),
  patchwork::area(3, 4, 4, 6),
  patchwork::area(5, 1, 6, 3),
  patchwork::area(5, 4, 6, 6),
  patchwork::area(1, 7, 2, 9), # g
  patchwork::area(1, 10, 2, 12), # h
  patchwork::area(3, 7, 6, 12) # i
)

p <- wrap_elements(plot = LinPopsGraphs[[1]]) +
  wrap_elements(plot = LinPopsGraphs[[2]]) +
  wrap_elements(plot = LinPopsGraphs[[3]]) +
  wrap_elements(plot = LinPopsGraphs[[4]]) +
  wrap_elements(plot = LinPopsGraphs[[5]]) +
  wrap_elements(plot = LinPopsGraphs[[6]]) +
  wrap_elements(plot = TotalLinPCAplotln) +
  wrap_elements(plot = CD4PCAplot) +
  wrap_elements(full = grid.grabExpr(
    draw(LinHeatmap,
      merge_legend = TRUE,
      heatmap_legend_side = "bottom",
      annotation_legend_side = "bottom",
      align_heatmap_legend = "global_center",
      background = "transparent",
      padding = unit(c(0.5,1.5,0.5,0.5), "lines"))
  ), clip = FALSE) +
  plot_annotation(tag_levels = list(LETTERS[1:9])) +
  plot_layout(design = layout) &
  plotTagTheme

saveFinalFigure(
  plot = p,
  prefixDir = "figures/outs",
  fn = "sfig3_v2_final",
  devices = c("pdf", "png"),
  addTimestamp = TRUE,
  gwidth = 8,
  gheight = 6)

