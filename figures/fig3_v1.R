# note that this code is written to be run from the project base directory

source("figures/genericFigureSettings.R")
source("scripts/dimPlots.R")
source("scripts/deg.R")
library(Seurat)
library(cowplot)
library(WGCNA)
library(hdWGCNA)
library(presto)
library(msigdbr)
library(fgsea)
library(ggtext)

set.seed(42)

parentDir <- "figures/greg_flow_data"

################################################################################
# A - cd4 differences in cd25 expression in pLN
################################################################################
dfDiseaseScales <- processGregFlowData(paste0(parentDir, "/rds/dfLineageFilter.rds"))

figA <- dfDiseaseScales %>%
  filter(LN_type == "pLN" & cd == "CD4") %>%
  filter(grepl("CD25+", metric) & !grepl("HLA-DR+", metric)) %>%
  ggplot(aes(`Disease Status`, value, color = `Disease Status`)) +
  labs(y = "% of CD4+ T cells") +
  geom_boxplot(
    width = 0.8,
    outlier.shape = NA,
    show.legend = FALSE) +
  geom_point(
    size = 1,
    stroke = 0.2,
    alpha = 0.4,
    show.legend = FALSE,
    position = position_jitterdodge(jitter.width = 1, dodge.width = 0.8)) +
  scale_color_manual(values = COLORS$disease) +
  theme_classic() +
  subplotTheme +
  theme(
    axis.title.x = element_blank(),
    strip.background = element_rect(fill = "#00000000", color = "#00000000"),
    strip.text = element_text(size = 8),
    axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1, color = "#000000"),
    axis.text.y = element_text(size = 6, color = "#000000"),
    axis.title.y = element_text(size = 6)) +
  facet_grid(cols = vars(tpop), switch = "y")


################################################################################
# B - cd4 differences in Treg-like
################################################################################
figB <- dfDiseaseScales %>%
  filter(LN_type == "pLN" & cd == "CD4" & metric == "CD4 Mem Tregs") %>%
  ggplot(aes(`Disease Status`, value, color = `Disease Status`)) +
  labs(
    y = "% of Mem CD4+ T cells",
    title = "Treg-like") +
  geom_boxplot(
    width = 0.8,
    outlier.shape = NA,
    show.legend = FALSE) +
  geom_point(
    size = 1,
    stroke = 0.2,
    alpha = 0.4,
    show.legend = FALSE,
    position = position_jitterdodge(jitter.width = 1, dodge.width = 0.8)) +
  scale_color_manual(values = COLORS$disease) +
  theme_classic() +
  subplotTheme +
  theme(
    axis.title.x = element_blank(),
    strip.background = element_rect(fill = "#00000000", color = "#00000000"),
    plot.title = element_text(size = 8, hjust = 0.5),
    plot.title.position = "panel",
    axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1, color = "#000000"),
    axis.text.y = element_text(size = 6, color = "#000000"),
    axis.title.y = element_text(size = 6))


################################################################################
# D - heatmap for treg cluster
################################################################################
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
# generate results for each cluster between nd vs t1d
################################################################################
commonGenesFn <- "outs/csv/pLN_findAllMarkers_byManualAnnot_inMoreThan12Clusters.csv"
if (!file.exists(commonGenesFn)) {
  allMarkers <- lapply(levels(manualClusterOrder), function(x) {
    result <- tryCatch(
      {
        message(paste0("working on: "), x)
        result <- subset(seu, subset = TissueCondensed == "pLN" & manualAnnot == x)
        
      }, error = function(e) {
        result <- NULL
      })
    
    message("doing DE testing now")
    if (!is.null(result) && ncol(result) > 0) {
      Idents(result) <- "Disease_Status"
      df <- FindAllMarkers(
        result,
        assay = "RNA",
        logfc.threshold = 0.1,
        only.pos = TRUE
      )
      
      rownames(df) <- NULL
      
      return(df)
    } else {
      return(NULL)
    }
  })
  
  names(allMarkers) <- levels(manualClusterOrder)
  
  allMarkersDf <- bind_rows(allMarkers, .id = "manualAnnot") %>%
    group_by(manualAnnot, cluster) %>%
    slice_min(p_val_adj, n = 250)
  
  # common genes that are expressed in all clusters for a particular disease state
  commonGenes <- allMarkersDf %>%
    group_by(cluster, gene) %>%
    dplyr::summarize(number = n()) %>%
    arrange(desc(number)) %>%
    group_by(cluster) %>%
    filter(number > 12)
  
  write.csv(commonGenes, file = commonGenesFn,
    row.names = FALSE,
    quote = FALSE)
  
} else {
  commonGenes <- read.csv(commonGenesFn)
}

commonGenesByDisease <- split(commonGenes, commonGenes$cluster)

################################################################################
# actual figure d
################################################################################
seu_pln_treg <- subset(seu, subset = manualAnnot == "CD4 Tcm/Treg" & TissueCondensed == "pLN")
treg_gene_list <- read.csv("figures/fig3_gene_lists/fig3_treg_genelist.csv", header = TRUE, sep = ",")$gene

treg_avg_exp <- AverageExpression(
  object = seu_pln_treg,
  assays = "RNA",
  return.seurat = FALSE,
  features = treg_gene_list,
  group.by = "Disease_Status",
  slot = "data")

treg_avg_exp_rna <- treg_avg_exp$RNA

Idents(seu_pln_treg) <- "Disease_Status"

# From FindMarkers doc: positive values indicate that the gene is more highly
# expressed in the first group
FindMarkers(
  object = seu_pln_treg,
  features = treg_gene_list,
  ident.1 = "ND",
  ident.2 = "AAb+",
  logfc.threshold	= 0
)

fig_treg <- Heatmap(
  matrix = t(scale(t(treg_avg_exp_rna))),
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  column_names_centered = TRUE,
  column_names_gp = gpar(fontsize = 6),
  row_dend_width = unit(3, "points"),
  column_names_rot = 0,
  name = "Scaled Average Expression",
  heatmap_legend_param = list(
    direction = "horizontal",
    title_position = "topcenter",
    title_gp = gpar(fontsize = 6, fontface = "plain"),
    labels_gp = gpar(fontsize = 6),
    grid_height = unit(1.25, "mm"),
    legend_height = unit(2, "mm"),
    legend_width = unit(10, "mm")),
  row_title_gp = gpar(fontsize = 6)
)



################################################################################
# in CD4 Tcm/treg subset, what are differential markers that are expressed?
################################################################################
seu_pln_treg_foxp3 <- subset(seu_pln_treg, subset = FOXP3 > 0.5)
foxp3_diseaseStatus_deg <- findMarkersCombinatorial(seu_pln_treg_foxp3, "Disease_Status")

figE <- plotCombinatorialDEGLollipop(
  foxp3_diseaseStatus_deg %>% filter(!(gene %in% commonGenesByDisease$`ND`$gene) &
      !(gene %in% commonGenesByDisease$`AAb+`$gene) &
      !(gene %in% commonGenesByDisease$`T1D`$gene)), title = "test") +
  coord_cartesian(clip = "off") +
  subplotTheme +
  theme(
    legend.position = "bottom",
    legend.box.margin = margin(t = -10, b = 0),
    axis.text.y = element_text(size = 5, color = "#000000"),
    axis.text.x = element_text(size = 6, color = "#000000"),
    axis.title = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    plot.title = element_blank(),
    strip.text = element_markdown(hjust = 0.5, size = 6)) +
  facet_wrap(~ matchup,
    scales = "free",
    labeller = as_labeller(c(
    "ND_vs_AAb+" = glue("<span style='color: {COLORS$disease['AAb+']};'>up in AAb+</span> | <span style='color: {COLORS$disease['ND']};'>up in ND</span>"),
    "ND_vs_T1D" = glue("<span style='color: {COLORS$disease['T1D']};'>up in T1D</span> | <span style='color: {COLORS$disease['ND']};'>up in ND</span>"),
    "T1D_vs_AAb+" = glue("<span style='color: {COLORS$disease['AAb+']};'>up in AAb+</span> | <span style='color: {COLORS$disease['T1D']};'>up in T1D</span>"))))
  


################################################################################
# potential figure - harmonizome
harmonizome_gene_list <- read.csv("figures/fig3_gene_lists/harmonizome_foxp3_targets.txt", header = FALSE, sep = ",")$V1
harmonizome_gene_list <- intersect(rownames(seu_pln_treg), harmonizome_gene_list)

harmonizome_avg_exp <- AverageExpression(
  object = seu_pln_treg,
  assays = "RNA",
  return.seurat = FALSE,
  features = harmonizome_gene_list,
  group.by = "Disease_Status",
  slot = "data")

harmonizome_avg_exp <- harmonizome_avg_exp$RNA

Idents(seu_pln_treg) <- "Disease_Status"

fig_harmonizome <- Heatmap(
  matrix = t(scale(t(harmonizome_avg_exp))),
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_names = FALSE,
  row_names_gp = gpar(fontsize = 6),
  column_names_centered = TRUE,
  column_names_gp = gpar(fontsize = 6),
  row_dend_width = unit(3, "points"),
  column_names_rot = 0,
  name = "Scaled Average Expression",
  heatmap_legend_param = list(
    direction = "vertical",
    title_position = "topcenter",
    title_gp = gpar(fontsize = 6, fontface = "plain"),
    labels_gp = gpar(fontsize = 6),
    grid_height = unit(1.25, "mm"),
    legend_width = unit(10, "mm")),
  row_title_gp = gpar(fontsize = 6)
)

# From FindMarkers doc: positive values indicate that the gene is more highly
# expressed in the first group
harmonizomeDeg <- findMarkersCombinatorial(
  seuratObj = seu_pln_treg,
  combVar = "Disease_Status",
  assay = "RNA",
  features = harmonizome_gene_list
)

harmonizomeDegNoCommon <- harmonizomeDeg %>%
  filter(!gene %in% commonGenes$gene)

tmp <- harmonizomeDegNoCommon %>%
  group_by(matchup, upregulated) %>%
  filter(p_val_adj_all < 0.05) %>%
  filter(matchup == "ND_vs_AAb+" & upregulated == "ND") %>%
  arrange(desc(avg_log2FC)) %>%
  # slice_max(abs(avg_log2FC), n = 10) %>%
  select(-p_val, -p_val_adj)


################################################################################
# mln (potential figure f)
################################################################################
seu_mln_treg <- subset(seu, subset = manualAnnot == "CD4 Tcm/Treg" & TissueCondensed == "mesLN")

treg_mln_avg_exp <- AverageExpression(
  object = seu_mln_treg,
  assays = "RNA",
  return.seurat = FALSE,
  features = treg_gene_list,
  group.by = "Disease_Status",
  slot = "data")

treg_mln_avg_exp_rna <- treg_mln_avg_exp$RNA

Idents(seu_mln_treg) <- "Disease_Status"


treg_mln_deg <- findMarkersCombinatorial(
  seuratObj = seu_mln_treg,
  combVar = "Disease_Status",
  assay = "RNA",
  features = treg_gene_list
)

treg_mln_deg %>%
  filter(p_val_adj_all < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  select(-p_val, -p_val_adj)



fig_treg_mln <- Heatmap(
  matrix = t(scale(t(treg_mln_avg_exp_rna))),
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  column_names_centered = TRUE,
  column_names_gp = gpar(fontsize = 6),
  row_dend_width = unit(3, "points"),
  column_names_rot = 0,
  name = "Scaled Average Expression",
  heatmap_legend_param = list(
    direction = "vertical",
    title_position = "topcenter",
    title_gp = gpar(fontsize = 6, fontface = "plain"),
    labels_gp = gpar(fontsize = 6),
    grid_height = unit(1.25, "mm"),
    legend_width = unit(15, "mm")),
  row_title_gp = gpar(fontsize = 6)
)





################################################################################
# Final layout and plot all
################################################################################
layout <- c(
  patchwork::area(1, 1, 2, 6), # a
  patchwork::area(3, 1, 4, 2), # b
  patchwork::area(3, 3, 4, 6), # c representative flow plots
  patchwork::area(1, 7, 4, 9), # d heatmap
  patchwork::area(5, 1, 9, 9) # e 
)

p <- wrap_elements(plot = figA) +
  wrap_elements(plot = figB) +
  wrap_elements(full = plot_spacer() + plot_annotation(theme = theme(plot.background = element_rect(fill = "#ffffff")))) +
  wrap_elements(full = grid.grabExpr(
    draw(fig_treg,
      heatmap_legend_side = "bottom",
      align_heatmap_legend = "global_center",
      background = "transparent",
      padding = unit(c(0.5,1.5,0.5,0.5), "lines"))
  ), clip = FALSE) +
  wrap_elements(full = figE) +
  plot_annotation(tag_levels = list(LETTERS[1:5])) +
  plot_layout(design = layout) &
  plotTagTheme

saveFinalFigure(
  plot = p,
  prefixDir = "figures/outs",
  fn = "fig3_final",
  devices = c("png"),
  addTimestamp = TRUE,
  gwidth = 6.5,
  gheight = 6)
  