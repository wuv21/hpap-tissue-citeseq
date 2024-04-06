# note that this code is written to be run from the project base directory
renv::load("/data/hpap-citeseq/hpap-citeseq-analysis")

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
treg_deg <- findMarkersCombinatorial(
  seuratObj = seu_pln_treg,
  features = treg_gene_list,
  combVar = "Disease_Status",
  logfc.threshold	= 0,
  min.pct = 0.05,
)

treg_deg_stats <- treg_deg %>%
  mutate(matchup = factor(matchup, levels = c("ND_vs_AAb+", "T1D_vs_AAb+", "ND_vs_T1D"))) %>%
  mutate(p_val_sym = pValSymnum(p_val_adj_all)) %>%
  group_by(gene) %>%
  arrange(matchup, .by_group = TRUE) %>%
  summarize(final = paste0(p_val_sym, collapse = " / "),
    matches = paste0(matchup, collapse = " / ")) %>%
  mutate(gene = factor(gene, levels = rownames(treg_avg_exp_rna))) %>%
  arrange(gene)

fig_treg <- Heatmap(
  matrix = t(scale(t(treg_avg_exp_rna))),
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  column_names_centered = TRUE,
  column_names_gp = gpar(fontsize = 6),
  row_dend_width = unit(3, "points"),
  column_names_rot = 45,
  width = unit(1.5, "cm"),
  name = "Scaled Average Expression",
  right_annotation = rowAnnotation(
    `tests` = anno_text(treg_deg_stats$final, gp = gpar(fontsize = 4), show_name = FALSE),
    show_annotation_name = TRUE,
    show_legend = FALSE,
    width = unit(9, "mm"),
    annotation_name_gp = gpar(fontsize = 4),
    annotation_name_rot = 45),
  heatmap_legend_param = list(
    direction = "horizontal",
    title_position = "topcenter",
    title_gp = gpar(fontsize = 4, fontface = "plain"),
    labels_gp = gpar(fontsize = 4),
    grid_height = unit(1.25, "mm"),
    legend_height = unit(2, "mm"),
    legend_width = unit(10, "mm")),
  row_title_gp = gpar(fontsize = 6)
)



################################################################################
# NEW: subcluster cd4 tcm/treg
################################################################################
# seu_pln_treg_sub <- seu_pln_treg
# 
# Idents(seu_pln_treg_sub) <- "manualAnnot"
# seu_pln_treg_sub <- FindNeighbors(seu_pln_treg_sub, assay = "RNA", reduction = "harmonyRNA", dims = 1:30)
# 
# seu_pln_treg_sub <- FindSubCluster(
#   seu_pln_treg_sub,
#   cluster = "CD4 Tcm/Treg",
#   graph.name = "RNA_nn",
#   subcluster.name = "subcluster",
#   resolution = 0.5)
# 
# 
# customCd4Genes <- read.csv("miscellaneous_gene_lists/CD4_genes.csv")
# 
# labellerAdt <- tsaCatalog$cleanName
# names(labellerAdt) <- paste0("adt_", tsaCatalog$DNA_ID)
# 
# subclusterAdtRidges <- RidgePlot(seu_pln_treg_sub,
#   features = customCd4Genes[customCd4Genes$modality == "ADT", "DNA_ID"],
#   group.by = "subcluster",
#   combine = FALSE)
# 
# subclusterAdtRidges <- lapply(subclusterAdtRidges, function(x) {
#   x$labels$title <- labellerAdt[x$labels$title]
#   
#   x <- x +
#     subplotTheme +
#     theme(
#       plot.margin = margin(3,3,3,3),
#       plot.title = element_text(size = BASEPTFONTSIZE),
#       axis.text = element_text(size = BASEPTFONTSIZE),
#       axis.title.y = element_blank(),
#       axis.title.x = element_blank(),
#       legend.position = "blank") 
#   
#   x$layers[[1]]$geom$default_aes$alpha <- 0.3
#   
#   return(x)
# })
# 
# subclusterRnaRidges <- RidgePlot(seu_pln_treg_sub,
#   features = toupper(customCd4Genes[customCd4Genes$modality == "RNA", "geneName"]),
#   group.by = "subcluster",
#   combine = FALSE)
# 
# subclusterRnaRidges <- lapply(subclusterRnaRidges, function(x) {
#   x <- x +
#     subplotTheme +
#     theme(
#       plot.margin = margin(3,3,3,3),
#       plot.title = element_text(size = BASEPTFONTSIZE),
#       axis.text = element_text(size = BASEPTFONTSIZE),
#       axis.title.y = element_blank(),
#       axis.title.x = element_blank(),
#       legend.position = "blank") 
#   
#   x$layers[[1]]$geom$default_aes$alpha <- 0.3
#   
#   return(x)
# })
# 
# subclusterRnaVln <- VlnPlot(seu_pln_treg_sub,
#   features = toupper(customCd4Genes[customCd4Genes$modality == "RNA", "geneName"]),
#   group.by = "subcluster",
#   pt.size = 0.2,
#   combine = TRUE)
# 
# subclusterRnaVln &
#     theme(legend.position = "blank",
#       axis.text = element_text(size = 5),
#       plot.title = element_text(size = 7),
#       axis.title = element_blank())



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
    "ND_vs_T1D" = glue("<span style='color: {COLORS$disease['T1D']};'>up in T1D</span> | <span style='color: {COLORS$disease['ND']};'>up in ND</span>"),
    "T1D_vs_AAb+" = glue("<span style='color: {COLORS$disease['AAb+']};'>up in AAb+</span> | <span style='color: {COLORS$disease['T1D']};'>up in T1D</span>"))))
  


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
  devices = c("pdf", "png"),
  addTimestamp = TRUE,
  gwidth = 6.5,
  gheight = 6)
  