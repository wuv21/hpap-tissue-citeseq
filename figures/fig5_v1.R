# note that this code is written to be run from the project base directory
renv::load("/data/hpap-citeseq/hpap-citeseq-analysis")

source("figures/genericFigureSettings.R")
source("scripts/dimPlots.R")
source("scripts/deg.R")
library(Seurat)
library(cowplot)
library(WGCNA)
library(hdWGCNA)

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

seu$Disease_Status <- factor(seu$Disease_Status, levels = c("ND", "AAb+", "T1D"))


parentDir <- "figures/greg_flow_data"

################################################################################
# A - cd8 CD25+ Tn differences in pLN
################################################################################
dfDiseaseScales <- processGregFlowData(paste0(parentDir, "/rds/dfLineageFilter.rds"))

figA <- dfDiseaseScales %>%
  filter(LN_type == "pLN" & cd == "CD8") %>%
  filter(grepl("CD25", metric)) %>%
  ggplot(aes(`Disease Status`, value, color = `Disease Status`)) +
  labs(
    y = "% CD25+ of\nCD8+ T cells") +
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
  facet_wrap(~ tpop, nrow = 1) +
  scale_color_manual(values = COLORS$disease) +
  theme_classic() +
  subplotTheme +
  theme(
    axis.title.x = element_blank(),
    strip.text = element_text(size = 8, color = "#000000"),
    strip.background = element_rect(color = NA, fill = NA),
    axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1, color = "#FFFFFF"),
    axis.text.y = element_text(size = 6, color = "#000000"),
    axis.title.y = element_text(size = 6, color = "#000000"))

################################################################################
# B - cd8 CD38+ Tn differences in pLN
################################################################################
figB <- dfDiseaseScales %>%
  filter(LN_type == "pLN" & cd == "CD8") %>%
  filter(grepl("CD38", metric)) %>%
  ggplot(aes(`Disease Status`, value, color = `Disease Status`)) +
  labs(
    y = "% CD38+ of\nCD8+ T cells") +
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
  facet_wrap(~ tpop, nrow = 1) +
  scale_color_manual(values = COLORS$disease) +
  theme_classic() +
  subplotTheme +
  theme(
    axis.title.x = element_blank(),
    strip.text = element_text(size = 8, color = "#FFFFFF"),
    strip.background = element_rect(color = NA, fill = NA),
    axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1, color = "#000000"),
    axis.text.y = element_text(size = 6, color = "#000000"),
    axis.title.y = element_text(size = 6, color = "#000000"))


################################################################################
# ??? - effector gene list heatmap
################################################################################
effectorGeneList <- read.csv("figures/fig5_gene_lists/CD8_genelist.csv")

effectorGeneListRna <- effectorGeneList$gene[effectorGeneList$modality == "RNA"]

seu_eff <- subset(seu,
  subset = TissueCondensed == "pLN" & manualAnnot %in% c("CD8 Tcm/Tem/Temra", "CD8 Tem/Trm #1", "CD8 Tem/Trm #2"))

# top 10 heatmap
effAvgExp <- AverageExpression(
  object = seu_eff,
  assays = "RNA",
  return.seurat = FALSE,
  features = effectorGeneListRna,
  group.by = c("Disease_Status", "manualAnnot"),
  slot = "data")

effAvgExp_df <- scale(t(effAvgExp$RNA)) %>%
  as.data.frame(.) %>%
  mutate(category = rownames(.)) %>%
  separate(category, sep = "_", into = c("disease", "cluster")) %>%
  pivot_wider(names_from = disease, values_from = c(1:length(effectorGeneListRna))) %>%
  relocate(c("cluster"), c(1))

effAvgExp_mat <- as.matrix(effAvgExp_df[, c(2:ncol(effAvgExp_df))])
rownames(effAvgExp_mat) <- effAvgExp_df$cluster

grid_size <- unit(0.2, "cm")

fig_effGeneList <- Heatmap(
  matrix = effAvgExp_mat,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  column_split = rep(effectorGeneListRna, each = 3),
  row_names_gp = gpar(fontsize = 6),
  row_dend_width = unit(3, "points"),
  column_title_rot = 90,
  name = "Avg Expression",
  heatmap_legend_param = list(
    direction = "horizontal",
    title_position = "topcenter",
    title_gp = gpar(fontsize = 6, fontface = "plain", hjust = 0.5),
    labels_gp = gpar(fontsize = 6),
    grid_height = grid_size,
    legend_width = unit(2, "cm"),
    legend_height = grid_size),
  row_title_gp = gpar(fontsize = 8),
  column_title_gp = gpar(fontsize = 6),
  bottom_annotation = ComplexHeatmap::HeatmapAnnotation(
    `Disease Status` = factor(str_split_fixed(colnames(effAvgExp_mat), "_", n = 2)[, 2], levels = c("ND", "AAb+", "T1D")),
    col = list(`Disease Status` = COLORS[["disease"]]),
    annotation_name_gp = gpar(fontsize = 6),
    simple_anno_size = grid_size,
    annotation_legend_param = list(
      `Disease Status` = list( 
        title_gp = gpar(fontsize = 6, fontface = "plain"),
        title_position = "topcenter",
        nrow = 1,
        grid_height = grid_size,
        grid_width = grid_size,
        legend_direction = "horizontal", 
        labels_gp = gpar(fontsize = 6)))
  )
)

################################################################################
# ??? - cxcr3, tox, gzmk investigation
################################################################################
seu_TemTemra <- subset(seu_eff, subset = manualAnnot == "CD8 Tcm/Tem/Temra")

fig_vlnCTG <- data.frame(
  CXCR3 = seu_TemTemra@assays$RNA@data["CXCR3", ],
  TOX = seu_TemTemra@assays$RNA@data["TOX", ],
  GZMK = seu_TemTemra@assays$RNA@data["GZMK", ],
  disease = seu_TemTemra$Disease_Status
) %>%
  pivot_longer(cols = -any_of(c("disease")), names_to = "gene", values_to = "data") %>%
  ggplot(aes(x = disease, y = data)) +
  geom_violin(aes(fill = disease)) +
  geom_jitter(color = "#00000020", size = 0.2, width = 0.3, height = 0) +
  geom_boxplot(fill = "#00000000", outlier.shape = NA, color = "#000000") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "#ff0000") +
  facet_wrap(~ gene) +
  scale_fill_manual(values = COLORS$disease) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(y = "Expression") +
  theme_classic() +
  subplotTheme +
  theme(
    legend.position = "blank",
    axis.title = element_blank(),
    strip.background = element_rect(color = NA, fill = NA),
    strip.text = element_text(size = 6, color = "#000000"),
    axis.text = element_text(size = 6, color = "#000000")
    )


findMarkersCombinatorial(
  seu_TemTemra,
  combVar = "Disease_Status",
  features = c("CXCR3", "TOX", "GZMK")
) %>%
  select(-p_val_adj)

seu_TemTemra$CXCR3pos <- seu_TemTemra@assays$RNA@data["CXCR3", ] > 0.5
seu_TemTemra$TOXpos <- seu_TemTemra@assays$RNA@data["TOX", ] > 0.5
seu_TemTemra$GZMKpos <- seu_TemTemra@assays$RNA@data["GZMK", ] > 0.5

sharedMarkers <- data.frame(
  cbc = Cells(seu_TemTemra),
  disease = seu_TemTemra$Disease_Status,
  donor = seu_TemTemra$DonorID,
  CXCR3 = seu_TemTemra$CXCR3pos,
  TOX = seu_TemTemra$TOXpos,
  GZMK = seu_TemTemra$GZMKpos
) %>%
  group_by(donor, disease, CXCR3, TOX, GZMK) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  group_by(donor) %>%
  mutate(proportion = n / sum(n)) %>%
  mutate(rowID = factor(row_number())) %>%
  mutate(nudge = n + max(n) * 0.05) %>%
  mutate(group = paste(CXCR3, TOX, GZMK, sep = "_"))

ggpubr::compare_means(data = sharedMarkers, formula = proportion ~ disease, group.by = "group", p.adjust.method = "BH") %>%
  select(group, group1, group2, p, p.adj, method) %>%
  print(n = 30)

sharedMarkers_boxPlot <- sharedMarkers %>%
  ggplot(aes(x = rowID, y = proportion, fill = disease)) +
  geom_boxplot(outlier.shape = 21) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(values = COLORS$disease) +
  coord_cartesian(clip = "off") +
  labs(y = "Proportion of CD8 Tem/Temra") +
  theme_classic() +
  theme(
    legend.position = "blank",
    axis.text = element_text(size = 6),
    axis.title.y = element_text(size = 6),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank())

sharedMarkers_axisPlot <- sharedMarkers %>%
  dplyr::select(-n) %>%
  pivot_longer(cols = c(CXCR3, TOX, GZMK), names_to = "modals", values_to = "presence") %>%
  mutate(lineGroup = ifelse(presence, rowID, paste(rowID, modals, NA))) %>%
  mutate(colorGroup = ifelse(presence, "#000000", "#CCCCCC")) %>%
  ggplot(aes(x = rowID, y = modals)) +
  geom_point(aes(color = colorGroup), size = 2) +
  geom_line(aes(group = lineGroup)) +
  theme_classic() +
  theme(axis.line = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(size = 6),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major.y = element_line(size = 7, color = "#EFEFEF80")) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_color_identity() +
  coord_cartesian(clip = "off")

sharedMarkers_finalPlot <- sharedMarkers_boxPlot / sharedMarkers_axisPlot + plot_layout(heights = c(6, 1))

seu_TemTemra$phenotype <- paste(seu_TemTemra$CXCR3pos, seu_TemTemra$TOXpos, seu_TemTemra$GZMKpos, sep = "_")
Idents(seu_TemTemra) <- "phenotype"
uniqADTByPhenotype <- FindAllMarkers(
  seu_TemTemra,
  logfc.threshold = 0.1,
  assay = "adt"
)

cNegTPosGPos_uniqAdt <- uniqADTByPhenotype %>%
  filter(p_val_adj < 0.05) %>%
  filter(cluster == "FALSE_TRUE_TRUE") %>%
  left_join(tsaCatalog %>% select(DNA_ID, cleanName), by = c("gene" = "DNA_ID"))

cPosTNegGNeg_uniqAdt <- uniqADTByPhenotype %>%
  filter(p_val_adj < 0.05) %>%
  filter(cluster == "TRUE_FALSE_FALSE") %>%
  left_join(tsaCatalog %>% select(DNA_ID, cleanName), by = c("gene" = "DNA_ID"))

# top 10 heatmap
seu$phenotypeCondensed2 <- seu$phenotypeCondensed
seu$phenotypeCondensed2[names(seu_TemTemra$phenotype)] <- seu_TemTemra$phenotype

combinationRenamer <- function(x, actualNames, sep = "_") {
  xNew <- sapply(x, function(s) {
    
    if (!grepl("(TRUE|FALSE|_)+", s)) {
      return(s)
    } 
    
    tmp <- stringr::str_split_fixed(s, pattern = sep, n = length(actualNames))
    newS <- ""
    for (i in c(1:length(actualNames))) {
      if (tmp[1, i] == "TRUE") {
        newS <- paste0(newS, " ", actualNames[i], "+")
      } else {
        newS <- paste0(newS, " ", actualNames[i], "-")
      }
    }
    
    return(trimws(newS))
  }) 
  
  return(xNew)
}

seu$phenotypeCondensed2 <- combinationRenamer(seu$phenotypeCondensed2, c("CXCR3", "TOX", "GZMK"))
 
cNegTPosGPos_avgExp <- AverageExpression(
  object = seu,
  assays = "adt",
  return.seurat = FALSE,
  features = cNegTPosGPos_uniqAdt$gene,
  group.by = c("phenotypeCondensed2"),
  slot = "data")

tsaCatalogDict <- tsaCatalog$cleanName
names(tsaCatalogDict) <- tsaCatalog$DNA_ID

cNegTPosGPos_adtMat <- t(scale(t(cNegTPosGPos_avgExp$adt)))
rownames(cNegTPosGPos_adtMat) <- tsaCatalogDict[rownames(cNegTPosGPos_adtMat)]

grid_size <- unit(0.2, "cm")

fig_cNegTPosGPos_adt <- Heatmap(
  matrix = cNegTPosGPos_adtMat,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  column_names_rot = 90,
  row_dend_width = unit(5, "points"),
  column_dend_height = unit(5, "points"),
  name = "Avg Expression",
  heatmap_legend_param = list(
    direction = "horizontal",
    title_position = "topcenter",
    title_gp = gpar(fontsize = 6, fontface = "plain", hjust = 0.5),
    labels_gp = gpar(fontsize = 6),
    grid_height = grid_size,
    legend_width = unit(2, "cm"),
    legend_height = grid_size),
  row_title_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6)
)

cPosTNegGNeg_avgExp <- AverageExpression(
  object = seu,
  assays = "adt",
  return.seurat = FALSE,
  features = cPosTNegGNeg_uniqAdt$gene,
  group.by = c("phenotypeCondensed2"),
  slot = "data")

cPosTNegGNeg_adtMat <- t(scale(t(cPosTNegGNeg_avgExp$adt)))
rownames(cPosTNegGNeg_adtMat) <- tsaCatalogDict[rownames(cPosTNegGNeg_adtMat)]

fig_cPosTNegGNeg_adt <- Heatmap(
  matrix = cPosTNegGNeg_adtMat,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  column_names_rot = 90,
  row_dend_width = unit(5, "points"),
  column_dend_height = unit(5, "points"),
  name = "Avg Expression",
  heatmap_legend_param = list(
    direction = "horizontal",
    title_position = "topcenter",
    title_gp = gpar(fontsize = 6, fontface = "plain", hjust = 0.5),
    labels_gp = gpar(fontsize = 6),
    grid_height = grid_size,
    legend_width = unit(2, "cm"),
    legend_height = grid_size),
  row_title_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6)
)

# ggpubr::compare_means(data = sharedMarkers, proportion ~ disease, group.by = c("group"), p.adjust.method = "BH")


uniqRNAByPhenotype <- FindAllMarkers(
  seu_TemTemra,
  logfc.threshold = 0.25,
  assay = "RNA"
)

cNegTPosGPos_uniqRNA <- uniqRNAByPhenotype %>%
  filter(p_val_adj < 0.05) %>%
  filter(cluster == "FALSE_TRUE_TRUE")

cPosTNegGNeg_uniqRNA <- uniqRNAByPhenotype %>%
  filter(p_val_adj < 0.05) %>%
  filter(cluster == "TRUE_FALSE_FALSE")


###
seu_TemTemra$phenotype2 <- combinationRenamer(seu_TemTemra$phenotype, actualNames = c("CXCR3", "TOX", "GZMK"))

cNegTPosGPos_avgExpRNA <- AverageExpression(
  object = seu_TemTemra,
  assays = "RNA",
  return.seurat = FALSE,
  features = cNegTPosGPos_uniqRNA$gene,
  group.by = c("phenotype2"),
  slot = "data")

cNegTPosGPos_rnaMat <- scale(t(cNegTPosGPos_avgExpRNA$RNA))

fig_cNegTPosGPos_rna <- Heatmap(
  matrix = cNegTPosGPos_rnaMat,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  column_names_rot = 90,
  row_dend_width = unit(5, "points"),
  column_dend_height = unit(5, "points"),
  name = "Avg Expression",
  heatmap_legend_param = list(
    direction = "horizontal",
    title_position = "topcenter",
    title_gp = gpar(fontsize = 6, fontface = "plain", hjust = 0.5),
    labels_gp = gpar(fontsize = 6),
    grid_height = grid_size,
    legend_width = unit(2, "cm"),
    legend_height = grid_size),
  row_title_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6)
)


###
cPosTNegGNeg_avgExpRNA <- AverageExpression(
  object = seu_TemTemra,
  assays = "RNA",
  return.seurat = FALSE,
  features = cPosTNegGNeg_uniqRNA$gene,
  group.by = c("phenotype2"),
  slot = "data")

cPosTNegGNeg_rnaMat <- scale(t(cPosTNegGNeg_avgExpRNA$RNA))

fig_cPosTNegGNeg_rna <- Heatmap(
  matrix = cPosTNegGNeg_rnaMat,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  column_names_rot = 90,
  row_dend_width = unit(5, "points"),
  column_dend_height = unit(5, "points"),
  name = "Avg Expression",
  heatmap_legend_param = list(
    direction = "horizontal",
    title_position = "topcenter",
    title_gp = gpar(fontsize = 6, fontface = "plain", hjust = 0.5),
    labels_gp = gpar(fontsize = 6),
    grid_height = grid_size,
    legend_width = unit(2, "cm"),
    legend_height = grid_size),
  row_title_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 6)
)



################################################################################
# Final layout and plot all
################################################################################
layout <- c(
  patchwork::area(1, 1, 4, 6), # a flow 1
  patchwork::area(5, 1, 6, 6), # c heatmap
  patchwork::area(7, 1, 8, 6), # d cxcr3/tox/gzmk
  patchwork::area(1, 7, 4, 12), # g shared marker
  patchwork::area(5, 7, 8, 9), # heatmap of adt differences between CXCR3- tox+ gzmk+ and all other cells
  patchwork::area(5, 10, 8, 12), # heatmap of adt differences between CXCR3+ tox- gzmk- and all other cells
  patchwork::area(9, 1, 12, 6), # j heatmap of cxcr3- tox+ gzmk+ markers
  patchwork::area(9, 7, 12, 12) # k heatmap of cxcr3- tox+ gzmk+ markers
)

p <- wrap_elements(full = figA / figB &
    plot_annotation(tag_levels = list(c("A", "B"))) &
    plotTagTheme, ignore_tag = TRUE) +
  wrap_elements(full = grid.grabExpr(
    draw(fig_effGeneList,
      heatmap_legend_side = "bottom",
      annotation_legend_side = "bottom",
      merge_legend = TRUE,
      align_heatmap_legend = "global_center",
      background = "transparent",
      padding = unit(c(0.5,1.5,1,0.5), "lines"))
  ), clip = FALSE) +
  wrap_elements(plot = fig_vlnCTG) +
  wrap_elements(plot = sharedMarkers_finalPlot) +
  wrap_elements(full = grid.grabExpr(
    draw(fig_cNegTPosGPos_adt,
      heatmap_legend_side = "bottom",
      align_heatmap_legend = "global_center",
      background = "transparent",
      padding = unit(c(0.5,1.5,1,0.5), "lines"))
  ), clip = FALSE) +
  wrap_elements(full = grid.grabExpr(
    draw(fig_cPosTNegGNeg_adt,
      heatmap_legend_side = "bottom",
      align_heatmap_legend = "global_center",
      background = "transparent",
      padding = unit(c(0.5,1.5,1,0.5), "lines"))
  ), clip = FALSE) +
  wrap_elements(full = grid.grabExpr(
    draw(fig_cNegTPosGPos_rna,
      heatmap_legend_side = "bottom",
      align_heatmap_legend = "global_center",
      background = "transparent",
      padding = unit(c(0.5,1.5,1,0.5), "lines"))
  ), clip = FALSE) +
  wrap_elements(full = grid.grabExpr(
    draw(fig_cPosTNegGNeg_rna,
      heatmap_legend_side = "bottom",
      align_heatmap_legend = "global_center",
      background = "transparent",
      padding = unit(c(0.5,1.5,1,0.5), "lines"))
  ), clip = FALSE) +
  plot_annotation(tag_levels = list(LETTERS[3:9])) +
  plot_layout(design = layout) &
  plotTagTheme


saveFinalFigure(
  plot = p,
  prefixDir = "figures/outs",
  fn = "fig5_v2_final",
  devices = c("pdf", "png"),
  addTimestamp = TRUE,
  gwidth = 8,
  gheight = 9)
