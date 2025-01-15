# note that this code is written to be run from the project base directory
# renv::load("/data/hpap-citeseq/hpap-citeseq-analysis")

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
    tmp <- readRDS("outs/rds/seuMergedPostHSP_forFigures_2025-01-12_04-07-24.rds")

    return(tmp)
})

manualClusterOrder <- unique(seu$manualAnnot)
manualClusterOrder <- factor(manualClusterOrder,
  levels = customSortAnnotation(manualClusterOrder),
  labels = stringr::str_trim(customSortAnnotation(manualClusterOrder)))

seu$Disease_Status <- factor(seu$Disease_Status, levels = c("ND", "AAb+", "T1D"))


parentDir <- "figures/greg_flow_data"

################################################################################
# A - cd4 Tn differences in pLN
################################################################################
dfDiseaseScales <- processGregFlowData(paste0(parentDir, "/rds/dfLineageFilter.rds"))

figA <- dfDiseaseScales %>%
  filter(LN_type == "pLN" & cd == "CD4") %>%
  filter(grepl("CD4 Tn$", metric)) %>%
  ggplot(aes(`Disease Status`, value, color = `Disease Status`)) +
  labs(
    y = "% of CD4+ T cells",
    title = "pLN CD4+ Tn") +
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
    plot.title = element_text(size = 8, hjust = 0.5),
    plot.title.position = "panel",
    axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1, color = "#000000"),
    axis.text.y = element_text(size = 6, color = "#000000"),
    axis.title.y = element_text(size = 6))


################################################################################
# B - cd8 Tn differences in pLN
################################################################################
figB <- dfDiseaseScales %>%
  filter(LN_type == "pLN" & cd == "CD8") %>%
  filter(grepl("CD8 Tn$", metric)) %>%
  ggplot(aes(`Disease Status`, value, color = `Disease Status`)) +
  labs(
    y = "% of CD8+ T cells",
    title = "pLN CD8+ Tn") +
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
    plot.title = element_text(size = 8, hjust = 0.5),
    plot.title.position = "panel",
    axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1, color = "#000000"),
    axis.text.y = element_text(size = 6, color = "#000000"),
    axis.title.y = element_text(size = 6))


################################################################################
# Figure C: Heatmap of top 10 genes in module 15
################################################################################
wcgnaCheckpointFile <- "rds/postNetworkPostModule_pLN_ND_T1D_v3.rds"
seuWcgna <- tryCatch(
  {
    print(head(seuWcgna$TissueCondensed))
    message("Seurat object already exists")
    return(seu)
  },
  error = function(cond) {
    message("Seurat object doesn't exist. Loading now.")
    tmp <- readRDS(wcgnaCheckpointFile)
    
    return(tmp)
  })

modules <- GetModules(seuWcgna, seuWcgna@misc$active_wgcna) %>% 
  filter(module != 'grey')

mods <- levels(modules$module)
mods <- mods[mods != 'grey']

mod_colors <- modules %>%
  filter(module %in% mods) %>%
  select(c(module, color)) %>%
  distinct()

n_hubs <- 10
hub_df <- do.call(rbind, lapply(mods, function(cur_mod){
  print(cur_mod)
  cur <- subset(modules, module == cur_mod)
  cur <- cur[,c('gene_name', 'module', paste0('kME_', cur_mod))]
  names(cur)[3] <- 'kME'
  cur <- dplyr::arrange(cur, kME)
  top_genes <- cur %>% dplyr::top_n(n_hubs, wt=kME) %>% .$gene_name
  cur$lab <- ifelse(cur$gene_name %in% top_genes, cur$gene_name, "")
  cur
}))

module14Genes <- hub_df %>%
  filter(module == "T1D-M14") %>%
  arrange(desc(kME))


seu_naiveT_pln <- subset(seu,
  subset = TissueCondensed == "pLN" & manualAnnot %in% c("CD4 naive #1", "CD4 naive #2", "CD8 naive #1", "CD8 naive #2"))

# top 10 heatmap
module14AvgExp <- AverageExpression(
  object = seu_naiveT_pln,
  assays = "RNA",
  return.seurat = FALSE,
  features = module14Genes[c(1:10), "gene_name"],
  group.by = c("Disease_Status", "manualAnnot"),
  slot = "data")

module14_df <- scale(t(module14AvgExp$RNA)) %>%
  as.data.frame(.) %>%
  mutate(category = rownames(.)) %>%
  separate(category, sep = "_", into = c("disease", "cluster")) %>%
  pivot_wider(names_from = disease, values_from = c(1:10)) %>%
  relocate(c("cluster"), c(1))

module14_mat <- as.matrix(module14_df[, c(2:ncol(module14_df))])
rownames(module14_mat) <- module14_df$cluster

fig_module14 <- Heatmap(
  matrix = module14_mat,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  column_split = rep(module14Genes[c(1:10), "gene_name"], each = 3),
  row_names_gp = gpar(fontsize = 6),
  row_dend_width = unit(3, "points"),
  name = "Avg Expression",
  heatmap_legend_param = list(
    direction = "horizontal",
    title_position = "topcenter",
    title_gp = gpar(fontsize = 6, fontface = "plain", hjust = 0.5),
    labels_gp = gpar(fontsize = 6),
    legend_width = unit(2, "cm"),
    legend_height = unit(0.5, "cm")),
  row_title_gp = gpar(fontsize = 8),
  column_title_gp = gpar(fontsize = 4),
  bottom_annotation = ComplexHeatmap::HeatmapAnnotation(
    `Disease Status` = factor(str_split_fixed(colnames(module14_mat), "_", n = 2)[, 2], levels = c("ND", "AAb+", "T1D")),
    col = list(`Disease Status` = COLORS[["disease"]]),
    annotation_name_gp = gpar(fontsize = 6),
    annotation_legend_param = list(
      `Disease Status` = list( 
        title_gp = gpar(fontsize = 6, fontface = "plain"),
        title_position = "topcenter",
        nrow = 1,
        legend_direction = "horizontal", 
        labels_gp = gpar(fontsize = 6)))
  )
)


################################################################################
# Fig D: frequency of naive CD4 in citeseq
################################################################################
# frequency of naive cd4 cells in pln only
cd4Clusters <- levels(manualClusterOrder)[grepl("^CD4", levels(manualClusterOrder))]
seu_cd4_pln <- subset(seu, subset = TissueCondensed == "pLN" & manualAnnot %in% cd4Clusters)

frequencyPlots <- list(
  geom_boxplot(fill = "#00000000", width = 0.8, outlier.shape = NA),
  geom_point(position = position_dodge(width = 0.75), size = 1, stroke = 0.2, alpha = 0.4),
  scale_color_manual(values = COLORS[["disease"]]),
  labs(x = "Cluster", color = ""),
  theme_classic(),
  subplotTheme,
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)),
  theme(
    axis.title.x = element_blank(),
    legend.position = "blank",
    plot.title = element_text(size = 8, hjust = 0.5),
    plot.title.position = "panel",
    axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1, color = "#000000"),
    axis.text.y = element_text(size = 6, color = "#000000"),
    axis.title.y = element_text(size = 6)))

naiveCD4Freq <- data.frame(
  donor = seu_cd4_pln$DonorID,
  disease = seu_cd4_pln$Disease_Status,
  manualAnnot = seu_cd4_pln$manualAnnot
) %>%
  group_by(donor, disease, manualAnnot) %>%
  summarize(nCells = n()) %>%
  mutate(prop = nCells / sum(nCells) * 100) %>%
  filter(grepl("CD4 naive", manualAnnot))


ggpubr::compare_means(data = naiveCD4Freq, prop ~ disease, group.by = "manualAnnot", p.adjust.method = "BH")

fig_naiveCD4Freq <- naiveCD4Freq %>%
  ggplot(aes(x = manualAnnot, y = prop, color = disease)) +
  frequencyPlots +
  labs(y = "% of pLN CD4+ cells")


################################################################################
# Fig E: frequency of naive CD8 in citeseq
################################################################################
# frequency of naive cd4 cells in pln only
cd8Clusters <- levels(manualClusterOrder)[grepl("^CD8", levels(manualClusterOrder))]
seu_cd8_pln <- subset(seu, subset = TissueCondensed == "pLN" & manualAnnot %in% cd8Clusters)

naiveCD8Freq <- data.frame(
  donor = seu_cd8_pln$DonorID,
  disease = seu_cd8_pln$Disease_Status,
  manualAnnot = seu_cd8_pln$manualAnnot
) %>%
  group_by(donor, disease, manualAnnot) %>%
  summarize(nCells = n()) %>%
  mutate(prop = nCells / sum(nCells) * 100) %>%
  filter(grepl("CD8 naive", manualAnnot))

ggpubr::compare_means(data = naiveCD8Freq, prop ~ disease, group.by = "manualAnnot", p.adjust.method = "BH")

fig_naiveCD8Freq <- naiveCD8Freq %>%
  ggplot(aes(x = manualAnnot, y = prop, color = disease)) +
  frequencyPlots +
  labs(y = "% of pLN CD8+ cells")


################################################################################
# Fig F:
# Make violin plots of TXK, FKBP5, ATM for the 4 
# naive clusters (CD4 Naive #1, CD4 Naive #2, CD8 Naive #1, CD4 Naive #2)
################################################################################
naiveGenesOfInterest <- c("TXK", "FKBP5", "ATM")
fig_naiveGenesOfInterest <- Reduce('+', VlnPlot(
  seu_naiveT_pln,
  features = naiveGenesOfInterest,
  group.by = "manualAnnot",
  split.by = "Disease_Status", pt.size = 0, combine = FALSE)) + plot_layout(nrow = 1) &
  scale_fill_manual(values = COLORS[["disease"]]) &
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) &
  geom_boxplot(outlier.shape = NA, size = 0.2, alpha = 0, width = 0.9) &
  subplotTheme &
  theme(
    plot.title = element_text(size = 6, face = "plain", hjust = 0.5),
    plot.title.position = "plot",
    axis.line = element_line(colour = 'black', size = 0.35),
    axis.ticks = element_line(colour = "black", size = 0.25),
    axis.ticks.length = unit(2, "pt"),
    axis.text.x = element_text(size = 4, margin = margin(t = 0)),
    axis.text.y = element_text(size = 4, margin = margin(r = 1)),
    axis.title = element_blank(),
    plot.margin = margin(t = 5, b = 0, l = 3, r = 3),
    legend.position = "none")


################################################################################
# Fig G: cd4 degs
################################################################################

seu_naiveT_pln$Disease_StatusWithManualAnnnot <- paste0(seu_naiveT_pln$Disease_Status, "_", seu_naiveT_pln$manualAnnot)

Idents(seu_naiveT_pln) <- "Disease_StatusWithManualAnnnot"
module14Deg_cd4Naive1 <- FindMarkers(
  object = seu_naiveT_pln,
  ident.1 = "ND_CD4 naive #1",
  ident.2 = "T1D_CD4 naive #1",
  features = module14Genes$gene_name,
  logfc.threshold = 0.1
)

module14Deg_cd4Naive2 <- FindMarkers(
  object = seu_naiveT_pln,
  ident.1 = "ND_CD4 naive #2",
  ident.2 = "T1D_CD4 naive #2",
  features = module14Genes$gene_name,
  logfc.threshold = 0.1
)

fig4DegLollipops <- function(deg, title) {
  deg_final <- deg %>%
    mutate(gene = rownames(.)) %>%
    filter(p_val_adj < 0.05) %>%
    mutate(upregulated = ifelse(avg_log2FC > 0, "ND", "T1D")) %>%
    group_by(upregulated) %>%
    slice_max(abs(avg_log2FC), n = 10) %>%
    arrange(desc(avg_log2FC), .by_group = TRUE)
    
  maxCounts <- deg_final %>%
    summarize(groupCounts = n())
  maxCounts <- ceiling(log10(max(maxCounts$groupCounts))) + 1
  
  deg_final <- deg_final %>%
    mutate(facet_gene_number = cur_group_id() * (10 ** maxCounts) + row_number()) %>%
    mutate(facet_gene_number = paste0(facet_gene_number, "_", gene)) %>%
    mutate(facet_gene_number = factor(facet_gene_number, levels = stringr::str_sort(facet_gene_number, numeric = TRUE)))
  
  p <- deg_final %>%
    ggplot(aes(y = facet_gene_number, x = avg_log2FC, color = upregulated)) +
    geom_point(size = 1) +
    geom_vline(xintercept = 0) +
    scale_color_manual(values = COLORS[["disease"]]) +
    scale_y_discrete(labels = function(x) gsub("\\d+_", "", x)) +
    theme_classic() +
    subplotTheme +
    theme(
      axis.line = element_line(color = "#000000"),
      panel.grid.major.y = element_line(color = "#dddddd"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(size = 8, hjust = 0.5),
      plot.title.position = "panel",
      axis.text = element_text(color = "#000000", size = 6),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "blank") +
    labs(
      y = "Gene",
      x = "Average log2 Fold Change",
      title = title
    )
  
  return(p)
}

fig_cd4Degs <- (fig4DegLollipops(module14Deg_cd4Naive1, "CD4 naive #1") +
  fig4DegLollipops(module14Deg_cd4Naive2, "CD4 naive #2")) / 
  ggplot(data.frame(l = "Average log2 Fold Change", x = 1, y = 1)) +
    geom_text(aes(x, y, label = l), angle = 0, size = 6 / ggplot2::.pt) + 
    theme_void() +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(t = 2.5, b = -2)) +
  plot_layout(heights = c(99, 1))


module14Deg_cd4Naive1 %>% 
  mutate(gene = rownames(.)) %>%
  filter(gene %in% naiveGenesOfInterest)


################################################################################
# Fig H: cd8 degs
################################################################################
module14Deg_cd8Naive1 <- FindMarkers(
  object = seu_naiveT_pln,
  ident.1 = "ND_CD8 naive #1",
  ident.2 = "T1D_CD8 naive #1",
  features = module14Genes$gene_name,
  logfc.threshold = 0.1
)

module14Deg_cd8Naive2 <- FindMarkers(
  object = seu_naiveT_pln,
  ident.1 = "ND_CD8 naive #2",
  ident.2 = "T1D_CD8 naive #2",
  features = module14Genes$gene_name,
  logfc.threshold = 0.1
)

fig_cd8Degs <- (fig4DegLollipops(module14Deg_cd8Naive1, "CD8 naive #1") +
    fig4DegLollipops(module14Deg_cd8Naive2, "CD8 naive #2")) / 
  ggplot(data.frame(l = "Average log2 Fold Change", x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), angle = 0, size = 6 / ggplot2::.pt) + 
  theme_void() +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 2.5, b = -2)) +
  plot_layout(heights = c(99, 1))


module14Deg_cd8Naive1 %>% 
  mutate(gene = rownames(.)) %>%
  filter(gene %in% naiveGenesOfInterest)


findMarkersCombinatorial(
  seuratObj = subset(seu_naiveT_pln, manualAnnot == "CD8 naive #2"),
  combVar = "Disease_Status",
  features = naiveGenesOfInterest
) %>%
  select(-p_val_adj)

################################################################################
# Final layout and plot all
################################################################################
layout <- c(
  patchwork::area(1, 1, 2, 3), # a
  patchwork::area(1, 4, 2, 6), # b representative flow plots
  patchwork::area(1, 7, 2, 12), # c heatmap
  patchwork::area(3, 1, 4, 3), # d citeseq frequency
  patchwork::area(3, 4, 4, 6), # e citeseq frequency
  patchwork::area(3, 7, 4, 12), # f violin
  patchwork::area(5, 1, 7, 6), # g cd4 module 15 degs
  patchwork::area(5, 7, 7, 12) # h cd8 module 15 degs
)

p <- wrap_elements(plot = figA) +
  wrap_elements(plot = figB) +
  wrap_elements(full = grid.grabExpr(
    draw(fig_module14,
      heatmap_legend_side = "bottom",
      annotation_legend_side = "bottom",
      merge_legend = TRUE,
      align_heatmap_legend = "global_center",
      background = "transparent",
      padding = unit(c(0.5,1.5,1,0.5), "lines"))
  ), clip = FALSE) +
  wrap_elements(plot = fig_naiveCD4Freq) +
  wrap_elements(plot = fig_naiveCD8Freq) +
  wrap_elements(full = fig_naiveGenesOfInterest & theme(plot.margin = margin(t = 8, b = 5, l = 5,))) +
  wrap_elements(plot = fig_cd4Degs) +
  wrap_elements(plot = fig_cd8Degs) +
  plot_annotation(tag_levels = list(LETTERS[1:8])) +
  plot_layout(design = layout) &
  plotTagTheme

saveFinalFigure(
  plot = p,
  prefixDir = "figures/outs",
  fn = "fig4_v3_final",
  devices = c("pdf", "png"),
  addTimestamp = TRUE,
  gwidth = 8,
  gheight = 6)
