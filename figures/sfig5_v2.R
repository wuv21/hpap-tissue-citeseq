# note that this code is written to be run from the project base directory
renv::load("/data/hpap-citeseq/hpap-citeseq-analysis")

source("figures/genericFigureSettings.R")
source("scripts/dimPlots.R")
source("scripts/deg.R")
library(Seurat)
library(cowplot)
library(presto)

set.seed(42)

parentDir <- "figures/greg_flow_data"

################################################################################
# A - cd4 differences in cd25 expression in mesln
################################################################################
dfDiseaseScales <- processGregFlowData(paste0(parentDir, "/rds/dfLineageFilter.rds"))

figA <- dfDiseaseScales %>%
  filter(LN_type == "mLN" & cd == "CD4") %>%
  filter(grepl("CD25+", metric) & !grepl("HLA-DR+", metric)) %>%
  ggplot(aes(`Disease Status`, value, color = `Disease Status`)) +
  labs(y = "% of CD4+ T cells", title = "mLN CD25+ CD4+ T cells") +
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
    strip.background = element_rect(fill = "#00000000", color = "#00000000"),
    strip.text = element_text(size = 8),
    axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1, color = "#000000"),
    axis.text.y = element_text(size = 6, color = "#000000"),
    axis.title.y = element_text(size = 6)) +
  facet_grid(cols = vars(tpop), switch = "y")

################################################################################
# b - cd4 differences in cd25 expression in spleen
################################################################################
figB <- dfDiseaseScales %>%
  filter(LN_type == "Spleen" & cd == "CD4") %>%
  filter(grepl("CD25+", metric) & !grepl("HLA-DR+", metric)) %>%
  ggplot(aes(`Disease Status`, value, color = `Disease Status`)) +
  labs(y = "% of CD4+ T cells", title = "Spleen CD25+ CD4+ T cells") +
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
    strip.background = element_rect(fill = "#00000000", color = "#00000000"),
    strip.text = element_text(size = 8),
    axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1, color = "#000000"),
    axis.text.y = element_text(size = 6, color = "#000000"),
    axis.title.y = element_text(size = 6)) +
  facet_grid(cols = vars(tpop), switch = "y")

################################################################################
# C - cd4 differences in Treg-like mesLN
################################################################################
figC <- dfDiseaseScales %>%
  filter(LN_type == "mLN" & cd == "CD4" & metric == "CD4 Mem Tregs") %>%
  ggplot(aes(`Disease Status`, value, color = `Disease Status`)) +
  labs(
    y = "% of Mem CD4+ T cells",
    title = "mLN Treg-like") +
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
# D - cd4 differences in Treg-like spleen
################################################################################
figD <- dfDiseaseScales %>%
  filter(LN_type == "Spleen" & cd == "CD4" & metric == "CD4 Mem Tregs") %>%
  ggplot(aes(`Disease Status`, value, color = `Disease Status`)) +
  labs(
    y = "% of Mem CD4+ T cells",
    title = "Spleen Treg-like") +
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
# E - citeseq treg/tcm phenotypes for mLN
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

seu$Disease_Status <- factor(seu$Disease_Status, levels = c("ND", "AAb+", "T1D"))

seu_mesLN_treg <- subset(seu, subset = manualAnnot == "CD4 Tcm/Treg" & TissueCondensed == "mesLN")
treg_gene_list <- read.csv("figures/fig3_gene_lists/fig3_treg_genelist.csv", header = TRUE, sep = ",")$gene

treg_mesln_avg_exp <- AverageExpression(
  object = seu_mesLN_treg,
  assays = "RNA",
  return.seurat = FALSE,
  features = treg_gene_list,
  group.by = "Disease_Status",
  slot = "data")

treg_mesln_avg_exp_rna <- treg_mesln_avg_exp$RNA

Idents(seu_mesLN_treg) <- "Disease_Status"

# From FindMarkers doc: positive values indicate that the gene is more highly
# expressed in the first group
treg_mesln_deg <- findMarkersCombinatorial(
  seuratObj = seu_mesLN_treg,
  features = treg_gene_list,
  combVar = "Disease_Status",
  logfc.threshold	= 0,
  min.pct = 0.05,
)

treg_mesln_deg_stats <- treg_mesln_deg %>%
  mutate(matchup = factor(matchup, levels = c("ND_vs_AAb+", "T1D_vs_AAb+", "ND_vs_T1D"))) %>%
  mutate(p_val_sym = pValSymnum(p_val_adj_all)) %>%
  group_by(gene) %>%
  arrange(matchup, .by_group = TRUE) %>%
  summarize(final = paste0(p_val_sym, collapse = " / "),
    matches = paste0(matchup, collapse = " / ")) %>%
  mutate(gene = factor(gene, levels = rownames(treg_mesln_avg_exp_rna))) %>%
  arrange(gene)

fig_treg_mesln <- Heatmap(
  matrix = t(scale(t(treg_mesln_avg_exp_rna))),
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
    `tests` = anno_text(treg_mesln_deg_stats$final, gp = gpar(fontsize = 4), show_name = FALSE),
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
# F - citeseq treg/tcm phenotypes for spleen
################################################################################
seu_spleen_treg <- subset(seu, subset = manualAnnot == "CD4 Tcm/Treg" & TissueCondensed == "Spleen")

treg_spleen_avg_exp <- AverageExpression(
  object = seu_spleen_treg,
  assays = "RNA",
  return.seurat = FALSE,
  features = treg_gene_list,
  group.by = "Disease_Status",
  slot = "data")

treg_spleen_avg_exp_rna <- treg_spleen_avg_exp$RNA

Idents(seu_spleen_treg) <- "Disease_Status"

# From FindMarkers doc: positive values indicate that the gene is more highly
# expressed in the first group
treg_spleen_deg <- findMarkersCombinatorial(
  seuratObj = seu_spleen_treg,
  features = treg_gene_list,
  combVar = "Disease_Status",
  logfc.threshold	= 0,
  min.pct = 0.05,
)

treg_spleen_deg_stats <- treg_spleen_deg %>%
  mutate(matchup = factor(matchup, levels = c("ND_vs_AAb+", "T1D_vs_AAb+", "ND_vs_T1D"))) %>%
  mutate(p_val_sym = pValSymnum(p_val_adj_all)) %>%
  group_by(gene) %>%
  arrange(matchup, .by_group = TRUE) %>%
  summarize(final = paste0(p_val_sym, collapse = " / "),
    matches = paste0(matchup, collapse = " / ")) %>%
  mutate(gene = factor(gene, levels = rownames(treg_spleen_avg_exp_rna))) %>%
  arrange(gene)

fig_treg_spleen <- Heatmap(
  matrix = t(scale(t(treg_spleen_avg_exp_rna))),
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
    `tests` = anno_text(treg_spleen_deg_stats$final, gp = gpar(fontsize = 4), show_name = FALSE),
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
# G - flow cd4 tn for mesln
################################################################################
figG <- dfDiseaseScales %>%
  filter(LN_type == "mLN" & cd == "CD4") %>%
  filter(grepl("CD4 Tn$", metric)) %>%
  ggplot(aes(`Disease Status`, value, color = `Disease Status`)) +
  labs(
    y = "% of CD4+ T cells",
    title = "mLN CD4+ Tn") +
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
# H - flow cd4 tn for spleen
################################################################################
figH <- dfDiseaseScales %>%
  filter(LN_type == "Spleen" & cd == "CD4") %>%
  filter(grepl("CD4 Tn$", metric)) %>%
  ggplot(aes(`Disease Status`, value, color = `Disease Status`)) +
  labs(
    y = "% of CD4+ T cells",
    title = "Spleen CD4+ Tn") +
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
# I - cd8 CD25+ Tn differences in mln
################################################################################
figI <- dfDiseaseScales %>%
  filter(LN_type == "mLN" & cd == "CD8") %>%
  filter(grepl("CD25", metric)) %>%
  ggplot(aes(`Disease Status`, value, color = `Disease Status`)) +
  labs(
    y = "% CD25+ of\nCD8+ T cells",
    title = "mLN CD25+ CD8+ T cells") +
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
    plot.title = element_text(size = 8, hjust = 0.5),
    plot.title.position = "panel",
    strip.background = element_rect(color = NA, fill = NA),
    axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1, color = "#FFFFFF"),
    axis.text.y = element_text(size = 6, color = "#000000"),
    axis.title.y = element_text(size = 6, color = "#000000"))


################################################################################
# J - cd8 CD25+ Tn differences in mln
################################################################################
figJ <- dfDiseaseScales %>%
  filter(LN_type == "Spleen" & cd == "CD8") %>%
  filter(grepl("CD25", metric)) %>%
  ggplot(aes(`Disease Status`, value, color = `Disease Status`)) +
  labs(
    y = "% CD25+ of\nCD8+ T cells",
    title = "Spleen CD25+ CD8+ T cells") +
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
    plot.title = element_text(size = 8, hjust = 0.5),
    plot.title.position = "panel",
    strip.background = element_rect(color = NA, fill = NA),
    axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1, color = "#FFFFFF"),
    axis.text.y = element_text(size = 6, color = "#000000"),
    axis.title.y = element_text(size = 6, color = "#000000"))


################################################################################
# K - cd8 cd38+ Tn differences in mln
################################################################################
figK <- dfDiseaseScales %>%
  filter(LN_type == "mLN" & cd == "CD8") %>%
  filter(grepl("CD38", metric)) %>%
  ggplot(aes(`Disease Status`, value, color = `Disease Status`)) +
  labs(
    y = "% CD38+ of\nCD8+ T cells",
    title = "mLN CD38+ CD8+ T cells") +
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
    plot.title = element_text(size = 8, hjust = 0.5),
    plot.title.position = "panel",
    strip.text = element_text(size = 8, color = "#FFFFFF"),
    strip.background = element_rect(color = NA, fill = NA),
    axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1, color = "#000000"),
    axis.text.y = element_text(size = 6, color = "#000000"),
    axis.title.y = element_text(size = 6, color = "#000000"))


################################################################################
# L - cd8 CD38+ Tn differences in mln
################################################################################
figL <- dfDiseaseScales %>%
  filter(LN_type == "Spleen" & cd == "CD8") %>%
  filter(grepl("CD38", metric)) %>%
  ggplot(aes(`Disease Status`, value, color = `Disease Status`)) +
  labs(
    y = "% CD38+ of\nCD8+ T cells",
    title = "Spleen CD38+ CD8+ T cells") +
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
    plot.title = element_text(size = 8, hjust = 0.5),
    plot.title.position = "panel",
    strip.text = element_text(size = 8, color = "#FFFFFF"),
    strip.background = element_rect(color = NA, fill = NA),
    axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1, color = "#000000"),
    axis.text.y = element_text(size = 6, color = "#000000"),
    axis.title.y = element_text(size = 6, color = "#000000"))


################################################################################
# M - cxcr3, tox, gzmk differences in mln
################################################################################
seu_mesln_TemTemra <- subset(seu, subset = manualAnnot == "CD8 Tem/Temra" & TissueCondensed == "mesLN")
fig_vlnCTG_mesln <- data.frame(
  CXCR3 = seu_mesln_TemTemra@assays$RNA@data["CXCR3", ],
  TOX = seu_mesln_TemTemra@assays$RNA@data["TOX", ],
  GZMK = seu_mesln_TemTemra@assays$RNA@data["GZMK", ],
  disease = seu_mesln_TemTemra$Disease_Status
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
  labs(y = "Expression", title = "mLN CD8 Tcm/Tem/Temra") +
  theme_classic() +
  subplotTheme +
  theme(
    legend.position = "blank",
    plot.title = element_text(size = 8, hjust = 0.5),
    plot.title.position = "panel",
    axis.title = element_blank(),
    strip.background = element_rect(color = NA, fill = NA),
    strip.text = element_text(size = 6, color = "#000000"),
    axis.text = element_text(size = 6, color = "#000000"),
    axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1, color = "#000000"),
  )

findMarkersCombinatorial(
  seu_mesln_TemTemra,
  combVar = "Disease_Status",
  features = c("CXCR3", "TOX", "GZMK")
) %>%
  select(-p_val_adj)


################################################################################
# N - cxcr3, tox, gzmk differences in spleen
################################################################################
seu_spleen_TemTemra <- subset(seu, subset = manualAnnot == "CD8 Tem/Temra" & TissueCondensed == "Spleen")
fig_vlnCTG_spleen <- data.frame(
  CXCR3 = seu_spleen_TemTemra@assays$RNA@data["CXCR3", ],
  TOX = seu_spleen_TemTemra@assays$RNA@data["TOX", ],
  GZMK = seu_spleen_TemTemra@assays$RNA@data["GZMK", ],
  disease = seu_spleen_TemTemra$Disease_Status
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
  labs(y = "Expression", title = "Spleen CD8 Tcm/Tem/Temra") +
  theme_classic() +
  subplotTheme +
  theme(
    legend.position = "blank",
    plot.title = element_text(size = 8, hjust = 0.5),
    plot.title.position = "panel",
    axis.title = element_blank(),
    strip.background = element_rect(color = NA, fill = NA),
    strip.text = element_text(size = 6, color = "#000000"),
    axis.text = element_text(size = 6, color = "#000000"),
    axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1, color = "#000000")
  )

findMarkersCombinatorial(
  seu_mesln_TemTemra,
  combVar = "Disease_Status",
  features = c("CXCR3", "TOX", "GZMK")
) %>%
  select(-p_val_adj)

################################################################################
# O - gzmb klrb1 in nk in mesln
################################################################################
seu_mesln_nk <- subset(seu, subset = manualAnnot %in% c("NK", "NK/ILC") & TissueCondensed == "mesLN")
fig_vlnGK_mesln <- data.frame(
  GZMB = seu_mesln_nk@assays$RNA@data["GZMB", ],
  KLRB1 = seu_mesln_nk@assays$RNA@data["KLRB1", ],
  disease = seu_mesln_nk$Disease_Status
) %>%
  pivot_longer(cols = -any_of(c("disease")), names_to = "gene", values_to = "data") %>%
  ggplot(aes(x = disease, y = data)) +
  geom_violin(aes(fill = disease)) +
  geom_jitter(color = "#00000020", size = 0.2, width = 0.3, height = 0) +
  geom_boxplot(fill = "#00000000", outlier.shape = NA, color = "#000000") +
  facet_wrap(~ gene) +
  scale_fill_manual(values = COLORS$disease) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(y = "Expression", title = "mLN NK") +
  theme_classic() +
  subplotTheme +
  theme(
    legend.position = "blank",
    axis.title = element_blank(),
    plot.title = element_text(size = 8, hjust = 0.5),
    plot.title.position = "panel",
    strip.background = element_rect(color = NA, fill = NA),
    strip.text = element_text(size = 6, color = "#000000"),
    axis.text = element_text(size = 6, color = "#000000"),
    axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1, color = "#000000")
  )

################################################################################
# P - gzmb klrb1 in nk in spleen
################################################################################
seu_spleen_nk <- subset(seu, subset = manualAnnot %in% c("NK", "NK/ILC") & TissueCondensed == "Spleen")
fig_vlnGK_spleen <- data.frame(
  GZMB = seu_spleen_nk@assays$RNA@data["GZMB", ],
  KLRB1 = seu_spleen_nk@assays$RNA@data["KLRB1", ],
  disease = seu_spleen_nk$Disease_Status
) %>%
  pivot_longer(cols = -any_of(c("disease")), names_to = "gene", values_to = "data") %>%
  ggplot(aes(x = disease, y = data)) +
  geom_violin(aes(fill = disease)) +
  geom_jitter(color = "#00000020", size = 0.2, width = 0.3, height = 0) +
  geom_boxplot(fill = "#00000000", outlier.shape = NA, color = "#000000") +
  facet_wrap(~ gene) +
  scale_fill_manual(values = COLORS$disease) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(y = "Expression", title = "Spleen NK") +
  theme_classic() +
  subplotTheme +
  theme(
    legend.position = "blank",
    axis.title = element_blank(),
    plot.title = element_text(size = 8, hjust = 0.5),
    plot.title.position = "panel",
    strip.background = element_rect(color = NA, fill = NA),
    strip.text = element_text(size = 6, color = "#000000"),
    axis.text = element_text(size = 6, color = "#000000"),
    axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1, color = "#000000")
  )

################################################################################
# Final layout and plot all
################################################################################
layout <- c(
  patchwork::area(1, 1, 2, 6), # a
  patchwork::area(3, 1, 4, 6), # b
  patchwork::area(5, 1, 6, 3), # c
  patchwork::area(5, 4, 6, 6), # d
  patchwork::area(7, 1, 10, 3), # e
  patchwork::area(7, 4, 10, 6), # f
  patchwork::area(1, 7, 2, 9), # g
  patchwork::area(1, 10, 2, 12), # h
  patchwork::area(3, 7, 4, 12), # i
  patchwork::area(5, 7, 6, 12), # j
  patchwork::area(7, 7, 8, 12), # k
  patchwork::area(9, 7, 10, 12), # l
  patchwork::area(11, 1, 12, 3), # m
  patchwork::area(11, 4, 12, 6), # n
  patchwork::area(11, 7, 12, 9), # o
  patchwork::area(11, 10, 12, 12) # p
)

p <- wrap_elements(plot = figA) +
  wrap_elements(plot = figB) +
  wrap_elements(plot = figC) +
  wrap_elements(plot = figD) +
  wrap_elements(full = grid.grabExpr(
    draw(fig_treg_mesln,
      heatmap_legend_side = "bottom",
      align_heatmap_legend = "global_center",
      background = "transparent",
      column_title = "mLN CD4 Treg/Tcm",
      column_title_gp=grid::gpar(fontsize = 8),
      padding = unit(c(0.5,1.5,0.5,0.5), "lines"))
  ), clip = FALSE) +
  wrap_elements(full = grid.grabExpr(
    draw(fig_treg_spleen,
      heatmap_legend_side = "bottom",
      align_heatmap_legend = "global_center",
      column_title = "Spleen CD4 Treg/Tcm",
      column_title_gp=grid::gpar(fontsize = 8),
      background = "transparent",
      padding = unit(c(0.5,1.5,0.5,0.5), "lines"))
  ), clip = FALSE) +
  wrap_elements(plot = figG) +
  wrap_elements(plot = figH) +
  wrap_elements(plot = figI) +
  wrap_elements(plot = figJ) +
  wrap_elements(plot = figK) +
  wrap_elements(plot = figL) +
  wrap_elements(plot = fig_vlnCTG_mesln) +
  wrap_elements(plot = fig_vlnCTG_spleen) +
  wrap_elements(plot = fig_vlnGK_mesln) +
  wrap_elements(plot = fig_vlnGK_spleen) +
  plot_annotation(tag_levels = list(LETTERS[1:15])) +
  plot_layout(design = layout) &
  plotTagTheme

saveFinalFigure(
  plot = p,
  prefixDir = "figures/outs",
  fn = "sfig5_v2_final",
  devices = c("pdf", "png"),
  addTimestamp = TRUE,
  gwidth = 8,
  gheight = 11)
