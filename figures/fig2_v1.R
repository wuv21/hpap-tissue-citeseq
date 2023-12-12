# note that this code is written to be run from the project base directory
renv::load("/data/hpap-citeseq/hpap-citeseq-analysis")


source("figures/genericFigureSettings.R")
source("scripts/dimPlots.R")
library(Seurat)
library(cowplot)
library(WGCNA)
library(hdWGCNA)
library(ComplexHeatmap)
library(circlize)
library(multcomp)
library(rstatix)
library(WRS2)
library(forcats)

set.seed(42)

parentDir <- "figures/greg_flow_data"

################################################################################
# A - heatmap of flow differences by disease and tissue
# code is modified from @Greg
################################################################################
dfLineageFilter <- readRDS(paste0(parentDir, "/rds/dfLineageFilter.rds"))

MetaFactors <- c("HPAP Donor", "Disease Status", "Tissue", "cold_ischemia",
  "Sex", "Age", "Staining_Flow_Researcher")
dfMetaCheck <- dfLineageFilter %>%
  dplyr::select(all_of(MetaFactors)) %>%
  filter(Tissue != "Spleen") %>%
  distinct()

dfAllPopsFreqStats <- dfLineageFilter %>%
  pivot_wider(names_from = metric, values_from = value)
colnames(dfAllPopsFreqStats) <- gsub("\\+ \\(\\%T ?cells\\)", "", colnames(dfAllPopsFreqStats))
colnames(dfAllPopsFreqStats) <- gsub(" ", "_", colnames(dfAllPopsFreqStats))

coldIsFreqStats <- dfAllPopsFreqStats %>%
  filter(!is.na(cold_ischemia)) %>%
  filter(Disease_Status == "ND") %>%
  filter(Tissue != "Spleen")
coldIsPops <- colnames(coldIsFreqStats[15:ncol(coldIsFreqStats)])

dfcoldIsCorrStats <- lapply(coldIsPops, function(x) {
  corrx <- cor_test(coldIsFreqStats, x, cold_ischemia, method = "spearman")
}) %>%
  bind_rows()

# for visualizing plots of interest (correlation value >= 0.3 or <= -0.3)
dfcoldIsPops2 <- dfcoldIsCorrStats %>%
  filter(cor >= 0.3 | cor <= -0.3) %>%
  arrange(desc(cor))
coldIsPops2 <- dfcoldIsPops2$var1


# @Greg showed that that cold ischemia time has an impact across disease state 
# on certain immune populations.
# 
# To control for this, ANCOVA will be used to account for cold ischemia time.
# Populations affected by cold ischemia are defined as having a significant p 
# value and an R value >= 0.3 or <= -0.3.
# 
# These thresholds are used because using ANCOVA when a factor has a small 
# effect can be detrimental, as over-correcting can lead to errant results.
# For immune populations that are not significantly affected by cold ischemia,
# robust ANOVA will be used.
# 
# Robust ANOVA trims 10% of the values in the extremes, thus being less
# sensitive to outliers and differences in variance than standard ANOVA.

# ANOVA cannot be computed when values have zero variance,
# hence CD27/CD127 in Tn must be removed.

TissueSepAnovapops <- dfAllPopsFreqStats %>%
  dplyr::select(15:ncol(dfAllPopsFreqStats)) %>%
  colnames()

ignorePops <- "(CD._Tn_CD127.CD27.|CD._Tn_CD27.|CD._Tn_CD127.|CD._Tnl_CD127\\+CD27\\+|CD._Mem_CD.*|^CD._Mem_.*\\+$)"
TissueSepAnovapops <- TissueSepAnovapops[!grepl(ignorePops, TissueSepAnovapops)]

dfAllPopsFreqStatsANOVA <- dfAllPopsFreqStats %>%
  pivot_longer(cols = 15:ncol(dfAllPopsFreqStats),
    names_to = "metric",
    values_to = "value")
alltissues <- unique(dfAllPopsFreqStats$Tissue)

# tissues pooled by general anatomical location (pLN/mesLN/spleen) 
LN_types <- unique(dfAllPopsFreqStatsANOVA$LN_type)

# modified fx from @Greg
run_ColisANOVA_withSubset <- function(
  anovapops,
  secondSubsetName,
  secondSubsetValues,
  df,
  selectpops
) {
  allStats <- lapply(anovapops, function(x) {
    subsetStats <- lapply(secondSubsetValues, function(s) {
      message(paste0("testing: ", s, " and ", x))
      
      dfAnova <- df %>%
        dplyr::filter(!!as.symbol(secondSubsetName) == s & metric == x) %>%
        mutate(Disease_Status = as_factor(Disease_Status))
      
        # message(paste0("Disease Status group is recorded in this order: ", paste0(levels(dfAnova$Disease_Status), collapse = " ")))
      
      if (x %in% selectpops) {
        pval.aov <- glht(
          aov(value ~ cold_ischemia + Disease_Status, data = dfAnova), 
          linfct = mcp(`Disease_Status` = "Tukey"))
        
        tablepval <- summary(pval.aov)
        
        comparisonOrder <- names(tablepval[[10]]$coefficients)
        comparisonOrder <- sapply(comparisonOrder, function(x) {
          tmp <- strsplit(x, " - ")[[1]]
          return(paste0(sort(tmp), collapse = "_"))
        })
        
        tablepval <- tablepval[[10]]$pvalues
        testUsed <- "ancova"
        
      } else {
        pval.aov <- lincon(value ~ Disease_Status, data = dfAnova, tr = 0.1)
        tablepval <- pval.aov[1]$comp
        
        comp1 <- as.character(levels(dfAnova$Disease_Status)[tablepval[, 1]])
        comp2 <- as.character(levels(dfAnova$Disease_Status)[tablepval[, 2]])
        comparisonOrder <- sapply(seq_along(comp1), function(i) {
          return(paste0(sort(c(comp1[i], comp2[i])), collapse = "_"))
        })
        
        tablepval <- tablepval[, 6]
      
        testUsed <- "robust anova"
      }
      
      statDf <- data.frame(
        immune_pop = x,
        secondSubset = s,
        pvalTest = testUsed
      )
      
      statDf[, comparisonOrder] <- tablepval
      
      return(statDf)
    })
    
    return(bind_rows(subsetStats))
  })
  
  return(bind_rows(allStats))
}

# by ln type
dfDiseasePoolStats <- run_ColisANOVA_withSubset(
  anovapops = TissueSepAnovapops,
  secondSubsetName = "LN_type",
  secondSubsetValues = unique(dfAllPopsFreqStatsANOVA$LN_type),
  df = dfAllPopsFreqStatsANOVA,
  selectpops = coldIsPops2
)

# tissues separated
dfDiseaseTissuesSepStats <- run_ColisANOVA_withSubset(
  anovapops = TissueSepAnovapops,
  secondSubsetName = "Tissue",
  secondSubsetValues = alltissues,
  df = dfAllPopsFreqStatsANOVA,
  selectpops = coldIsPops2
)

dfNormTissuePop <- dfLineageFilter %>%
  group_by(metric, LN_type) %>%
  mutate(z_score = scale(value)) %>%
  dplyr::select("HPAP Donor":metric, z_score)

pValSymnum <- function(x, showNs = TRUE) {
  tmp <- sapply(x, function(y) {
    if (is.na(y)) {
      return(NA)
    }
    
    if (y < 0.001) {
      return("***")
    } else if (y >= 0.001 & y < 0.01) {
      return("**")
    } else if (y >= 0.01 & y < 0.05) {
      return("*")
    } else {
      return(ifelse(showNs, "ns", ""))
    }
  })
  
  return(tmp)
}

dfsigPopsDisease <- dfDiseasePoolStats %>%
  rename(tissue = secondSubset) %>%
  mutate(immune_pop = str_replace_all(immune_pop, "_", " ")) %>%
  mutate(across(all_of(c("AAb+_T1D", "ND_T1D", "AAb+_ND")),
    .fns = list(pval = ~ pValSymnum(.)),
    .names = "{fn}{col}")) %>%
  filter_at(all_of(c("AAb+_T1D", "ND_T1D", "AAb+_ND")), all_vars(!is.na(.))) %>%
  filter(!if_all(all_of(c("AAb+_T1D", "ND_T1D", "AAb+_ND")), function(x) x > 0.05))
  
dfsigPopsSamples <- dfsigPopsDisease %>%
  dplyr::select(immune_pop, tissue) %>%
  filter(!grepl("CD127.CD27.", immune_pop)) %>% # for large heatmap, CD27 and CD127 will be separated for viewing ease
  filter(!grepl("Tregs$", immune_pop)) # for large heatmap, Tregs represented by CD4+ Tem CD25+ (mostly)

dfsigPopsHeatmap <- dfNormTissuePop %>%
  group_by(metric, `Disease Status`, LN_type) %>%
  summarise(zMean = mean(z_score)) %>%
  pivot_wider(names_from = `Disease Status`, values_from = zMean) %>% 
  rename("immune_pop" = "metric") %>%
  rename("tissue" = "LN_type") %>%
  right_join(dfsigPopsSamples) %>%
  ungroup() %>%
  mutate(Sample_ID = paste("S", as.character(row_number()), sep = "")) %>%
  relocate(Sample_ID, .before = immune_pop) %>%
  left_join(dfsigPopsDisease) %>%
  relocate(pvalND_T1D:`pvalAAb+_T1D`, .before = ND_T1D) %>%
  as.data.frame()

rownames(dfsigPopsHeatmap) <- dfsigPopsHeatmap$Sample_ID
dfsigPopsHeatmap <- dplyr::select(dfsigPopsHeatmap, immune_pop:ncol(dfsigPopsHeatmap))
immune_pops_names <- structure(dfsigPopsHeatmap$immune_pop, names = rownames(dfsigPopsHeatmap))

dfsigPopsHeatmapPlot <- dfsigPopsHeatmap %>%
  dplyr::select(ND, `AAb+`, T1D) %>%
  as.matrix()

col_fun_DiseaseHm <- colorRamp2(c(-0.75, 0, 0.75), c("blue", "white", "red")) 
DiseaseAnn <- HeatmapAnnotation(
  `LNs` = dfsigPopsHeatmap$tissue,
  `NDAAb` = dfsigPopsHeatmap$`pvalAAb+_ND`,
  `AAbT1D` = dfsigPopsHeatmap$`pvalAAb+_T1D`,
  `NDT1D` = dfsigPopsHeatmap$`pvalND_T1D`,
  which = "row",
  col = list(
    `LNs` = COLORS$tissue,
    `NDT1D` = COLORS$`pval-heatmap`, 
    `NDAAb` = COLORS$`pval-heatmap`,
    `AAbT1D` = COLORS$`pval-heatmap`),
  simple_anno_size = unit(0.1, "in"),
  show_annotation_name = TRUE,
  annotation_label = c(
    `LNs` = "Tissue",
    `NDAAb` = "ND v AAb+", 
    `AAbT1D` = "AAb+ v T1D",
    `NDT1D` = "ND v T1D"),
  annotation_name_side = "top",
  annotation_name_rot = 45,
  annotation_name_offset = unit(0.05, "in"),
  annotation_name_gp = gpar(fontsize = 6),
  gap = unit(c(0.1, 0, 0), "in"),
  show_legend = c(`LNs` = TRUE, `NDAAb` = TRUE,
    `AAbT1D` = FALSE, `NDT1D` = FALSE),
  annotation_legend_param = list(
    `LNs` = list(title = "Tissue",
      at = names(COLORS$tissue),
      labels = names(COLORS$tissue),
      title_position = "topcenter",
      grid_width = unit(0.1, "in"),
      grid_height = unit(0.1, "in"),
      labels_gp = gpar(fontsize = 6),
      title_gp = gpar(fontsize = 6)),
    `NDAAb` = list(
      title = "p-value",
      at = c("ns", "*", "**", "***"),
      labels = c("ns", "*", "**", "***"),
      ncol = 2,
      title_position = "topcenter",
      grid_width = unit(0.1, "in"),
      grid_height = unit(0.1, "in"),
      labels_gp = gpar(fontsize = 6),
      title_gp = gpar(fontsize = 6)) 
  ))

figA <- Heatmap(
  dfsigPopsHeatmapPlot,
  name = "Mean Z-score",
  col = col_fun_DiseaseHm,
  cluster_columns = FALSE,
  show_column_names = TRUE,
  column_names_side = "top",
  column_names_rot = 0,
  column_names_centered = TRUE,
  column_names_gp = gpar(fontsize = 6, fontface = "plain"),
  row_km = 3,
  row_gap = unit(0.1, "in"),
  right_annotation = DiseaseAnn,
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 6),
  row_labels = immune_pops_names,
  row_names_gp = gpar(fontsize = 6),
  row_dend_width = unit(0.15, "in"),
  heatmap_legend_param = list(
    legend_direction = "horizontal",
    legend_width = unit(0.6, "in"),
    title_position = "topcenter",
    grid_height = unit(0.1, "in"),
    labels_gp = gpar(fontsize = 6),
    title_gp = gpar(fontsize = 6)
  ))


################################################################################
# load seurat object for downstream...
################################################################################

wcgnaCheckpointFile <- "rds/postNetworkPostModule_pLN_ND_T1D_v2.rds"
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

manualClusterOrder <- unique(seuWcgna$manualAnnot)
manualClusterOrder <- factor(manualClusterOrder,
  levels = customSortAnnotation(manualClusterOrder),
  labels = stringr::str_trim(customSortAnnotation(manualClusterOrder)))


################################################################################
# dendrogram
################################################################################
# network <- GetNetworkData(seuWcgna, seuWcgna@misc$active_wgcna)
# modules <- GetModules(seuWcgna, seuWcgna@misc$active_wgcna)

# gg_dend <- dendextend::as.ggdend(as.dendrogram(network$dendrograms[[1]]))
# ggplot(gg_dend, labels = F) + coord_flip() + scale_y_reverse()
# 


# figDend <- WGCNA::plotDendroAndColors(
#   network$dendrograms[[1]],
#   as.character(modules$color),
#   groupLabels = "Module colors",
#   dendroLabels = FALSE,
#   hang = 0.03,
#   addGuide = TRUE,
#   guideHang = 0.05,
#   main = ""
# )
# gridGraphics::grid.echo()
# figB <- gridGraphics::grid.grab

# PlotDendrogram(seuWcgna, main='ND/T1D pLN hdWGCNA Dendrogram')



################################################################################
# modules of interest
################################################################################

# code taken and modified from PlotKMEs() in hdWGCNA
modules <- GetModules(seuWcgna, seuWcgna@misc$active_wgcna) %>% 
  filter(module != 'grey')

mods <- levels(modules$module)
mods <- mods[mods != 'grey']

mod_colors <- modules %>%
  filter(module %in% mods) %>%
  dplyr::select(c(module, color)) %>%
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


mods_of_interest <- c(6,7,8,15)
mods_of_interest <- paste0("T1D-M", mods_of_interest)

plot_list_kmes <- lapply(mods_of_interest, function(x) {
  cur_color <- subset(mod_colors, module == x) %>% .$color
  cur_df <- subset(hub_df, module == x)
  top_genes <- cur_df %>% 
    dplyr::top_n(n_hubs, wt=kME) %>% 
    arrange(desc(kME))
  
  p <- cur_df %>% ggplot(aes(x = reorder(gene_name, kME), y = kME)) +
    geom_point(size = 1, color = cur_color, fill = cur_color) +
    ggtitle(x) +
    coord_cartesian(clip = "off", expand = FALSE) +
    scale_y_continuous(limits = c(0.1, NA), expand = expansion(0, 0.05)) +
    annotate(
      geom = "label",
      x = -Inf,
      y = Inf,
      label = paste0(top_genes$gene_name, collapse="\n"),
      size = 4 / ggplot2:::.pt,
      hjust = 0,
      vjust = 1,
      label.size = 0,
      fill = "#ffffff30"
    ) +
    # annotate(geom = "rect", xmin = 0, ymin = 0.1, ymax = 1, xmax = 10, fill = "#33333380") +
    labs(
      x = "Gene rank"
    ) +
    subplotTheme +
    theme(
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_text(size = 6),
      axis.title.y = element_text(size = 6),
      axis.text.y = element_text(size = 6, color = "#000000"),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(size = 8, hjust = 0.5),
      plot.title.position = "panel",
      axis.line = element_line(size = 0.5, colour = "#000000")
    )

  return(p)
})

################################################################################
# modules of interest heatmap
################################################################################
seuWcgna <- ModuleTraitCorrelation(
  seuWcgna,
  traits = c("Disease_Status"),
  group.by = "manualAnnot"
)

moduleTraitCorRes <- GetModuleTraitCorrelation(seuWcgna)

corMatrix <- do.call("rbind", moduleTraitCorRes$cor)[, mods_of_interest]
fdrMatrix <- do.call("rbind", moduleTraitCorRes$fdr)[, mods_of_interest]

rownames(corMatrix) <- gsub("all_cells", "Aggregate (Total cells)", rownames(corMatrix))
rownames(fdrMatrix) <- gsub("all_cells", "Aggregate (Total cells)", rownames(fdrMatrix))

tmp <- Heatmap(
  matrix = corMatrix,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_row_names = TRUE,
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 6),
  column_names_rot = 0,
  column_names_centered = TRUE,
  column_names_gp = gpar(fontsize = 6, fontface = "plain"),
  row_gap = unit(5, "pt"),
  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
    fdrVal <- fdrMatrix[i,j]
    fdrText <- pValSymnum(fdrVal)
    
    if(fdrText == "ns") {
      fdrText <- ""
    }
    
    grid.text(fdrText, x, y, gp = gpar(fontsize = 6))
  },
  name = "Module Trait Correlation",
  heatmap_legend_param = list(
    legend_direction = "horizontal",
    legend_width = unit(0.6, "in"),
    title_position = "topcenter",
    grid_height = unit(0.1, "in"),
    labels_gp = gpar(fontsize = 6),
    title_gp = gpar(fontsize = 6))
)


################################################################################
# Final layout and plot all
################################################################################
layout <- c(
  patchwork::area(1, 1, 10, 4), # a
  patchwork::area(1, 5, 2, 12), # b dendrogram
  patchwork::area(3, 5, 10, 12) # d correlation
)


p <- wrap_elements(full = grid.grabExpr(
  draw(figA,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom",
    align_heatmap_legend = "global_center",
    background = "transparent",
    padding = unit(c(0, 0.75, 0.5, 0.5), "lines"),
    merge_legend = TRUE)
), clip = FALSE) +
  wrap_elements(full = wrap_plots(plot_list_kmes, nrow = 1)) +
  wrap_elements(full = grid.grabExpr(
    draw(tmp,
      heatmap_legend_side = "bottom",
      annotation_legend_side = "bottom",
      background = "transparent",
      padding = unit(c(1, 3, 1, 4), "lines"))
  ), clip = FALSE)  +
  plot_annotation(tag_levels = list(LETTERS[1:3])) +
  plot_layout(design = layout) &
  plotTagTheme

saveFinalFigure(
  plot = p,
  prefixDir = "figures/outs",
  fn = "fig2_final",
  devices = c("pdf", "png"),
  addTimestamp = TRUE,
  gwidth = 8.5,
  gheight = 6.75)

