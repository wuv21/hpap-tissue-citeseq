# note that this code is written to be run from the project base directory
renv::load("/data/hpap-citeseq/hpap-citeseq-analysis")

source("figures/genericFigureSettings.R")
source("scripts/dimPlots.R")
source("scripts/deg.R")
source("scripts/generic.R")
library(ComplexHeatmap)
library(dplyr)
library(Seurat)
library(cowplot)
library(WGCNA)
library(hdWGCNA)
library(patchwork)
library(stringr)

set.seed(42)


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

# one of the clusters has an extra space so just removing it here
seuWcgna$manualAnnot <- stringr::str_trim(seuWcgna$manualAnnot)

# fix cluster naming
seuWcgna$manualAnnot <- gsub("CD8 Tem/Temra", "CD8 Tcm/Tem/Temra", seuWcgna$manualAnnot)

# cluster order
manualClusterOrder <- unique(seuWcgna$manualAnnot)
sortedClust <- customSortAnnotation(manualClusterOrder)
manualClusterOrder <- factor(manualClusterOrder,
  levels = sortedClust,
  labels = sortedClust)


################################################################################
# S5A - other wcgna modules
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


mods_of_interest <- c(1:4,7:13, 15, 16)
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
    scale_y_continuous(limits = c(0, NA), expand = expansion(0, 0.05)) +
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
      plot.margin = margin(t = 5, b = 5),
      plot.title = element_text(size = 8, hjust = 0.5),
      plot.title.position = "panel",
      axis.line = element_line(size = 0.5, colour = "#000000")
    )
  
  return(p)
})

################################################################################
# S5B - other wcgna modules heatmaps
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
  area(1, 1, 6, 12), #a
  area(7, 1, 10, 12) #b
)

p <- wrap_elements(full = wrap_plots(plot_list_kmes, nrow = 4)) +
  wrap_elements(full = grid.grabExpr(
    draw(tmp,
      heatmap_legend_side = "bottom",
      annotation_legend_side = "bottom",
      background = "transparent",
      padding = unit(c(1, 3, 1, 4), "lines"))
  ), clip = FALSE)  +
  plot_annotation(tag_levels = list(LETTERS[1:2])) +
  plot_layout(design = layout) &
  plotTagTheme

saveFinalFigure(
  plot = p,
  prefixDir = "figures/outs",
  fn = "sfig4_v3_final",
  devices = c("pdf", "png"),
  addTimestamp = TRUE,
  gwidth = 8.5,
  gheight = 11)

