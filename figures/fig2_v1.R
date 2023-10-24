# note that this code is written to be run from the project base directory

source("figures/genericFigureSettings.R")
source("scripts/dimPlots.R")
library(Seurat)
library(cowplot)
library(WGCNA)
library(hdWGCNA)

set.seed(42)

# seu <- tryCatch(
#   {
#     print(head(seu$TissueCondensed))
#     message("Seurat object already exists")
#     return(seu)
#   },
#   error = function(cond) {
#     message("Seurat object doesn't exist. Loading now.")
#     tmp <- readRDS("outs/rds/seuMergedPostHSP_forFigures_2023-09-17_09-03-10.rds")
#     
#     return(tmp)
# })

wcgnaCheckpointFile <- "rds/postNetworkPostModule_pLN_ND_T1D_v2.rds"
# wcgnaCheckpointImage <- "rds/postNetworkPostModule_pLN_ND_T1D_v2.RData"
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
# 
# gg_dend <- dendextend::as.ggdend(as.dendrogram(network$dendrograms[[1]]))
# ggplot(gg_dend, labels = F) + coord_flip() + scale_y_reverse()
# 
figDend <- PlotDendrogram(seuWcgna, main='ND/T1D pLN hdWGCNA Dendrogram')


################################################################################
# modules of interest
################################################################################
PlotKMEs(seuWcgna, ncol = 7, text_size = 3.5)

# code taken and modified from PlotKMEs() in hdWGCNA
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


mods_of_interest <- c(6,7,8,15,16)
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
      size = 8 / ggplot2:::.pt,
      hjust = 0,
      vjust = 1,
      label.size = 0,
      fill = "#ffffff30"
    ) +
    theme(
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(size = 0.5, colour = "#000000")
    )

  return(p)
})

cowplot::plot_grid(plotlist = plot_list_kmes, nrow = 1)

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

tmp <- Heatmap(
  matrix = corMatrix,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_row_names = TRUE,
  row_title_rot = 0,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  row_gap = unit(5, "pt"),
  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
    fdrVal <- fdrMatrix[i,j]
    fdrText <- symnum(
      fdrVal,
      corr = FALSE, na = FALSE,
      cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
      symbols = c("***", "**", "*", " ", ""))
    
    grid.text(fdrText, x, y, gp = gpar(fontsize = 8))
  },
  name = "Module Trait Correlation",
  heatmap_legend_param = list(
    direction = "vertical"),
  row_title_gp = gpar(fontsize = 10)
)
draw(tmp, padding = unit(c(1, 3, 1, 1), "lines"))



################################################################################
# Final layout and plot all
################################################################################
layout <- c(
  area(1, 1, 2, 7), #abc
  area(1, 8, 2, 12), #legend
  area(3, 1, 4, 12), #d
  area(5, 1, 11, 6), #e
  area(5, 7, 11, 12) #f
)

figABC_final <- figABC +
  plot_annotation(tag_levels = list(LETTERS[1:3])) &
  plotTagTheme

p <- wrap_elements(full = figABC_final, ignore_tag = TRUE) +
  wrap_elements(full = figABCLegend, ignore_tag = TRUE, clip = FALSE) +
  wrap_elements(full = figD + theme(plot.margin = margin(l = 10, t = 10, b = 30)), clip = FALSE) +
  wrap_elements(full = 
    grid.grabExpr(
      draw(figE,
        heatmap_legend_side = "bottom",
        annotation_legend_side = "bottom",
        background = "transparent",
        padding = unit(c(0, 1.5, 0.5, -1), "lines"),
        merge_legend = TRUE)), clip = FALSE) +
  wrap_elements(full = 
      grid.grabExpr(
        draw(figF,
          heatmap_legend_side = "bottom",
          annotation_legend_side = "bottom",
          background = "transparent",
          padding = unit(c(0, 1.5, 0.5, -1), "lines"),
          merge_legend = TRUE)), clip = FALSE) +
  plot_annotation(tag_levels = list(c("D", "E", "F"))) +
  plot_layout(design = layout) &
  plotTagTheme

saveFinalFigure(plot = p, prefixDir = "outs", fn = "fig3_citeseq_intro", devices = c("png", "pdf"), gwidth = 8.5, gheight = 11)
