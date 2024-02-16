# note that this code is written to be run from the project base directory
# renv::load("/data/hpap-citeseq/hpap-citeseq-analysis")

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
# S8A ridgline of cd38 and cd25 for cd8 t cell
# include plasmablast and cd4 tcm/treg as pos controls, respectively.
################################################################################
labellerAdt <- tsaCatalog$cleanName
names(labellerAdt) <- paste0("adt_", tsaCatalog$DNA_ID)

specClusters <- levels(manualClusterOrder)[grepl("^(CD8|CD4 Tcm/Treg|B plasma)", levels(manualClusterOrder))]
seu_spec_pln <- subset(seu, subset = TissueCondensed == "pLN" & manualAnnot %in% specClusters)

featuresOfInterest <- c(
  "A0389",
  "A0085"
)

tmp <- RidgePlot(seu_spec_pln,
  features = featuresOfInterest, group.by = "manualAnnot", combine = FALSE)

tmp <- lapply(tmp, function(x) {
  x$labels$title <- labellerAdt[x$labels$title]
  
  x <- x +
    stat_summary(aes(group = ident, color = ident),
      fun = median, fun.min = "median", fun.max= "median", geom = "crossbar", linetype = "dotted", show.legend = FALSE) +
    subplotTheme +
    theme(
      plot.margin = margin(3,3,3,3),
      axis.text = element_text(size = BASEPTFONTSIZE),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "blank") 
  
  x$layers[[1]]$geom$default_aes$alpha <- 0.3
  
  return(x)
})



################################################################################
# S8B Heatmap of significant DEGs within module 6 between ND -> AAb+ 
# or ND -> T1D.
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

n_hubs <- 20
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

module6Top <- hub_df %>%
  filter(module == "T1D-M6") %>%
  slice_max(kME, n = 20)

seu_cd8Mem_pln <- subset(
  seu,
  subset = TissueCondensed == "pLN" & manualAnnot %in% c("CD8 Tem/Trm #1", "CD8 Tem/Trm #2", "CD8 Tem/Temra"))

avgExp <- AverageExpression(
  object = seu_cd8Mem_pln,
  assays = "RNA",
  return.seurat = FALSE,
  features = module6Top$lab,
  group.by = c("Disease_Status", "manualAnnot"),
  slot = "data")

avgExp_df <- scale(t(avgExp$RNA)) %>%
  as.data.frame(.) %>%
  mutate(category = rownames(.)) %>%
  separate(category, sep = "_", into = c("disease", "cluster")) %>%
  pivot_wider(names_from = disease, values_from = c(1:length(module6Top$lab))) %>%
  relocate(c("cluster"), c(1))

avgExp_mat <- as.matrix(avgExp_df[, c(2:ncol(avgExp_df))])
rownames(avgExp_mat) <- avgExp_df$cluster

grid_size <- unit(0.2, "cm")

figB <- Heatmap(
  matrix = avgExp_mat,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  column_split = rep(module6Top$lab, each = 3),
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
    `Disease Status` = factor(str_split_fixed(colnames(avgExp_mat), "_", n = 2)[, 2], levels = c("ND", "AAb+", "T1D")),
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
# Final layout and plot all
################################################################################
layout <- c(
  area(1, 1, 6, 12), #a
  area(7, 1, 10, 12) #b
)

p <- wrap_elements(plot = wrap_plots(tmp, nrow = 1)) +
  wrap_elements(full = grid.grabExpr(
    draw(figB,
      heatmap_legend_side = "bottom",
      annotation_legend_side = "bottom",
      merge_legend = TRUE,
      align_heatmap_legend = "global_center",
      background = "transparent",
      padding = unit(c(0.5,1.5,1,0.5), "lines"))
  ), clip = FALSE) +
  plot_annotation(tag_levels = list(LETTERS[1:2])) +
  plot_layout(design = layout) &
  plotTagTheme

saveFinalFigure(
  plot = p,
  prefixDir = "figures/outs",
  fn = "sfig8_final",
  devices = c("pdf", "png"),
  addTimestamp = TRUE,
  gwidth = 8.5,
  gheight = 6)

