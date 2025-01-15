# note that this code is written to be run from the project base directory
# renv::load("/data/hpap-citeseq/hpap-citeseq-analysis")

source("figures/genericFigureSettings.R")
source("scripts/deg.R")
library(ComplexHeatmap)
library(Seurat)
library(cowplot)
library(presto)
library(flowCore) #for working with FCS files
library(flowWorkspace) #for working with FCS files
library(cyCombine) #for turning FCS files into CS
library(ggridges)

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

sfigClusterOrder <- read.csv("figures/sfig4_ordering/cluster_order.csv")

sfigClusterOrder <- factor(sfigClusterOrder$clusterOrder,
  levels = sfigClusterOrder$clusterOrder,
  labels = stringr::str_trim(sfigClusterOrder$clusterOrder))

seu$Disease_Status <- factor(seu$Disease_Status, levels = c("ND", "AAb+", "T1D"))

# one of the clusters has an extra space so just removing it here
seu$manualAnnot <- stringr::str_trim(seu$manualAnnot)

# fix cluster naming
seu$manualAnnot <- gsub("CD8 Tem/Temra", "CD8 Tcm/Tem/Temra", seu$manualAnnot)


################################################################################
# read in fcs files
# code is modified from @greg
################################################################################
extendedFlowDir <- "figures/extended_greg_flow_data/fcs/"
finalRds <- "figures/extended_greg_flow_data/uncorrected.rds"

if (!file.exists(finalRds)) {
  # Channel names
  template <- flowCore::read.FCS(paste0(extendedFlowDir, "HPAP099_Spleen_B cells.fcs"),
    transformation = FALSE, truncate_max_range = FALSE)
  template <-  as.data.frame(template@parameters@data)
  
  # Import MonoDC FCS file, only keep template parameters, then re-export
  monoDCfcs <- flowCore::read.FCS(paste0(extendedFlowDir, "HPAP099_Spleen_MonoDC.fcs"),
    transformation = FALSE, truncate_max_range = FALSE)
  monoDCfcs <- monoDCfcs[,!colnames(monoDCfcs) %in% c("SampleID")]
  
  write.FCS(monoDCfcs, paste0(extendedFlowDir, "HPAP099_Spleen_MonoDC.fcs"))
  
  # Create a flowSet of transformed FCS files, using the biexponential transform in flowJo
  fs <- flowCore::read.flowSet(
    list.files(extendedFlowDir), 
    path = extendedFlowDir, 
    transformation = FALSE,
    truncate_max_range = FALSE)
  
  normchnls <- as.data.frame(template$name[4:27])
  normchnls <- as.character(normchnls[[1]])
  
  biex <- flowWorkspace::flowjo_biexp(widthBasis = -100)
  
  # flowCore based transform method
  tl <- flowCore::transformList(from = normchnls, 
    to = normchnls, 
    tfun = biex, 
    transformationId = "biex transformation")
  fs <- flowCore::transform(fs, tl)
  
  ### Prepare table from FCS files and export it
  # Extract markers from panel
  panel_file <- "figures/extended_greg_flow_data/metadata/panel.csv"
  metadata_file <- "figures/extended_greg_flow_data/metadata/CyCombine_metadata.csv"
  
  markers <- read.csv(panel_file) %>%
    pull(Antigen)
  
  # put all FCS files into an uncorrected tibble format
  uncorrected <- cyCombine::prepare_data(
    flowset = fs,
    metadata = metadata_file, 
    filename_col = "Filename",
    batch_ids = "CellPop",
    condition = "Disease_Status",
    down_sample = FALSE,
    markers = markers,
    transform = FALSE
  )
  
  # remove blank columns
  uncorrected <- uncorrected[-c(2:4, 29)]
  
  # rename columns to match `markers` list
  uncorrected <- uncorrected %>%
    rename(CD11c = CD11cBUV395) %>% 
    rename(`CD279 (PD-1)` = PD1BV421) %>% 
    rename(CD14 = CD14BV480) %>% 
    rename(CD15 = CD15FITC) %>%
    rename(LiveDead = LDAquaBV510) %>%
    rename(CD8 = CD8BUV496) %>%
    rename(CD45RA = CD45RABUV563) %>%
    rename(`CD56 (NCAM)` = CD56BV570) %>%
    rename(CD123 = CD123PE) %>%
    rename(`HLA-DR` = HLADRBV605) %>%
    rename(`CD127 (IL-7Rα)` = CD127PECF594) %>%
    rename(`CD185 (CXCR5)` = CXCR5AF647) %>%
    rename(CD38 = CD38BUV661) %>%
    rename(CD27 = CD27BV650) %>%
    rename(CD69 = CD69PECy5) %>%
    rename(CD16 = CD16BV711) %>%
    rename(CD34 = CD34PECy55) %>%
    rename(CD45 = CD45AF700) %>%
    rename(CD25 = CD25BUV737) %>%
    rename(CD4 = CD4BB790) %>%
    rename(CCR7 = CCR7APCCy7) %>%
    rename(CD21 = CD21PECy7) %>%
    rename(CD19 = CD19BV785) %>%
    rename(CD3 = CD3BUV805)
  
  saveRDS(uncorrected, finalRds)
} else {
  uncorected <- readRDS(finalRds)
}


# matched ADT numbers
adtIds <- tsaCatalog[tsaCatalog$cleanName %in% colnames(uncorrected)[2:25], c("DNA_ID", "cleanName")]
uncorrected <- as.data.frame(uncorrected)

allPlots <- lapply(seq_along(adtIds$DNA_ID), function(i) {
  maxQuantile <- 0.99
  minQuantile <- 0.01
  rescaleBounds <- c(0, 5)
  
  message(paste0("Working on: ", adtIds[i, "cleanName"]))
  
  tmp <- data.frame(
    marker = seu@assays$adt@data[adtIds[i, "DNA_ID"], ],
    phenotypeCondensed = seu$phenotypeCondensed,
    dataType = "CITEseq"
  )
  
  rownames(tmp) <- NULL
  
  tmp <- tmp %>%
    dplyr::mutate(phenotypeCondensed = case_when(
      phenotypeCondensed == "Mono/APC cell" ~ "Mono-DC",
      phenotypeCondensed == "NK/ILC cell" ~ "NK cells",
      TRUE ~ paste0(phenotypeCondensed, "s")
    ))
  
  tmpQuantiles <- quantile(tmp$marker, probs = c(minQuantile, maxQuantile))
  
  tmp <- tmp %>%
    dplyr::filter(marker < tmpQuantiles[2] & marker > tmpQuantiles[1]) %>%
    dplyr::mutate(marker = scales::rescale(marker, to = rescaleBounds))

  tmp2 <- data.frame(
    marker = uncorrected[, adtIds[i, "cleanName"]],
    phenotypeCondensed = uncorrected$batch
  )
  
  tmp2Quantiles <- quantile(tmp2$marker, probs = c(minQuantile, maxQuantile))
  
  tmp2 <- tmp2 %>%
    dplyr::filter(marker < tmp2Quantiles[2] & marker > tmp2Quantiles[1]) %>%
    mutate(dataType = "Flow") %>%
    mutate(marker = scales::rescale(marker, to = rescaleBounds))
  
  p <- bind_rows(tmp, tmp2) %>%
    ggplot(aes(x = marker, y = phenotypeCondensed, fill = dataType)) +
    geom_density_ridges(alpha = 0.6, scale = 1) +
    theme_bw() +
    labs(
      title = adtIds[i, "cleanName"],
      x = "Scaled expression",
      fill = "Method") +
    theme(
      axis.text = element_text(size = BASEPTFONTSIZE, color = "#000000"),
      axis.title.x = element_text(size = BASEPTFONTSIZE, color = "#000000"),
      plot.title = element_text(hjust = 0.5, size = BASEPTFONTSIZE + 2),
      plot.title.position = "panel",
      legend.position = "bottom",
      axis.title.y = element_blank(),
      panel.grid.major.y = element_blank())
  
  return(p)
})


# temporary for innate stuff
library(ggtext)

seuSpleen <- subset(seu, subset = TissueCondensed == "Spleen")
commonGenesFn <- "outs/csv/pLN_findAllMarkers_byManualAnnot_inMoreThan12Clusters.csv"
commonGenes <- read.csv(commonGenesFn)
commonGenesByDisease <- split(commonGenes, commonGenes$cluster)

dcDeg <- findMarkersCombinatorial(
  subset(seuSpleen, subset = manualAnnot == "Monocyte #1"),
  combVar = "Disease_Status"
) %>%
  select(-p_val_adj)


plotCombinatorialDEGLollipop(dcDeg %>%
  dplyr::filter(!(gene %in% commonGenesByDisease$`ND`$gene) &
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
      "T1D_vs_ND" = glue("<span style='color: {COLORS$disease['ND']};'>up in ND</span> | <span style='color: {COLORS$disease['T1D']};'>up in T1D</span>"),
      "T1D_vs_AAb+" = glue("<span style='color: {COLORS$disease['AAb+']};'>up in AAb+</span> | <span style='color: {COLORS$disease['T1D']};'>up in T1D</span>"))))


################################################################################
# Final layout and plot all
################################################################################
saveFinalFigure(
  plot = wrap_plots(allPlots, ncol = 4) / guide_area() + plot_layout(guides = "collect", heights = c(24, 1)),
  prefixDir = "figures/outs",
  fn = "sfigX_v3_final",
  devices = c("pdf", "png"),
  addTimestamp = TRUE,
  gwidth = 8,
  gheight = 10)

