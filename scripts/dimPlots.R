BASEPTFONTSIZE <- 8
BASEFONTSIZE <- BASEPTFONTSIZE / ggplot2:::.pt


minimalDimPlotTheme <- theme(
  legend.position = "bottom",
  legend.text = element_text(size = BASEPTFONTSIZE),
  legend.title = element_blank(),
  axis.title = element_blank(),
  axis.text = element_text(size = BASEPTFONTSIZE),
  legend.key.size = unit(BASEPTFONTSIZE * 1.1, 'points'),
  legend.background = element_rect(fill = "transparent", colour = NA),
  panel.background = element_rect(fill = "transparent", colour = NA),
  plot.background = element_rect(fill = "transparent", colour = NA)
)

DimPlotCustom <- function(
  seu,
  reduction = "wnn.umap",
  groupBy,
  groupByTitles = groupBy,
  minimalTheme = TRUE,
  nCols = length(groupBy) %/% 2
) {
  if (!(is.vector(groupBy) && is.atomic(groupBy))) {
    groupBy = c(groupBy)
  }
  
  names(groupByTitles) <- groupBy
  
  listP <- lapply(groupBy, function(x) {
    p <- DimPlot(seu, reduction = 'wnn.umap', group.by = x) +
      labs(title = groupByTitles[x]) +
      minimalDimPlotTheme
    
    return(p)
  })
  
  return(wrap_plots(full = listP, ncol = nCols))
}

FeaturePlotCustom <- function(
  seu,
  reduction = "wnn.umap",
  markers,
  modality,
  graphsPerRow = 6,
  tsa_catalog,
  ...
) {
  
  if (modality == "adt") {
    featColors <- c("lightgrey", "darkgreen")
  } else {
    featColors <- c("lightgrey", "purple")
  }
  
  
  featurePs <- FeaturePlot(
    seu,
    features = markers,
    cols = featColors,
    combine = FALSE,
    ...)
  
  featurePs <- lapply(featurePs, function(x) {
    if (modality == "adt") {
      tsaID <- x$labels$title
      tsaID <- gsub(paste0(modality, "_"), "", tsaID)
      tsaCDName <- tsa_catalog[tsa_catalog$DNA_ID == tsaID, "cleanName"]
      
      x <- x +
        labs(title = glue("{tsaCDName} ({tsaID})"))
    }
  
    return(x +
      theme(
        legend.key.size = unit(BASEPTFONTSIZE * 1.1, 'points'),
        legend.text = element_text(size = BASEPTFONTSIZE),
        axis.text = element_text(size = BASEPTFONTSIZE),
        axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(size = BASEPTFONTSIZE)))
  })
  
  return(wrap_plots(full = featurePs, ncol = graphsPerRow))
}