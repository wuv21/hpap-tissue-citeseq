source("scripts/generic.R")

minimalDimPlotTheme <- theme(
  legend.position = "bottom",
  legend.text = element_text(size = BASEPTFONTSIZE),
  legend.direction = "horizontal",
  legend.justification = "center",
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
  reduction = "rna.umap",
  groupBy,
  groupByTitles = groupBy,
  minimalTheme = TRUE,
  nLegendCols = 5,
  nCols = NULL,
  ...
) {
  if (!(is.vector(groupBy) && is.atomic(groupBy))) {
    groupBy = c(groupBy)
  }
  
  if (is.null(nCols) & length(groupBy) != 1) {
    nCols <- length(groupBy) %/% 2
  } else if (is.null(nCols) & length(groupBy) == 1) {
    nCols <- 1
  }
  
  names(groupByTitles) <- groupBy
  
  listP <- lapply(groupBy, function(x) {
    p <- DimPlot(seu, reduction = reduction, group.by = x, ...) +
      labs(title = groupByTitles[x]) +
      minimalDimPlotTheme +
      guides(color = guide_legend(ncol = nLegendCols))
    
    return(p)
  })
  
  return(wrap_plots(full = listP, ncol = nCols))
}

smallMultipleUmaps <- function(
  seu,
  parameter,
  umapReduction = "rna.umap",
  raster = TRUE,
  devices = c("png"),
  filename = paste0("smallMultipleUmap_", parameter),
  height = 6,
  width = 10,
  ncol = 8,
  titleWrapNChars = 20,
  return_figures = FALSE
) {
  
  colOfInterest <- FetchData(seu, parameter)
  uniqVals <- str_sort(unique(colOfInterest[, parameter]), numeric = TRUE)
  
  smallUmaps <- lapply(uniqVals, function(x) {
    message(paste0("Generating small plot for: ", x))
    
    Idents(seu) <- parameter
    cellsOfInterest <- rownames(colOfInterest)[colOfInterest[, parameter] == x]
    
    p <- DimPlot(
      seu,
      raster = TRUE,
      cells.highlight = cellsOfInterest,
      reduction = umapReduction) +
      labs(title = stringr::str_wrap(x, titleWrapNChars)) +
      theme(
        legend.position = "none",
        axis.text = element_text(size = BASEPTFONTSIZE - 3),
        axis.title = element_blank(),
        plot.title = element_text(size = BASEPTFONTSIZE - 3)
      )
    
    return(p)
  })
  
  if (!return_figures) {
    savePlot(
      smallUmaps,
      fn = filename,
      customSavePlot = patchwork::wrap_plots(smallUmaps, ncol = ncol),
      devices = devices,
      gheight = height,
      gwidth = width)
  } else {
    message("Returning list of plots.")
    return(smallUmaps, ncol = ncol)
  }
}


FeaturePlotCustom <- function(
  seu,
  reduction = "rna.umap",
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

FeaturePlotCustomPaged <- function(
  seu,
  reduction = "rna.umap",
  markers,
  modality,
  graphsPerRow = 4,
  rowsPerPage = 4,
  tsa_catalog,
  ...) {
  
  page_breaks = seq(1, length(markers), (graphsPerRow*rowsPerPage))
  
  plotIdxByPage = lapply(page_breaks, function(pb) {
    r=seq(pb, pb+((graphsPerRow*rowsPerPage)-1), 1)
    if (max(r) > length(markers)) {
      r = seq(pb, length(markers), 1)
    }
    return(r)
  })
  markersByPage = lapply(plotIdxByPage, function(pp) markers[pp])
  lapply(markersByPage, function(m) {
    
    FeaturePlotCustom(
        seu = seu,
        reduction = "rna.umap",
        markers = m,
        modality = modality,
        graphsPerRow = graphsPerRow,
        tsa_catalog = tsa_catalog,
        ...
    )
  })
}

