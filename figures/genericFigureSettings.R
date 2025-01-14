library(ggplot2)
library(patchwork)
library(png)
library(ggrastr)
library(glue)
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(ComplexHeatmap))

COLORS <- list(
  "disease" = c("ND" = "#949494", "AAb+" = "#2E96FF", "T1D" = "#FF8B3D"),
  "tissue" = c("pLN" = "#0777B6", "mLN" = "#73A942", "Spleen" = "#BC4749"),
  "tissue_detailed" = 
    c("Spleen" = "#BC4749",
      "SMA mLN" = "#228B22",
      "MES mLN" = "#73A942",
      "Head pLN" = "#023E8A",
      "Body pLN" = "#0777B6",
      "Tail pLN" = "#00B4D8"),
  "expanded_palette" = c(
    '#D0D586', '#687FE0','#E2B2DD','#7B9280','#7AE551','#85ADE0','#DE5AB9',
    '#C4E8C5','#D6DF4F','#E35D76','#8B6684','#CFD9E2','#A965DB','#C88BD8',
    '#DEB5A6','#DF9151','#72DE9A','#D846E3','#743BDF','#72DCDF'),
  "pval-heatmap" = c("ns" = "#dadaeb", "*" = "#b3cde3", "**" = "#8c96c6", "***" = "#88419d")
)

plotTagTheme <- theme(
  plot.tag = element_text(size = 10, face = "plain", family = "sans"))

subplotTheme <- theme(
  plot.title.position = "plot",
  plot.title = element_text(size = 12, margin = margin(0,0,0,0), hjust = 0),
  plot.margin = unit(c(0,0,0,0), "pt"),
  panel.background = element_rect(fill = "transparent", colour = NA_character_),
  plot.background = element_rect(fill = "transparent", colour = NA_character_),
  # text = element_text(family = "sans", size = 8),
  axis.title = element_text(family = "sans", size = 8),
  rect = element_rect(fill = "transparent", colour = NULL)
)

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

umapPlotThemeNoLeg <- list(
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    plot.title.position = "panel",
    plot.title = element_text(size = 8, hjust = 0.5, face = "plain"))
)

customSortAnnotation <- function(x) {
  priority <- c("B", "CD4", "CD8", "NK", "DC", "Monocyte")
  
  x <- sort(unique(x))
  newX <- c()
  
  for (p in priority) {
    pIndices <- grepl(p, x)
    newX <- append(newX, x[pIndices])
    x <- x[!pIndices]
  }
  
  if (length(newX) > 0) {
    newX <- append(newX, x)
  }
  
  
  return(newX)
}


wrapNewAnnotLevel <- function(p) {
  return(p + plot_layout(tag_level = "new"))
}

unclipStripText <- function(g) {
  g2 <- patchworkGrob(g)
  for (i in which(grepl("strip-r", g2$layout$name))){
    
    g2$grobs[[i]]$grobs[[1]]$layout$clip <- "off"
  }
  
  return(g2)
}

processGregFlowData <- function(rdsFn) {
  dfLineageFilter <- readRDS(rdsFn)
  
  dfDiseaseScales <- dfLineageFilter %>%
    filter(!grepl("^CD. Mem CD.*$", metric)) %>%
    filter(!grepl("^CD. Mem .*\\+$", metric)) %>%
    filter(grepl("CD. ", metric)) %>%
    mutate(`Disease Status` = factor(`Disease Status`, levels = c("ND", "AAb+", "T1D"))) %>%
    mutate(LN_type = factor(LN_type, levels = c("pLN", "mLN", "Spleen"))) %>%
    mutate(cd = case_when(
      str_detect(metric, "CD4") ~ "CD4",
      str_detect(metric, "CD8") ~ "CD8")) %>%
    mutate(tpop = case_when(
      str_detect(metric, "Tn ") ~ "Tn",
      str_detect(metric, "Tnl") ~ "Tnl",
      str_detect(metric, "Tcm") ~ "Tcm",
      str_detect(metric, "Tem ") ~ "Tem",
      str_detect(metric, "Temra") ~ "Temra",
      str_detect(metric, "Tn$") ~ "Tn",
      str_detect(metric, "Tem$") ~ "Tem",
      str_detect(metric, "CD4 Mem") ~ "Mem")) %>%
    mutate(tpop = factor(tpop, levels = c("Tn", "Tnl", "Tcm", "Tem", "Temra", "Mem")))
  
  return(dfDiseaseScales)
}

saveFinalFigure <- function(
  plot,
  prefixDir = "../outs",
  fn,
  devices = c("png", "pdf", "pdf_base"),
  addTimestamp = FALSE,
  gheight,
  gwidth) {
  
  if (!is.vector(devices)) {
    devices <- c(devices)
  }
  

  for (d in devices) {
    
    if (addTimestamp) {
      currentTime <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
      if (d == "pdf_base") { 
        gfn <- glue("{prefixDir}/{d}/{fn}_{currentTime}.pdf")
      } else {
        gfn <- glue("{prefixDir}/{d}/{fn}_{currentTime}.{d}")
      }
      
    } else {
      if (d == "pdf_base") { 
        gfn <- glue("{prefixDir}/{d}/{fn}.pdf")
      } else {
        gfn <- glue("{prefixDir}/{d}/{fn}.{d}")
      }
    }
    
    
    if (d == "rds") {
      saveRDS(plot, gfn)
    } else if (d == "pdf") {
      ggsave(gfn, plot = plot, dpi = "retina", device = cairo_pdf, width = gwidth, height = gheight, units = "in")
    } else if (d == "pdf_base") {
      ggsave(gfn, plot = plot, dpi = "retina", device = "pdf", width = gwidth, height = gheight, units = "in")
    } else {
      ggsave(gfn, plot = plot, dpi = "retina", device = d, width = gwidth, height = gheight, units = "in")  
    }
  }
}
