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
  "tissue" = c("pLN" = "#0777B6", "mesLN" = "#73A942", "Spleen" = "#BC4749")
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

umapPlotThemeNoLeg <- list(
  # coord_fixed(),
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    plot.title.position = "panel",
    plot.title = element_text(size = 8, hjust = 0.5, face = "plain"))
)

# umapPlotThemeLeg <- theme(
#   legend.position = "bottom",
#   axis.ticks = element_blank(),
#   axis.text = element_blank(),
#   legend.margin = margin(0,0,0,0),
#   plot.title.position = "panel"
# )


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

saveFinalFigure <- function(
    plot,
  prefixDir = "../outs",
  fn,
  devices = c("png", "pdf"),
  gheight,
  gwidth) {
  
  if (!is.vector(devices)) {
    devices <- c(devices)
  }
  
  for (d in devices) {
    gfn <- glue("{prefixDir}/{d}/{fn}.{d}")
    
    if (d == "rds") {
      saveRDS(plot, gfn)
    } else if (d == "pdf") {
      ggsave(gfn, plot = plot, dpi = "retina", device = cairo_pdf, width = gwidth, height = gheight, units = "in")
    } else {
      ggsave(gfn, plot = plot, dpi = "retina", device = d, width = gwidth, height = gheight, units = "in")  
    }
  }
}