source("scripts/generic.R")
set.seed(42)

findMarkersCombinatorial <- function(
  seuratObj,
  combVar,
  assay = NULL,
  logfc.threshold = 0.1,
  pAdjustMethod = "bonferroni",
  ...
) {
  if (!is.null(assay) && DefaultAssay(seuratObj) != assay) {
    message(sprintf("Switching Assay to: %s", assay)) 
    DefaultAssay(seuratObj) = assay
  }

  Idents(seuratObj) <- combVar
  combs <- combn(unique(seuratObj@meta.data[, combVar]), 2)
  
  comparisons <- lapply(c(1:ncol(combs)), function(i) {
    group1 <- as.character(combs[1, i])
    group2 <- as.character(combs[2, i])
    
    comp <- paste0(group1, "_vs_", group2)
    message(paste0("testing ", comp))
    
    deg <- FindMarkers(seuratObj,
      ident.1 = group1,
      ident.2 = group2,
      logfc.threshold = logfc.threshold,
      ...)
    
    deg <- deg %>%
      mutate(
        gene = rownames(.),
        upregulated = as.character(ifelse(avg_log2FC > 0, group1, group2)),
        matchup = comp
      )
    
    rownames(deg) <- NULL
    
    return(deg)
  }) 

  #Remove any row-less elements

  
  message("collating and correcting for multiple tests...")
  comparisons <- bind_rows(comparisons) %>%
    mutate(p_val_adj_all = p.adjust(p_val, method = pAdjustMethod))
  
  return(comparisons)
}

getPrestoResult <- function(
  seu,
  posPop,
  negPop,
  group
) {
  res <- wilcoxauc(
    X = seu,
    groups_use = c(negPop, posPop),
    group_by = group)
  
  return(res)
}

cleanPrestoResult <- function(
  df,
  posPop,
  negPop,
  minPctIn = 60,
  fromPosPopPerspectiveOnly = TRUE
) {
  cleanDf <- df %>%
    filter(pct_in >= minPctIn) %>%
    mutate(log2FC = log2(10 ** logFC)) %>%
    select(-logFC, -statistic, -pval) %>%
    mutate(padj_trans = -1 * log10(padj)) %>%
    mutate(piScore = padj_trans * log2FC) %>%
    mutate(interestGroup = case_when(
      padj < 0.05 & log2FC > 0 ~ paste0("Upregulated in ", posPop),
      padj < 0.05 & log2FC < 0 ~ paste0("Upregulated in ", negPop),
      TRUE ~ "n.s."
    )) %>%
    group_by(group) %>%
    arrange(desc(abs(piScore)))
  
  
  if (fromPosPopPerspectiveOnly) {
    cleanDf <- cleanDf %>%
      filter(group == posPop)
  }
  
  return(cleanDf)
}

getPrestoResultForAllGroupComb <- function(
  seu,
  pairedGroups,
  groupInterest
) {
  pairedGroupNames <- sapply(pairedGroups, function(x) {
    return(paste(x, collapse = " vs "))
  })
  
  resAll <- lapply(pairedGroups, function(x) {
    message(paste0("specifically comparing ", paste(x, collapse = " vs ")))
    
    res <- getPrestoResult(
      seu = seu,
      posPop = x[2],
      negPop = x[1],
      group = groupInterest
    )
    
    res <- cleanPrestoResult(
      df = res,
      posPop = x[2],
      negPop = x[1]
    )
    
    return(res)
  })
  
  names(resAll) <- pairedGroupNames
  return(resAll)
}


plotCleanMA <- function(
  cleanedPrestoDf,
  fn,
  colorScale = NULL,
  labelTopN = 15,
  height = 5,
  width = 5,
  devices = "png",
  title = NULL,
  returnPlot = FALSE
) {
  
  if (is.null(colorScale)) {
    colorScale <- DISEASESTATUSCOLORS
    names(colorScale) <- paste0("Upregulated in ", names(colorScale))
    colorScale <- c("n.s" = "#CCCCCC", colorScale)
  }
  
  labelDf <- cleanedPrestoDf %>%
    filter(padj < 0.05) %>%
    arrange(desc(log2FC))
  
  labelDf <- labelDf[c(1:labelTopN, (nrow(labelDf) - labelTopN):nrow(labelDf)), ]
  
  p <- cleanedPrestoDf %>%
    arrange(interestGroup) %>%
    ggplot(aes(x = avgExpr, y = log2FC, color = interestGroup)) +
    geom_point() +
    theme_classic() +
    scale_color_manual(values = colorScale) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.key.size = unit(BASEPTFONTSIZE * 1.1, 'points'),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.margin = margin(t = -5),
      plot.title = element_text(hjust = 0.5, size = BASEPTFONTSIZE + 2),
      plot.title.position = "panel",
      axis.text = element_text(size = BASEPTFONTSIZE),
      axis.title = element_text(size = BASEPTFONTSIZE),
      legend.text = element_text(size = BASEPTFONTSIZE),
      panel.background = element_rect(fill = "transparent", colour = "#000000", size = 1),
      plot.background = element_rect(fill = "transparent", colour = NA),
      axis.line = element_blank()
    ) +
    ggrepel::geom_text_repel(
      data = labelDf,
      aes(label = feature),
      max.overlaps = labelTopN,
      force = labelTopN,
      size = (BASEPTFONTSIZE) / ggplot2:::.pt,
      min.segment.length = 0,
      color = "black") +
    labs(x = "Mean expression", y = "log2FC")
  
  if (!is.null(title)) {
    p <- p + labs(title = title)
  }
  
  savePlot(p, fn = fn, gheight = height, gwidth = width, devices = devices)
  
  if (returnPlot) {
    return(p)
  }
}



plotCombinatorialDEGLollipop <- function(combMarkersRes, title) {
  df <- combMarkersRes %>%
    group_by(matchup, upregulated) %>%
    # mutate(piscore = -1 * log10(p_val_adj_all) * avg_log2FC) %>%
    filter(p_val_adj_all < 0.05) %>%
    slice_max(abs(avg_log2FC), n = 20) %>%
    group_by(matchup) %>%
    arrange(desc(avg_log2FC), .by_group = TRUE)
  
  maxCounts <- df %>%
    summarize(groupCounts = n())
  maxCounts <- ceiling(log10(max(maxCounts$groupCounts))) + 1
  
  df <- df %>%
    mutate(facet_gene_number = cur_group_id() * (10 ** maxCounts) + row_number()) %>%
    mutate(facet_gene_number = paste0(facet_gene_number, "_", gene)) %>%
    mutate(facet_gene_number = factor(facet_gene_number, levels = stringr::str_sort(facet_gene_number, numeric = TRUE)))
  
  print(df %>%
      select(avg_log2FC, facet_gene_number, gene, matchup) %>%
      slice_max(abs(avg_log2FC), n = 5) %>%
      arrange(avg_log2FC, .by_group = TRUE))
  
  p <- df %>%
    {
      ggplot(., aes(y = facet_gene_number, x = avg_log2FC, color = upregulated)) +
        geom_point(size = 1) +
        geom_vline(xintercept = 0) +
        scale_color_manual(values = COLORS[["disease"]]) +
        scale_y_discrete(labels = function(x) gsub("\\d+_", "", x)) +
        facet_wrap(~ matchup, nrow = 1, scales = "free") +
        theme_bw() +
        theme(
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          strip.background = element_blank()) +
        labs(
          color = "Disease",
          y = "Gene",
          x = "Average log2 Fold Change",
          title = title
        )
    }
  
  return(p)
}


volcanoPlotTheme <- theme(
  legend.position = "bottom",
  legend.title = element_blank(),
  legend.key.size = unit(BASEPTFONTSIZE * 1.1, 'points'),
  legend.background = element_rect(fill = "transparent", colour = NA),
  legend.margin = margin(t = -5),
  plot.title = element_text(hjust = 0.5, size = BASEPTFONTSIZE + 2),
  plot.title.position = "panel",
  axis.text = element_text(size = BASEPTFONTSIZE - 2),
  axis.title = element_text(size = BASEPTFONTSIZE - 2),
  legend.text = element_text(size = BASEPTFONTSIZE - 2),
  panel.background = element_rect(fill = "transparent", colour = NA),
  plot.background = element_rect(fill = "transparent", colour = NA)
)


plotDegVolcano <- function(
  df,
  posPopulation,
  xvar = "log2FC",
  yvar = "padj_trans",
  colorvar = "interestGroup",
  title = "",
  xlabel = xvar,
  ylabel = yvar,
  topNLabel = 15,
  customTheme = NULL,
  rasterise = FALSE,
  manualColorPalette = NULL,
  fn = NULL,
  returnPlot = is.null(fn)
) {
  
  dfOnlyPos <- df %>%
    filter(group == posPopulation)
  
  labelDf <- dfOnlyPos %>%
    arrange(piScore)
  
  labelDf <- labelDf[c(1:topNLabel, (nrow(labelDf) - topNLabel):nrow(labelDf)), ]
  
  p <- ggplot(dfOnlyPos, aes(x = .data[[xvar]], y = .data[[yvar]], color = .data[[colorvar]]))
  
  if (rasterise) {
    p <- p +
      ggrastr::rasterise(geom_point(size = 0.5, alpha = 0.5), dpi = 300)
  } else {
    p <- p +
      geom_point(size = 0.5, alpha = 0.5)
  }
  
  p <- p +
    ggrepel::geom_text_repel(
      data = labelDf,
      aes(label = feature),
      max.overlaps = topNLabel,
      force = 15,
      size = (BASEPTFONTSIZE - 2) / ggplot2:::.pt,
      min.segment.length = 0,
      color = "black") +
    geom_hline(yintercept = -1 * log10(0.05), linetype = "dotted") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(
      title = title,
      x = xlabel,
      y = ylabel) +
    theme_classic() +
    volcanoPlotTheme
  
  if (!is.null(manualColorPalette)) {
    p <- p +
      scale_color_manual(values = manualColorPalette)
  } else {
    p <- p +
      scale_color_manual(values = c("grey", "blue", "red"))
  }
  
  if (returnPlot) {
    return(p)
  }
}


plotGseaEnrich <- function(
  df,
  groupOfInterest,
  feature = "auc",
  msigdbCategory = "C2",
  msigdbSubcategory = "CP:REACTOME",
  removePrefix = "REACTOME ",
  showTopN = 20,
  returnDf = FALSE,
  barColor = NULL,
  fn = NULL,
  returnPlot = is.null(fn)
) {
  
  mSigs <- msigdbr(species = "Homo sapiens", category = msigdbCategory, subcategory = msigdbSubcategory)
  fgseaSets <- mSigs %>% split(x = .$gene_symbol, f = .$gs_name)
  
  df <- df %>%
    ungroup() %>%
    filter(group == groupOfInterest) %>%
    arrange(desc(auc)) %>%
    dplyr::select(feature, auc)
  
  ranks <- tibble::deframe(df)
  
  message(glue("Computing pathway enrichment for {msigdbCategory} with sub-category {msigdbSubcategory}"))
  fgseaRes <- fgsea(fgseaSets, stats = ranks, scoreType = "pos")
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
  
  if (returnDf) {
    return(fgseaResTidy)
  }
  
  fgseaResTidy <- fgseaResTidy %>%
    filter(padj < 0.05)

  
  if (nrow(fgseaResTidy) == 0) {
    message("No significant pathways found.")
    return(NULL)
  }
  
  fgseaResTidy <- fgseaResTidy %>%
    slice_max(NES, n = showTopN) %>%
    mutate(pathway = gsub("_", " ", pathway)) %>%
    mutate(pathway = gsub(paste0("^", removePrefix), "", pathway)) %>%
    arrange(NES) %>%
    mutate(pathway = factor(pathway, levels = pathway))
  
  p <- ggplot(fgseaResTidy, aes(y = pathway, x = NES)) +
    geom_col(fill = ifelse(is.null(barColor), "#444444", barColor)) +
    scale_y_discrete(labels = function(x) {return(str_wrap(tolower(x), width = 60))}) +
    scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    theme_classic() +
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = BASEPTFONTSIZE - 2),
      panel.background = element_blank(),
      axis.text = element_text(size = BASEPTFONTSIZE - 2),
      axis.line = element_line(colour = "#000000"),
      panel.grid = element_blank()
    )
  
  if (returnPlot) {
    return(p)
  }
}

plotVolcanoAndGsea <- function(
  df,
  posPopulation,
  negPopulation,
  category = "C2",
  subcategory = "CP:REACTOME",
  gseaTopN = 10,
  height = 3.75,
  width = 7,
  devices = "png",
  fn = NULL,
  returnPlot = is.null(fn)
) {
  
  negPlotGseaEnrich <- plotGseaEnrich(
    df,
    groupOfInterest = negPopulation,
    showTopN = gseaTopN,
    barColor = "#0000ff",
    msigdbCategory = category,
    msigdbSubcategory = subcategory)
  
  posPlotGseaEnrich <- plotGseaEnrich(
    df,
    groupOfInterest = posPopulation,
    showTopN = gseaTopN,
    barColor = "#ff0000",
    msigdbCategory = category,
    msigdbSubcategory = subcategory)
  
  volcanoPlot <- plotDegVolcano(df, posPopulation = posPopulation)
  
  if (is.null(negPlotGseaEnrich)) {
    negPlotGseaEnrich <- plot_spacer()
    nRowsForNegPlot <- 0
    
  } else {
    nRowsForNegPlot <- nrow(negPlotGseaEnrich$data)
  }
  
  if (is.null(posPlotGseaEnrich)) {
    posPlotGseaEnrich <- plot_spacer()
    nRowsForPosPlot <- 0
    
  } else {
    nRowsForPosPlot <- nrow(posPlotGseaEnrich$data)
  }
  
  maxRows <- max(nRowsForNegPlot, nRowsForPosPlot)

  patchP <- volcanoPlot + 
    (negPlotGseaEnrich / 
    posPlotGseaEnrich /
    plot_spacer() +
    plot_layout(heights = c(nRowsForNegPlot, nRowsForPosPlot, gseaTopN * 2 - nRowsForPosPlot - nRowsForNegPlot)))

  
  patchP <- patchP +
    plot_layout(widths = c(7, 3))
  
  if (returnPlot) {
    return(patchP)
  }
  
  
  if (!is.null(fn)) {
    savePlot(
      plot = patchP,
      fn = fn,
      devices = devices,
      gheight = height,
      gwidth = width
    )
  }
}
