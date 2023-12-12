library(Seurat)

percent_expressing = function(sotib, compareVar, annotVar, zero = 0, barcodeVar = "cell") {
  sotib = sotib[,c(barcodeVar, colnames(sotib[,attr(sotib, "datacol")]), compareVar, annotVar)]
  sotib %>%
    select(-barcodeVar) %>%
    group_by(!!sym(compareVar), !!sym(annotVar)) %>%
    summarize_all(function(x) length(which(x > zero))/length(x)) %>%
    pivot_longer(cols=-c(compareVar, annotVar), names_to = "feature", values_to = "pctexp") %>%
    pivot_wider(names_from = compareVar, values_from = "pctexp") %>%
    select(annotVar, feature, unique(sotib[[compareVar]]))
}
expression_stderr_bars = function(sotib, compareVar, annotVar, barcodeVar = "cell", method = mean, bar_method = function(x) sd(x)/sqrt(n())) {
  sotib = sotib[,c(colnames(sotib[,attr(sotib, "datacol")]), compareVar, annotVar)]
  sotib %>%
    pivot_longer(cols = !c(compareVar, annotVar), names_to = "Gene", values_to = "normalized_expression") %>% 
    group_by(!!sym(compareVar), !!sym(annotVar), Gene) %>%
    summarize_all(list(mid = method, barinc = bar_method)) %>%
    mutate(ymax = mid + barinc) %>%
    mutate(ymin = mid - barinc)
}

mean_expression = function(sotib, compareVar, annotVar=NULL, barcodeVar = "cell") {
  if (!is.null(annotVar)) {
    sotib = sotib[,c(colnames(sotib[,attr(sotib, "datacol")]), compareVar, annotVar)]
    sotib %>%
      group_by(!!sym(compareVar), !!sym(annotVar)) %>%
      summarize_all(mean)
  } else {
    sotib = sotib[,c(colnames(sotib[,attr(sotib, "datacol")]), compareVar)]
    sotib %>%
      group_by(!!sym(compareVar)) %>%
      summarize_all(mean)
  }
}

scale_expression = function(sotib, compareVar, by, barcodeVar = "cell") {
  lapply(unique(sotib[[by]]), function(xx) {
           out = sotib[which(sotib[[by]] == xx),]
           rnames = sprintf("(%s) %s", out[[compareVar]], out[[by]])
           out = out[,-c(1,2)]
           out = t(scale(as.matrix(out)))
           colnames(out) = rnames

           #Check and see if this is necessary
           attr(out, "cluster") = xx

           out
  })
}

fillPvalsFromFindFeaturesCombinatorial = function(scaledExp, wuv_compres, modality, pvalue_column, na_fill_value, goi=NULL) {
  wuv_compres = wuv_compres[[attr(scaledExp, "cluster")]]
  wuv_compres = cbind(wuv_compres, do.call(rbind, strsplit(gsub("[_]vs[_]", "-",wuv_compres$matchup), "-")))
  features = rownames(scaledExp)
  pltcomps = matrix(nrow = length(features), ncol = 3)
  pltcomps[,] = na_fill_value
  rownames(pltcomps) = features
  colnames(pltcomps) = c("ND-AAb+", "ND-T1D", "AAb+-T1D")
  for(comp in strsplit(colnames(pltcomps), '-')) {
    wuv_matchup =  sprintf("%s_vs_%s", comp[1], comp[2])
    rele_res = wuv_compres[which(wuv_compres$matchup == wuv_matchup),]
    if (length(seq_along(rele_res[,1])) == 0) {
      wuv_matchup_swapped = sprintf("%s_vs_%s", comp[2], comp[1])
      rele_res = wuv_compres[which(wuv_compres$matchup == wuv_matchup_swapped),]
      if (length(seq_along(rele_res[,1])) == 0) {
        print(sprintf("No rele_res for comparison betwixt: %s - %s", comp[1], comp[2]))
        next
      } else {
        wuv_matchup = wuv_matchup_swapped
      }
    }
    rownames(rele_res) = rele_res$gene
    if (!is.null(goi) && modality == "adt") {
      rownames(rele_res) = sapply(rownames(rele_res), convert_name, named_vector(goi$adt_number, goi$gene))
    }
    rows = which(rownames(pltcomps) %in% rownames(rele_res))
    pvalues = rele_res[rownames(pltcomps)[rows],]
    if (length(seq_along(rele_res[,1])) > 0) {
      pltcomps[,sprintf("%s-%s", comp[1], comp[2])] = pvalues[rownames(pltcomps),][[pvalue_column]]
    }
  }
  apply(pltcomps, 2, function(x) ifelse(is.na(x), na_fill_value, x))
}

pvalues2pubr = function(pvalues, gene, xlvl, bracketx) {
  offset = 0
  out = list()
  for (x in xlvl) {
    recast = as_tibble(do.call(rbind, lapply(names(pvalues[[x]][g,]), function(n) c(strsplit(n, "[-]")[[1]], pvalues[[x]][g,n]))))
    colnames(recast) = c("group1", "group2", "p.signif")
    recast$group1x = sapply(recast$group1, function(x) bracketx(offset, x))
    recast$group2x = sapply(recast$group2, function(x) bracketx(offset, x))
    recast$p.signif = as.numeric(recast$p.signif)
    recast$p.sym = pValSymnum(recast$p.signif)
    recast$p.sym = ifelse(is.na(recast$p.sym), "ns", recast$p.sym)
    offset=offset+1
    out[[offset]] = recast
  }
  return(do.call(rbind, out))
}
