library(Seurat)


soAddGroupedAnnotVar = function(so, annotVar, groupedAnnotVar, convert) {
  so[[groupedAnnotVar]] = unname(sapply(so[[annotVar]][,1], function(old) convert[old]))
  so
}

split_column_list = function(cc) { 
  data.frame(id1 = sapply(cc, function(x) x[1]), id2 = sapply(cc, function(x) x[2]))
}
expand_comparisons_tidy = function(genelist) {
  genelist %>% 
    dplyr::mutate(to_split = comparison) %>%
    dplyr::mutate(to_split = lapply(stringr::str_split(to_split, '_'), combn, 2)) %>%
    dplyr::mutate(to_split = lapply(to_split, apply, 2, as.list)) %>% 
    tidyr::unnest(to_split) %>% 
    dplyr::mutate(to_split = lapply(to_split, unlist)) %>%
    dplyr::mutate(split_column_list(to_split)) %>%
    dplyr::select(-to_split) %>%
    tidyr::unite(col = comparison, id1, id2, remove=FALSE)
}

expand_comparisons_base = function(genelist) {
  combos<-lapply(strsplit(genelist$comparison,'_'),combn,2)
  orig_row_range = 1:nrow(genelist)
  ncombos = sapply(combos,ncol)
  expanded_row_indices=rep(orig_row_range,ncombos)
  out<-cbind(genelist[expanded_row_indices,],'id'=t(do.call(cbind,combos)))
  out
}

expand_comparisons_base_short = function(genelist, clusters) {
  
  # Expand rows by semicolon-separated comparison column
  expand_combinations_by_row=by(genelist,1:nrow(genelist),function(yy)cbind(yy,'id'=t(combn(strsplit(yy$comparison,'_')[[1]],2))))
  out<-do.call(rbind,expand_combinations_by_row)
  out$comparison = paste(out$id.1, out$id.2, sep = "_")
  # Expand rows by clusters to compare
  expand_comparisons_by_cluster = by(out,seq_along(rownames(out)),function(cc) {
                                       rm_cluster_col = which(colnames(cc)=="cluster")
                                       print(cc)
                                       print(clusters[clusters$cgroup == cc$cluster,])
                                       cbind(cc,clusters[clusters$cgroup == cc$cluster,])[,-rm_cluster_col]
                                      }
  )
  print(expand_comparisons_by_cluster)
  do.call(rbind, expand_comparisons_by_cluster)
}

new_comparison = function(markers = c(), id1, id2, modality, cluster, results = NULL) {
  structure(markers, id1 = id1, id2 = id2, modality = modality, cluster = cluster, class = "comparison")
}

Idents.comparison = function(comparison) {
  list(id1 = attr(comparison, "id1"), id2 = attr(comparison, "id2"))
}

comparisons_from_genelist = function(genelist) {
  genelist$modal_compar_clustid = paste(genelist$modality, genelist$comparison,  genelist$cidx, sep = "__")
  by(genelist, genelist$modal_compar_clustid, function(gg)new_comparison(gg$feature, gg$id.1[1], gg$id.2[1], modality=gg$modality[1], cluster=gg$cluster[1]))
}

named_vector = function(keys, values) {
  names(values) = keys
  values
}

convert_name = function(name, convert) {
  convert[[name]]
}

fill_pvalues_from_wuv = function(scaledExp, wuv_compres, modality, goi=NULL) {
  wuv_compres = cbind(wuv_compres_rna, do.call(rbind, strsplit(gsub("[_]vs[_]", "-",wuv_compres_rna$matchup), "-")))
  features = rownames(scaledExp)
  print(features)
  pltcomps = matrix(nrow = length(features), ncol = 3)
  print(pltcomps)
  rownames(pltcomps) = features
  colnames(pltcomps) = c("ND-AAb+", "ND-T1D", "AAb+-T1D")
  for(comp in strsplit(colnames(pltcomps), '-')) {
    rele_res = wuv_compres[which(wuv_compres$matchup == sprintf("%s_vs_%s", comp[1], comp[2])),]
    rownames(rele_res) = rele_res$gene
    if (length(seq_along(rele_res[,1])) == 0) {
      print(sprintf("No rele_res for comparison betwixt: %s - %s", comp[1], comp[2]))
      next
    }
    if (!is.null(goi) && modality == "adt") {
      rownames(rele_res) = sapply(rownames(rele_res), convert_name, named_vector(goi$adt_number, goi$gene))
    }
    rows = which(rownames(pltcomps) %in% rownames(rele_res))
    print(rele_res)
    print(rownames(pltcomps)[rows])
    pvalues = rele_res[rownames(pltcomps)[rows],]
    print(pvalues)
    if (length(seq_along(rele_res[,1])) > 0) {
      pltcomps[,sprintf("%s-%s", comp[1], comp[2])] = pvalues[rownames(pltcomps),]$p_val_adj_all
    }
    #pltcomps[is.na(pltcomps)] = 1
  }
  pltcomps
}

fill_pvalues = function(scaledExp, modality, goi=NULL) {
  features = rownames(scaledExp)
  pltcomps = matrix(nrow = length(features), ncol = 3)
  rownames(pltcomps) = features
  colnames(pltcomps) = c("ND-AAb+", "ND-T1D", "AAb+-T1D")
  pltcomps
  compres[[1]]
  for(comp in strsplit(colnames(pltcomps), '-')) {
    rele_res = compres[unname(which(sapply(compres, function(cc) all(comp %in% Idents(cc)) & attr(cc, "cluster") == attr(scaledExp, "cluster") & attr(cc, "modality") == modality)))]
    if (length(rele_res) == 0) {
      print(sprintf("No rele_res for comparison betwixt: %s - %s", comp[1], comp[2]))
      next
    }
    stopifnot(length(rele_res) <= 1)
    rele_res = attr(rele_res[[1]], "results")
    if (!is.null(goi)) {
      rownames(rele_res) = sapply(rownames(rele_res), convert_name, named_vector(goi$adt_number, goi$gene))
    }
    rows = which(rownames(pltcomps) %in% rownames(rele_res))
    print(rele_res)
    print(rownames(pltcomps)[rows])
    pvalues = rele_res[rownames(pltcomps)[rows],]
    if (length(seq_along(rele_res[,1])) > 0) {
      pltcomps[,sprintf("%s-%s", comp[1], comp[2])] = pvalues[rownames(pltcomps),]$p_val
    }
    #pltcomps[is.na(pltcomps)] = 1
  }
  pltcomps
}

get_genes = function(compres, so, modality, target_clusters) {
  genes = compres[unname(which(sapply(compres, function(cc)attr(cc, "modality") == modality & attr(cc, "cluster") %in% target_clusters)))]
  print(genes)
  genes = unique(unname(unlist(lapply(genes, function(cc) as.character(cc)))))
  genes = genes[which(genes %in% rownames(so))]
  genes
}

pval2pch = function(pvalues) {
  pchd = ifelse(is.na(pvalues), "", "*")
  pchd = ifelse(pvalues>=10, "**", pchd)
  pchd = ifelse(pvalues>=50, "***", pchd)
  pchd
}

make_hm = function(scaled_data,
                   pvalues,
                   colmin, colmax,
                   rect_label,
                   colscale,
                   annotVar,
                   pctexp = NULL,
                   rannot_labels = TRUE,
                   showlegends=FALSE) {
  column_order = unname(sapply(colnames(scaled_data), function(cn) which(grepl(strsplit(cn, " ")[[1]][1], c("ND", "AAb+", "T1D")))))
  if (!is.null(pctexp)) {
    ct =gsub("[(][^)]+[)][ ]", "", colnames(scaled_data))
    bars = pctexp[pctexp[annotVar] == ct[1],]
    bars_mat = as.matrix(bars[,-c(1,2)])
    hma_colscale = list("percent_expressing" = circlize::colorRamp2(breaks=c(0,1), rev(hcl.colors(n=24, palette = "Purples"))[c(1,24)]),
                        "pvalues" = circlize::colorRamp2(breaks=c(0,10,50,100,323), hcl.colors(n=24, palette = "Earth")[c(6,12,16,20,24)]))
    print(colscale)
#rannot = ComplexHeatmap::HeatmapAnnotation(which="row", bar = ComplexHeatmap::anno_barplot(bars_mat, beside = TRUE))
    rannot = ComplexHeatmap::rowAnnotation(percent_expressing = ComplexHeatmap::anno_simple(as.matrix(bars_mat),
                                                                                 col=hma_colscale[["percent_expressing"]], 
                                                                                 na_col="white",
                                                                                 simple_anno_size=unit(0.2, "in"), 
                                                                                 gp=gpar(col="white")),
                                           pvalues = ComplexHeatmap::anno_simple(as.matrix(pvalues),
                                                                                 col=hma_colscale[["pvalues"]], 
                                                                                 na_col="white",
                                                                                 pch=pval2pch(pvalues),
                                                                                 pt_size=unit(0.4, "in"),
                                                                                 simple_anno_size=unit(0.4, "in"), 
                                                                                 gp=gpar(col="white")),
                                           col=hma_colscale,
                                           show_annotation_name = rannot_labels, 
                                           simple_anno_size= unit(0.2, "in"), 
                                           gp=gpar(col="white"),
                                           annotation_name_side=c("bottom","bottom"))
  } else {
    rannot = NULL
  }
  hm = ComplexHeatmap::Heatmap(
                          scaled_data,
                          na_col = "grey90", 
                          cluster_rows=T, 
                          cluster_columns=F, 
                          rect_gp = gpar(col = "white", lwd = 2), 
                          col = colscale,
                          row_title = rect_label, 
                          row_title_gp = gpar(fill = "black", col = "white", border = "black"),
                          show_heatmap_legend = showlegends,
                          heatmap_legend_param = list(title="Norm Expr"),
                          border = FALSE,
                          column_order = column_order,
                          right_annotation=rannot
  )
  hm
}
