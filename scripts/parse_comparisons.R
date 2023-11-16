library(Seurat)

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
  pchd = ifelse(is.na(pvalues), "", "")
  pchd = ifelse(pvalues>=1.30103,"*",pchd)
  pchd = ifelse(pvalues>=10, "**", pchd)
  pchd = ifelse(pvalues>=50, "***", pchd)
  pchd
}

