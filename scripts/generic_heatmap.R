# This function needs to be further generalized
# Functionality has to be moved from collar and into function (eg legends)
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
