library(Seurat)

soAddGroupedAnnotVar = function(so, annotVar, groupedAnnotVar, convert) {
  so[[groupedAnnotVar]] = unname(sapply(so[[annotVar]][,1], function(old) convert[old]))
  so
}

seuratObjMetaTibble = function(so, assay = NULL, slot = "data", barcodeVar = "cell") {
  if (!is.null(assay)) {
    DefaultAssay(so) = assay
  }
  ret = as_tibble(t(as.matrix(GetAssayData(so))), rownames = "cell") %>% left_join(as_tibble(so[[]], rownames="cell"))
  attr(ret, "datacol") = seq_along(rownames(so))+1
  ret
}
