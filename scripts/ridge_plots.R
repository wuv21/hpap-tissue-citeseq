library(magrittr)
library(colorspace)
library(scales)

goi = tibble::as_tibble(read.table("B_NK_CITEseq_gene_list_allComparisons.csv", sep = ',', header=TRUE, row.names=NULL))
goi$comparison = gsub("AAb", "AAb+", goi$comparison)

goi[goi$gene %in% c("CD27"),]$adt_number
so = readRDS("../seuMergedPostHSP_forFigures_2023-09-17_09-03-10.rds")
unique(so[["manualAnnot"]])

bcell_clusters = unique(so[["manualAnnot"]][,1][which(startsWith(so[["manualAnnot"]][,1], "B"))])
bcell_clusters
bcells = Cells(subset(so, manualAnnot %in% bcell_clusters))
so_bcells = subset(so, manualAnnot %in% bcell_clusters)
DefaultAssay(so_bcells) = "adt"
Idents(so_bcells) = "Disease_Status"

pdf("figures/CD27_ridges.pdf", width=8.5, height=17)
plts = list()
fcol = as.vector(sapply(c("red", "blue", "purple"), function(c) sapply(seq(1,8,1), function(x,c) colorspace::lighten(colorspace::desaturate(c, 0.10*x), 0.10*x), c)))
idxs = list(1:8,9:16,17:24)
fcol[1:5]
i=1
for (ds in unique(so[["Disease_Status"]][,1])) {
  plts[[i]] = Seurat::RidgePlot(subset(so_bcells, Disease_Status == ds), features="A0154", group.by = "DonorID", col = fcol[idxs[[i]]]) +
    ggplot2::ggtitle(sprintf("%s - %s", ds, "CD27")) +
    ggplot2::xlim(0,6.5)
  i=i+1
}
patchwork::wrap_plots(plts, ncol=1)
dev.off()


########################################
convert=goi[goi$gene %in% c("CD16", "CD56"),][,c("gene", "adt_number")]
convert=named_vector(convert$gene, convert$adt_number)
nk_cells = Cells(subset(so, manualAnnot %in% c("NK", "NK/ILC")))
so_nk = so[,nk_cells]
Idents(so_nk) = "Disease_Status"
DefaultAssay(so_nk) = "adt"

pdf("figures/CD56_CD16_NK_ridges.pdf", width=8.5, height=17)
for (mark in unique(names(convert))) {
  if (mark == "CD56") {
    intercept=2.3
  } else if (mark == "CD16") {
    intercept=1.2
  } else {
    stop("bad marker")
  }
  plts = list()
  fcol = as.vector(sapply(c("red", "blue", "purple"), function(c) sapply(seq(1,8,1), function(x,c) colorspace::lighten(colorspace::desaturate(c, 0.10*x), 0.10*x), c)))
  idxs = list(1:8,9:16,17:24)
  fcol[1:5]
  i=1
  for (ds in unique(so_nk[["Disease_Status"]][,1])) {
    plts[[i]] = Seurat::RidgePlot(subset(so_nk, Disease_Status == ds), features=convert[[mark]], group.by = "DonorID", col = fcol[idxs[[i]]]) +
      ggplot2::ggtitle(sprintf("%s - %s : %s", ds, mark, "NK cells")) +
      ggplot2::xlim(0,6.5)+
      ggplot2::geom_vline(xintercept=intercept, linetype="dashed")
    i=i+1
  }
  print(patchwork::wrap_plots(plts, ncol=1))
}
dev.off()

###############

convert=goi[goi$gene %in% c("NKG2D", "CD161"),][,c("gene", "adt_number")]
convert=named_vector(convert$gene, convert$adt_number)
nk_cells = Cells(subset(so, manualAnnot %in% c("NK", "NK/ILC")))
so_nk = so[,nk_cells]
Idents(so_nk) = "Disease_Status"
DefaultAssay(so_nk) = "adt"

pdf("figures/NKG2D_CD161_NK_ridges.pdf", width=8.5, height=17)
for (mark in unique(names(convert))) {
  if (mark == "NKG2D") {
    intercept=1.5
  } else if (mark == "CD161") {
    intercept=1.0
  } else {
    stop("bad marker")
  }
  plts = list()
  fcol = as.vector(sapply(c("red", "blue", "purple"), function(c) sapply(seq(1,8,1), function(x,c) colorspace::lighten(colorspace::desaturate(c, 0.10*x), 0.10*x), c)))
  idxs = list(1:8,9:16,17:24)
  fcol[1:5]
  i=1
  for (ds in unique(so_nk[["Disease_Status"]][,1])) {
    plts[[i]] = Seurat::RidgePlot(subset(so_nk, Disease_Status == ds), features=convert[[mark]], group.by = "DonorID", col = fcol[idxs[[i]]]) +
      ggplot2::ggtitle(sprintf("%s - %s : %s", ds, mark, "NK cells")) +
      ggplot2::xlim(0,4.5)+
      ggplot2::geom_vline(xintercept=intercept, linetype="dashed")
    i=i+1
  }
  print(patchwork::wrap_plots(plts, ncol=1))
}
dev.off()




