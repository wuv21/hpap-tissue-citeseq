library(Seurat)
library(magrittr)
library(dplyr)
library(ggplot2)
library(patchwork)
setwd("/home/ubuntu/mnt/betts/hpap/hpap-tissue-citeseq")
source("scripts/dimPlots.R")
source("scripts/parse_comparisons.R")
source("scripts/deg.R")
source("scripts/so_helpers.R")

set.seed(42)

tsa_catalog = readRDS("rds/tsa_catalog.rds")

so = readRDS("../seuMergedPostHSP_forFigures_2023-09-17_09-03-10.rds")

so_pln_only = readRDS("rds/so_pln_only.rds")
clusters = unique(so_pln_only[["manualAnnot"]])
bcell_clusters = clusters[which(startsWith(clusters[,1], "B")),1]
nk_clusters = c("NK", "NK/ILC")

convert = c(rep("All B cells combined", length(bcell_clusters)), rep("All NK Cells combined", length(nk_clusters)))
names(convert)=c(bcell_clusters, nk_clusters)

so_pln_only = soAddGroupedAnnotVar(so_pln_only, "manualAnnot", "groupedAnnot", convert)
unique(so_pln_only[["groupedAnnot"]])

goi = tibble::as_tibble(read.table("miscellaneous_gene_lists/HPAP_CITEseq_gene_list_V3.csv", sep = ',', header=TRUE, row.names=NULL, fill=NA))
goi$comparison = gsub("AAb", "AAb+", goi$comparison)
goi$modality
unique(goi[goi$modality == "RNA",]$Heatmap)
goi

#Bcells v3
annotVar = "groupedAnnot"
compsout_adt = findMarkersCombinatorial(subset(so_pln_only, groupedAnnot == "All B cells combined"), combVar = "Disease_Status", assay = "adt")
compsout_adt = compsout_adt[compsout_adt$gene %in% goi[goi$Heatmap == 1,]$adt_number,]
saveRDS(compsout_adt, "rds/wuv_compres_adt_genelist_V3.rds")

compsout_rna = findMarkersCombinatorial(subset(so_pln_only, groupedAnnot == "All B cells combined"), combVar = "Disease_Status", assay = "RNA")
compsout_rna = compsout_rna[compsout_rna$gene %in% goi[goi$Heatmap == 2,]$gene,]
saveRDS(compsout_rna, "rds/wuv_compres_allbcells_rna_genelist_V3.rds")

#Bcells v3 separate clusters
annotVar = "manualAnnot"
compsout_rna = lapply(unique(bcell_clusters), function(xx) {
                        df = findMarkersCombinatorial(subset(so_pln_only, !!sym(annotVar) == xx), combVar = "Disease_Status", assay = "RNA")
                        attr(df, "cluster") = xx
                        print(df)
})
names(compsout_rna) = bcell_clusters
names(compsout_rna)

saveRDS(compsout_rna, "rds/wuv_compres_bcellsep_rna_genelist_V3.rds")

#NK cells v1, for ranked plots
so_pln_only = readRDS("rds/so_pln_only.rds")
colnames(so_pln_only[[]])
unique(so_pln_only[["groupedAnnot"]])
goiv1 = tibble::as_tibble(read.table("B_NK_CITEseq_gene_list.csv", sep = ',', header=TRUE, row.names=NULL, fill=NA))
goiv1$comparison = gsub("AAb", "AAb+", goiv1$comparison)
goiv1

compsout_nk_rna = findMarkersCombinatorial(subset(so_pln_only, groupedAnnot == "All NK Cells combined"), combVar = "Disease_Status", assay = "RNA")
compsout_nk_rna
saveRDS(compsout_nk_rna, "rds/wuv_compres_rna_nkcells.rds")
