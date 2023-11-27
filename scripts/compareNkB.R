library(Seurat)
library(magrittr)
library(dplyr)
library(ggplot2)
library(patchwork)
setwd("/home/ubuntu/mnt/betts/hpap/hpap-tissue-citeseq")
source("scripts/dimPlots.R")
source("scripts/parse_comparisons.R")
source("scripts/deg.R")

tsa_catalog = readRDS("rds/tsa_catalog.rds")

so = readRDS("../seuMergedPostHSP_forFigures_2023-09-17_09-03-10.rds")

so_pln_only = readRDS("rds/so_pln_only.rds")
clusters = unique(so_pln_only[["manualAnnot"]])
bcell_clusters = clusters[which(startsWith(clusters[,1], "B")),1]
nk_clusters = c("NK", "NK/ILC")

convert = c(rep("All B cells combined", length(bcell_clusters)), rep("All NK Cells combined", length(nk_clusters)))
names(convert)=c(bcell_clusters, nk_clusters)

# Make another column in so metadata for combined cluster groups
# so_pln_only[["groupedAnnot"]] = ifelse(so_pln_only[["manualAnnot"]][,1] %in% bcell_clusters, "All B Cells combined", NA)
#so_pln_only[["groupedAnnot"]] = ifelse(so_pln_only[["manualAnnot"]][,1] %in% nk_clusters, "All NK Cells combined", so_pln_only[["groupedAnnot"]][,1])
#unique(so_pln_only[["groupedAnnot"]])

so_pln_only = soAddGroupedAnnotVar(so_pln_only, "manualAnnot", "groupedAnnot", convert)
unique(so_pln_only[["groupedAnnot"]])

goi = tibble::as_tibble(read.table("HPAP_CITEseq_gene_list_V3.csv", sep = ',', header=TRUE, row.names=NULL, fill=NA))
goi$comparison = gsub("AAb", "AAb+", goi$comparison)
goi$modality
unique(goi[goi$modality == "RNA",]$Heatmap)
goi

annotVar = "manualAnnot"

#Bcells v3
compsout_adt = findMarkersCombinatorial(subset(so_pln_only, groupedAnnot == "All B Cells combined"), combVar = "Disease_Status", assay = "adt")
compsout_adt = compsout_adt[compsout_adt$gene %in% goi[goi$Heatmap == 1,]$adt_number,]
saveRDS(compsout_adt, "rds/wuv_compres_adt_genelist_V3.rds")

compsout_rna = findMarkersCombinatorial(subset(so_pln_only, groupedAnnot == "All B Cells combined"), combVar = "Disease_Status", assay = "RNA")
compsout_rna = compsout_rna[compsout_rna$gene %in% goi[goi$Heatmap == 2,]$gene,]
saveRDS(compsout_rna, "rds/wuv_compres_rna_genelist_V3.rds")

#Bcells v3 separate clusters
compsout_rna = lapply(unique(bcell_clusters), function(xx) {
                        df = findMarkersCombinatorial(subset(so_pln_only, !!sym(annotVar) == xx), combVar = "Disease_Status", assay = "RNA")
                        df["cluster"] = xx
})
compsout_rna
compsout_rna = lapply(compsout_rna, function(compsout_rna) compsout_rna[compsout_rna$gene %in% goi[goi$modality == "RNA",]$gene,])
saveRDS(compsout_rna, "rds/wuv_compres_rna_bcellsep_genelist_V3.rds")
readRDS( "rds/wuv_compres_rna_bcellsep_genelist_V3.rds")

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


#################### other abstractions

# Generate the comparisons list object for comparing all of the individual clusters
clusters = list()
clusters[["all B cell clusters"]] = data.frame(cgroup = "all B cell clusters", cluster = bcell_clusters)
clusters[["NK, NK/ILC"]] = data.frame(cgroup = "NK, NK/ILC", cluster = nk_clusters)
so = readRDS("../seuMergedPostHSP_forFigures_2023-09-17_09-03-10.rds")
clusters = do.call(rbind, clusters)
clusters$cidx =seq_along(rownames(clusters))
rownames(clusters) = seq_along(rownames(clusters))
clusters

exp_goi=expand_comparisons_base_short(goi, clusters)

genelist = exp_goi %>%
  mutate(feature = ifelse(modality == "ADT", adt_number, gene)) %>%
  mutate(modality = ifelse(modality == "ADT", "adt", modality))

comparisons = comparisons_from_genelist(genelist)

Seurat::Idents(so) = "Disease_Status"
so_pln_only = subset(so, Tissue %in% c("pLN-H", "pLN-T"))
comparisons_with_results = lapply(comparisons, function(cc, so) {
        print(sprintf("%s : %s, %s - %s", attr(cc, "modality"), attr(cc, "id1"), attr(cc, "id2"), attr(cc, "cluster")))
        DefaultAssay(so) = attr(cc, "modality")
        so_sub = subset(so, manualAnnot == attr(cc, "cluster"))
        if (length(which(Idents(so_sub) == attr(cc, "id1"))) == 0 | length(which(Idents(so_sub) == attr(cc, "id2"))) == 0) {
          attr(cc, "results") = NULL
        } else {
          features = as.character(cc)[which(cc %in% rownames(so_sub))]
          attr(cc, "results") = Seurat::FindMarkers(so_sub, ident.1 = attr(cc, "id1"), ident.2 = attr(cc, "id2"), features = features, logfc.threshold=0.1)
          cc
        }
}, so_pln_only)
saveRDS(comparisons_with_results, "rds/individual_clusters_comparisons_and_results_V3.rds")
saveRDS(clusters, "rds/individual_clusters_nk_bcell_clusters_V3.rds")

######################
# Generate the comparisons list object for comparing the cluster groups
clusters = unique(so[["combined_clust_grp"]])
bcell_clusters = c("All B Cells combined")
nk_clusters = c("All NK Cells combined")
clusters = list()
clusters[["all B cell clusters"]] = data.frame(cgroup = "all B cell clusters", cluster = bcell_clusters)
clusters[["NK, NK/ILC"]] = data.frame(cgroup = "NK, NK/ILC", cluster = nk_clusters)
clusters = do.call(rbind, clusters)
clusters$cidx =seq_along(rownames(clusters))
rownames(clusters) = seq_along(rownames(clusters))
clusters

goi = tibble::as_tibble(read.table("HPAP_CITEseq_gene_list_V3.csv", sep = ',', header=TRUE, row.names=NULL, fill=NA))
goi$comparison = gsub("AAb", "AAb+", goi$comparison)
exp_goi=expand_comparisons_base_short(goi, clusters)

genelist = exp_goi %>%
  mutate(feature = ifelse(modality == "ADT", adt_number, gene)) %>%
  mutate(modality = ifelse(modality == "ADT", "adt", modality))

genelist[1:20,]
genelist$cgroup
comparisons = comparisons_from_genelist(genelist)
comparisons

Seurat::Idents(so) = "Disease_Status"
so_pln_only = subset(so, Tissue %in% c("pLN-H", "pLN-T"))

comparisons_with_results = lapply(comparisons, function(cc, so) {
        print(cc)
        print(sprintf("%s : %s, %s - %s", attr(cc, "modality"), attr(cc, "id1"), attr(cc, "id2"), attr(cc, "cluster")))
        DefaultAssay(so) = attr(cc, "modality")
        print(attr(cc,"cluster"))
        so_sub = subset(so, combined_clust_grp == attr(cc, "cluster"))
        print(so_sub)
        if (length(which(Idents(so_sub) == attr(cc, "id1"))) == 0 | length(which(Idents(so_sub) == attr(cc, "id2"))) == 0) {
          attr(cc, "results") = NULL
        } else {
          features = as.character(cc)[which(cc %in% rownames(so_sub))]
          print(so_sub)
          attr(cc, "results") = Seurat::FindMarkers(so_sub, ident.1 = attr(cc, "id1"), ident.2 = attr(cc, "id2"), features = features, logfc.threshold=0.1)
          cc
        }
}, so_pln_only)
methods(Idents)
Idents(comparisons_with_results[[1]])
saveRDS(comparisons_with_results, "rds/cluster_groups_comparisons_and_results_V3.rds")
saveRDS(clusters, "rds/cluster_groups_bcell_clusters_V3.rds")
comparisons_with_results

############################################

sub
names(comparisons)
DefaultAssay(so_pln_only) = "RNA"
meh = subset(so_pln_only, manualAnnot == "NK") 
meh
FindMarkers(meh, ident.1 = "T1D", ident.2 = "ND", features = as.character(comparisons[["RNA__ND_T1D__13"]])[which(comparisons[["RNA__ND_T1D__13"]] %in% rownames(meh))], logfc.threshold=0.05)
which(!comparisons[["RNA__ND_T1D__13"]] %in% rownames(meh)) 
unique(Idents(meh))

allres = tibble()
for (comp in unique(goi$comparison)) {
  sub = goi %>% filter(comparison == comp)
  if (comp == "ND_AAb_T1D") {
    combinations =combn(c("ND", "T1D", "AAb+"), 2)
  } else {
    combinations = combn(c("ND", "T1D"),2)
  }
  for (m in c("RNA", "ADT")) {
    ss = sub %>% filter(modality == m)
    if (dim(ss)[1] == 0) {next}
    for (c in c("all B cell clusters", "NK, NK/ILC")) {
      sss = ss %>% filter(cluster == c)
      if (dim(sss)[1] == 0) {next}
      so_sub = subset(so, manualAnnot %in% clusters[[c]])
      if (m == "RNA") {
        DefaultAssay(so_sub) = "RNA"
        features = sss %>% filter(gene %in% rownames(so_sub)) %>% pull(gene)
      } else {
        DefaultAssay(so_sub) = "adt"
        features = sss %>% filter(adt_number %in% rownames(so_sub)) %>% pull(adt_number)
      }
      for (i in seq_along(combinations[1,])) {
        print(sprintf("%s -- %s", m, c))
        res = FindMarkers(so_sub, ident.1 = combinations[1,i], ident.2 = combinations[2,i], features = features, logfc.threshold=0.05)
        if (dim(res)[1] == 0) {next}
        res=as_tibble(res, rownames="feature")
        res["comparison"] = sprintf("%s_%s", combinations[1,i], combinations[2,i])
        res["modality"]=m
        res["cluster"]=c
        allres=bind_rows(allres, res)
        print(res)
      }
    }
  }
}
write.table(allres, file = "B_NK_CITEseq_gene_list_DEG.tsv", quote=FALSE, row.names=FALSE, sep='\t')


