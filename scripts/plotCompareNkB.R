library(magrittr)
library(ggplot2)
library(languageserver)
library(dplyr)
library(tidyr)
source("scripts/dimPlots.R")
source("scripts/parse_comparisons.R")
source("figures/genericFigureSettings.R")
library(Seurat)
library(grid)

#### Plotting some umaps, etc
deg = read.table("B_NK_CITEseq_gene_list_DEG.tsv", sep = '\t', header=TRUE)

bcell_rna = deg %>%
  filter(modality == "RNA") %>%
  filter(cluster == "all B cell clusters")
bcell_rna

colnames(so[[]])
pdf("figures/manualAnnot.pdf", width=11, height=11)
DimPlotCustom(so, groupBy = "manualAnnot")
dev.off()

pdf("figures/disease_status.pdf", width=11, height=11)
DimPlotCustom(so, groupBy = "Disease_Status", cols = COLORS[["disease"]])
dev.off()

pdf("figures/disease_status_hl.pdf", width=11, height=11, raster=FALSE)
patchwork::wrap_plots(smallMultipleUmaps(so, "Disease_Status", return_figures = TRUE), ncol=8)
dev.off()

pdf("figures/bcell_rna_deg.pdf", width=11, height=11)
FeaturePlotCustomPaged(so, markers = bcell_rna$feature, modality = "RNA", graphsPerRow = 4, rowsPerPage = 4, tsa_catalog=tsa_catalog)
dev.off()

sig = deg %>%
  filter(modality == "RNA") %>%
  filter(cluster == "all B cell clusters") %>%
  pivot_wider(id_cols = c("feature"), names_from="comparison", values_from="avg_log2FC") 
for (r in seq_along(sig[,1])) {
  m = matrix(nrow=3, ncol=3)
  rownames(m) = c("ND", "AAb+", "T1D")
  colnames(m) = c("ND", "AAb+", "T1D")
  for (c in seq_along(sig[1,])[-1]) {
   id1 = str_split(colnames(sig)[i],'_')[[1]][1] 
   id2 = str_split(colnames(sig)[i],'_')[[1]][2] 
   if (!is.na(sig[r,c])) {
     m[id1,id2] = 1
   } else {
     m[id1,id2] = 0
   }
  }
  print(m)
}

ht_opt$TITLE_PADDING = unit(c(8.5, 8.5), "points")


##### B cells #####
goi = tibble::as_tibble(read.table("HPAP_CITEseq_gene_list_V3.csv", sep = ',', header=TRUE, row.names=NULL, fill = NA))

so = readRDS("../seuMergedPostHSP_forFigures_2023-09-17_09-03-10.rds")

so_pln_only = readRDS("rds/so_pln_only.rds")
so_pln_only = soAddGroupedAnnotVar(so_pln_only, "manualAnnot", "groupedAnnot", convert)

clusters = unique(so[["manualAnnot"]])
bcell_clusters = clusters[which(startsWith(clusters[,1], "B")),1]
nk_clusters = c("NK", "NK/ILC")


### All b cell as a single cluster
#compres = readRDS("rds/cluster_groups_comparisons_and_results_V3.rds")
#compres
#clusters = readRDS("rds/cluster_groups_nk_bcell_clusters.rds")

compareVar = "Disease_Status"
annotVar = "manualAnnot"

### All b cell clusters handled separately 
wuv_compres_rna = readRDS("rds/wuv_compres_rna_bcellsep_genelist_V3.rds")
clusters = readRDS("rds/individual_clusters_nk_bcell_clusters.rds")

bcell_clust = clusters %>% filter(cgroup == "all B cell clusters")
nk_clust = clusters %>% filter(cgroup == "NK, NK/ILC")

Seurat::DefaultAssay(so_pln_only) = "RNA"
wuv_compres_rna = readRDS("rds/wuv_compres_rna_bcellsep_genelist_V3.rds")
rna_genes = goi[goi$Heatmap == 2,]$gene
#rna_genes = get_genes(compres, so, "RNA", bcell_clust$cluster)
rna_genes

### For Seurat object...
# colnames are cells
# rownames are features
###

# These should both do the same thing, but just keep as a Seurat oblect for now (below)... if it works
hmplt = Seurat::GetAssayData(so_pln_only, slot="data")[which(rownames(so_pln_only) %in% rna_genes),]
hmplt = Seurat::GetAssayData(so_pln_only, slot="data")[rna_genes,]

# Keep as Seurat object to make downstream simpler
hmplt = so_pln_only[rna_genes,]

# Are these really the same thing?? That sure would make me look dumb
hmplt = subset(hmplt, {{ annotVar }} %in% bcell_clust$cluster)
hmplt = hmplt[,which(colnames(so_pln_only) %in% rownames(so_pln_only[[]])[which(so_pln_only[[annotVar]][,1] %in% bcell_clust$cluster)])]

# These should both do the same thing
hmplt_data = as_tibble(t(as.matrix(hmplt)), rownames="cell") %>% left_join(as_tibble(so_pln_only[[c(annotVar, "Disease_Status")]], rownames="cell")) 
hmplt_data = seuratObjMetaTibble(hmplt, assay = "RNA")

# These should both do the same thing
rna_pctexp = percent_expressing(hmplt_data, compareVar, annotVar, zero = 0.0)
rna_pctexp = hmplt_data %>%
  select(-c(cell)) %>%
  group_by(Disease_Status, !!sym(annotVar)) %>%
  summarize_all(function(x) length(which(x != 0))/length(x)) %>%
  pivot_longer(cols=-c(annotVar, "Disease_Status"), names_to = "feature", values_to = "pctexp") %>%
  pivot_wider(names_from = "Disease_Status", values_from = "pctexp") %>%
  select(!!sym(annotVar), feature, ND, 'AAb+', T1D)

# These should both do the same thing
rna_meanExp = as.data.frame(mean_expression(hmplt_data))
hmplt_meanExp = hmplt_data %>%
  select(-c(cell)) %>%
  group_by(Disease_Status, !!sym(annotVar)) %>%
  summarize_all(mean) %>%
  as.data.frame
hmplt_meanExp

# These should both do the same thing (the first line compared to the rest of the lines
hmplt_scaledExp = scale_expression(sotib, compareVar = compareVar, by = annotVar)
ds = c("AAb+", "ND", "T1D")
bc = bcell_clust$cluster
bc
hmplt_scaledExp
hmplt_scaledExp = lapply(bc, function(xx) {
                           out = hmplt_meanExp[which(hmplt_meanExp$groupedAnnot == xx),]
                           rownames(out) = sprintf("(%s) %s", out$Disease_Status, out$groupedAnnot)
                           out = out[,-c(1,2)]
                           out = t(scale(as.matrix(out)))
                           attr(out, "cluster") = xx
                           out
})

hmplt_scaledExp

rna_pvalues = lapply(hmplt_scaledExp, fill_pvalues_from_wuv, modality="RNA")
rna_pvalues

#### ADT scaling
Seurat::DefaultAssay(so_pln_only) = "adt"
rownames(so_pln_only)
#adt_genes = get_genes(compres, so_pln_only, "adt", bcell_clust$cluster)
adt_genes = goi[goi$Heatmap == 1,]$adt_number
which(adt_genes %in% rownames(so_pln_only))

hmplt = Seurat::GetAssayData(so_pln_only, slot = "data")[which(rownames(so_pln_only) %in% adt_genes),]
hmplt = hmplt[,which(colnames(so_pln_only) %in% rownames(so_pln_only[[]])[which(so_pln_only[["groupedAnnot"]][,1] %in% bcell_clust$cluster)])]
hmplt_data = as_tibble(t(as.matrix(hmplt)), rownames="cell") %>% left_join(as_tibble(so_pln_only[[c("groupedAnnot", "Disease_Status")]], rownames="cell")) 
hmplt_data
hmplt_meanExp = hmplt_data %>%
  select(-c(cell)) %>%
  group_by(Disease_Status, groupedAnnot) %>%
  summarize_all(mean) %>%
  as.data.frame
adt_pctexp = hmplt_data %>%
  select(-c(cell)) %>%
  group_by(Disease_Status, groupedAnnot) %>%
  summarize_all(function(x) length(which(x != 0))/length(x)) %>%
  pivot_longer(cols=-c("groupedAnnot", "Disease_Status"), names_to = "feature", values_to = "pctexp") %>%
  pivot_wider(names_from = "Disease_Status", values_from = "pctexp") %>%
  select(groupedAnnot, feature, ND, 'AAb+', T1D) %>%
  mutate(feature = sapply(feature, convert_name, named_vector(goi$adt_number, goi$gene)))
hmplt_meanExp

ds = c("AAb+", "ND", "T1D")
bc = bcell_clust$cluster
bc
adt_scaledExp = lapply(bc, function(xx) {
                           out = hmplt_meanExp[which(hmplt_meanExp$groupedAnnot == xx),]
                           #rnames = hmplt_meanExp[keep,2]
                           #print(rnames)
                           rownames(out) = sprintf("(%s) %s", out$Disease_Status, out$groupedAnnot)
                           out = out[,-c(1,2)]
                           #out = out[,keep]
                           #colnames(out) = rnames
                           out = t(scale(as.matrix(out)))
                           attr(out, "cluster") = xx
                           out
})
adt_scaledExp = lapply(adt_scaledExp, function(x) {
         rownames(x) = sapply(rownames(x), function(xx) goi[goi$adt_number == xx,"gene"])
         x
})
adt_scaledExp
colmax = max(do.call(rbind, c(hmplt_scaledExp, adt_scaledExp)))
colmin = min(do.call(rbind, c(hmplt_scaledExp, adt_scaledExp)))
adt_pvalues = lapply(adt_scaledExp, fill_pvalues, modality="adt", goi=goi[goi$adt_number != "",])
adt_pvalues


####### Plot big ol patchworked heatmaps
pctexp = list()
colscale = list()
for (i in 1:length(hmplt_scaledExp)) {
  pctexp[[i]] = rna_pctexp
  colscale[[i]] = circlize::colorRamp2(breaks = c(colmin,0,colmax), rev(hcl.colors(n=24, palette = "RdBu"))[c(1,12,24)])
}
colscale
pctexp
log(3**2)
lapply(rna_pvalues, function(x) -log10(x+1e-323))
rna_plts = mapply(make_hm, hmplt_scaledExp, lapply(rna_pvalues, function(x) apply(x, 2, function(xx) -log10(xx+1e-300))), colmin= colmin, colmax=colmax, rect_label="RNA", colscale=colscale, pctexp = pctexp, rannot_labels = FALSE)  

for (i in 1:length(adt_scaledExp)) {
  pctexp[[i]] = adt_pctexp
  colscale[[i]] = circlize::colorRamp2(breaks = c(colmin,0,colmax), rev(hcl.colors(n=24, palette = "YlOrBr"))[c(1,12,24)])
}
adt_plts = mapply(make_hm, adt_scaledExp, lapply(adt_pvalues, function(x) apply(x, 2, function(xx) -log10(xx+1e-300))), colmin, colmax, "ADT", colscale, pctexp = pctexp, rannot_labels = FALSE)  
#grid.grabExpr(draw(adt_plts[[1]]))
pdf("pdf_tmp")
hm_list = lapply(seq_along(hmplt_scaledExp), function(x) {
                   ComplexHeatmap::draw(rna_plts[[x]] %v% adt_plts[[x]], padding=unit(c(1,0,1,0), "in"))
                   ComplexHeatmap::decorate_annotation("pvalues", {grid::grid.xaxis(at=1:3, label=colnames(rna_pvalues[[1]]), main=FALSE, edits = gEdit(gPath="labels", rot=90, y=unit(1.01,"npc"), just = "left"))})
                   ComplexHeatmap::decorate_annotation("percent_expressing", {grid::grid.xaxis(at=1:3, label=colnames(pctexp[[1]])[3:5], main=FALSE, edits = gEdit(gPath="labels", rot=90, y=unit(1.01,"npc"), just = "left"))})
                   print("meh")
                   grid.grab()
  })
dev.off()
#cairo_pdf("figures/bcell_diseaseStatus_hm_cluster_groups.pdf", width = 5*length(hmplt_scaledExp), height = 15)
cairo_pdf("figures/bcell_diseaseStatus_hm_individual_clusters.pdf", width = 6*length(hmplt_scaledExp), height = 25)
patchwork::wrap_plots(hm_list, nrow=1)
dev.off()

##### Second edition of heatmaps (smaller gene set, split up into multiple heatmaps
goiV2 = tibble::as_tibble(read.table("HPAP_CITEseq_gene_list_V2.csv", sep = ',', header=TRUE, row.names=NULL, fill=NA))

rna_genes = goiV2 %>% filter(modality == "RNA")
adt_genes = goiV2 %>% filter(modality == "ADT")

hma_colscale = list("percent_expressing" = circlize::colorRamp2(breaks=c(0,1), rev(hcl.colors(n=24, palette = "Purples"))[c(1,24)]),
                    "pvalues" = circlize::colorRamp2(breaks=c(0,10,50,100,323), hcl.colors(n=24, palette = "Earth")[c(6,12,16,20,24)]))
lgds = list("percent_expressing" = ComplexHeatmap::Legend(col_fun = hma_colscale[["percent_expressing"]], title = "% expr"),
            "pvalues" = ComplexHeatmap::Legend(col_fun = hma_colscale[["pvalues"]], title = "-log10(p)", at = c(0,10,50,100,323), break_dist=c(10,5,2,3.23)))
heatmaps =sort(as.numeric(unique(rna_genes$Heatmap)))
for (i in seq_along(heatmaps)) {
  toplt = hmplt_scaledExp[[1]][-which(is.na(hmplt_scaledExp[[1]][,1])),]
  toplt = toplt[which(rownames(toplt) %in% (rna_genes %>% filter(Heatmap == heatmaps[i]) %>% pull(gene))),]
  toplt_p = rna_pvalues[[1]][which(rownames(rna_pvalues[[1]]) %in% rownames(toplt)),]
  toplt_pe = rna_pctexp[which(rna_pctexp$feature %in% rownames(toplt)),]
  colmax = max(toplt)
  colmin = min(toplt)
  colscale=circlize::colorRamp2(breaks = c(colmin,0,colmax), rev(hcl.colors(n=24, palette = "RdBu"))[c(1,12,24)])
  cairo_pdf(sprintf("figures/bcells_singleClust_diseaseStatus_RNA_%s.pdf", i), width = 8.5, height = 8.5*(length(seq_along(rownames(toplt)))/30))
  draw(make_hm(toplt, apply(toplt_p, 2, function(xx) -log10(xx+1e-300)), colmin= colmin, colmax=colmax, rect_label="RNA", colscale=colscale, pctexp = toplt_pe, rannot_labels = FALSE, showlegends=TRUE))
  draw(lgds[["percent_expressing"]], x = unit(0.94, "npc"), y = unit(0.8, "npc"))
  draw(lgds[["pvalues"]], x = unit(0.94, "npc"), y = unit(0.47, "npc"))
  dev.off()
}

heatmaps =sort(as.numeric(unique(adt_genes$Heatmap)))
for (i in seq_along(heatmaps)) {
  torm = which(is.na(adt_scaledExp[[1]][,1]))
  length(torm)
  if (length(torm)>0) {
    toplt = adt_scaledExp[[1]][-torm,]
  } else {
    toplt = adt_scaledExp[[1]]
  }
  toplt = toplt[which(rownames(toplt) %in% (adt_genes %>% filter(Heatmap == heatmaps[i]) %>% pull(gene))),]
  toplt_p = adt_pvalues[[1]][which(rownames(adt_pvalues[[1]]) %in% rownames(toplt)),]
  toplt_pe = adt_pctexp[which(adt_pctexp$feature %in% rownames(toplt)),]
  colmax = max(toplt)
  colmin = min(toplt)
  colscale=circlize::colorRamp2(breaks = c(colmin,0,colmax), rev(hcl.colors(n=24, palette = "YlOrBr"))[c(1,12,24)])
  cairo_pdf(sprintf("figures/bcells_singleClust_diseaseStatus_ADT_%s.pdf", i), width = 8.5, height = 8.5*(length(seq_along(rownames(toplt)))/20))
  draw(make_hm(toplt, apply(toplt_p, 2, function(xx) -log10(xx+1e-300)), colmin= colmin, colmax=colmax, rect_label="RNA", colscale=colscale, pctexp = toplt_pe, rannot_labels = TRUE, showlegends=TRUE), padding=unit(c(10,10,10,10), "in"))
  draw(lgds[["percent_expressing"]], x = unit(0.94, "npc"), y = unit(0.20, "npc"))
  draw(lgds[["pvalues"]], x = unit(0.94, "npc"), y = unit(0.47, "npc"))
  dev.off()
}

##### Third edition of heatmaps (smaller gene set, split up into multiple heatmaps
goiV3 = tibble::as_tibble(read.table("HPAP_CITEseq_gene_list_V3.csv", sep = ',', header=TRUE, row.names=NULL, fill=NA))
unique(goiV3$Heatmap)
rna_genes = goiV3 %>% filter(modality == "RNA")
rna_genes
adt_genes = goiV3 %>% filter(modality == "ADT")

hma_colscale = list("percent_expressing" = circlize::colorRamp2(breaks=c(0,1), rev(hcl.colors(n=24, palette = "Purples"))[c(1,24)]),
                    "pvalues" = circlize::colorRamp2(breaks=c(0,10,50,100,323), hcl.colors(n=24, palette = "Earth")[c(6,12,16,20,24)]))
lgds = list("percent_expressing" = ComplexHeatmap::Legend(col_fun = hma_colscale[["percent_expressing"]], title = "% expr"),
            "pvalues" = ComplexHeatmap::Legend(col_fun = hma_colscale[["pvalues"]], title = "-log10(p)", at = c(0,10,50,100,323), break_dist=c(10,5,2,3.23)))
heatmaps =sort(as.numeric(unique(rna_genes$Heatmap)))
heatmaps
hmplt_scaledExp
for (i in seq_along(heatmaps)) {
  torm = which(is.na(hmplt_scaledExp[[1]][,1]))
  length(torm)
  if (length(torm)>0) {
    toplt = hmplt_scaledExp[[1]][-torm,]
  } else {
    toplt = hmplt_scaledExp[[1]]
  }
  toplt
  toplt = toplt[which(rownames(toplt) %in% (rna_genes %>% filter(Heatmap == heatmaps[i]) %>% pull(gene))),]
  toplt
  toplt_p = rna_pvalues[[1]][which(rownames(rna_pvalues[[1]]) %in% rownames(toplt)),]
  toplt_pe = rna_pctexp[which(rna_pctexp$feature %in% rownames(toplt)),]
  colmax = max(toplt)
  colmin = min(toplt)
  colscale=circlize::colorRamp2(breaks = c(colmin,0,colmax), rev(hcl.colors(n=24, palette = "RdBu"))[c(1,12,24)])
  cairo_pdf(sprintf("/srv/http/betts/hpap/figures/bcells_singleClust_diseaseStatus_RNA_V3_%s.pdf", i), width = 8.5, height = 8.5*(length(seq_along(rownames(toplt)))/10))
  draw(make_hm(toplt, apply(toplt_p, 2, function(xx) -log10(xx+1e-300)), colmin= colmin, colmax=colmax, annotVar = "groupedAnnot", rect_label="RNA", colscale=colscale, pctexp = toplt_pe, rannot_labels = FALSE, showlegends=TRUE), padding=unit(c(1,0.1,1,0.1), "in"))
  ComplexHeatmap::decorate_annotation("pvalues", {grid::grid.xaxis(at=1:3, label=colnames(toplt_p), main=FALSE, edits = gEdit(gPath="labels", rot=90, y=unit(1.01,"npc"), just = "left"))})
  ComplexHeatmap::decorate_annotation("percent_expressing", {grid::grid.xaxis(at=1:3, label=colnames(toplt_pe)[3:5], main=FALSE, edits = gEdit(gPath="labels", rot=90, y=unit(1.01,"npc"), just = "left"))})
  draw(lgds[["percent_expressing"]], x = unit(0.94, "npc"), y = unit(0.2, "npc"))
  draw(lgds[["pvalues"]], x = unit(0.94, "npc"), y = unit(0.47, "npc"))
  dev.off()
}

heatmaps =sort(as.numeric(unique(adt_genes$Heatmap)))
for (i in seq_along(heatmaps)) {
  torm = which(is.na(adt_scaledExp[[1]][,1]))
  length(torm)
  if (length(torm)>0) {
    toplt = adt_scaledExp[[1]][-torm,]
  } else {
    toplt = adt_scaledExp[[1]]
  }
  toplt = toplt[which(rownames(toplt) %in% (adt_genes %>% filter(Heatmap == heatmaps[i]) %>% pull(gene))),]
  toplt_p = adt_pvalues[[1]][which(rownames(adt_pvalues[[1]]) %in% rownames(toplt)),]
  toplt_pe = adt_pctexp[which(adt_pctexp$feature %in% rownames(toplt)),]
  colmax = max(toplt)
  colmin = min(toplt)
  colscale=circlize::colorRamp2(breaks = c(colmin,0,colmax), rev(hcl.colors(n=24, palette = "YlOrBr"))[c(1,12,24)])
  cairo_pdf(sprintf("/srv/http/betts/hpap/figures/bcells_singleClust_diseaseStatus_ADT_V3_%s.pdf", i), width = 8.5, height = 8.5*(length(seq_along(rownames(toplt)))/10))
  draw(make_hm(toplt, apply(toplt_p, 2, function(xx) -log10(xx+1e-300)), colmin= colmin, colmax=colmax, rect_label="ADT", colscale=colscale, annotVar = "groupedAnnot", pctexp = toplt_pe, rannot_labels = FALSE, showlegends=TRUE), padding = unit(c(1, 0.1, 1, 0.1), "in"))
  ComplexHeatmap::decorate_annotation("pvalues", {grid::grid.xaxis(at=1:3, label=colnames(toplt_p), main=FALSE, edits = gEdit(gPath="labels", rot=90, y=unit(1.01,"npc"), just = "left"))})
  ComplexHeatmap::decorate_annotation("percent_expressing", {grid::grid.xaxis(at=1:3, label=colnames(toplt_pe)[3:5], main=FALSE, edits = gEdit(gPath="labels", rot=90, y=unit(1.01,"npc"), just = "left"))})
  draw(lgds[["percent_expressing"]], x = unit(0.94, "npc"), y = unit(0.20, "npc"))
  draw(lgds[["pvalues"]], x = unit(0.94, "npc"), y = unit(0.47, "npc"))
  dev.off()
}

###### NK Cells #####
compres = readRDS("rds/cluster_groups_comparisons_and_results.rds")
Seurat::DefaultAssay(so) = "RNA"
rna_genes = get_genes(compres, so, "RNA", nk_clust$cluster)
rna_genes

compres
hmplt = Seurat::GetAssayData(so, slot="data")[which(rownames(so) %in% rna_genes),]
hmplt = hmplt[,which(colnames(so) %in% rownames(so[[]])[which(so[["manualAnnot"]][,1] %in% nk_clust$cluster)])]
hmplt_data = as_tibble(t(as.matrix(hmplt)), rownames="cell") %>% left_join(as_tibble(so[[c("manualAnnot", "Disease_Status")]], rownames="cell")) 
hmplt_data
rna_pctexp = hmplt_data %>%
  select(-c(cell)) %>%
  group_by(Disease_Status, manualAnnot) %>%
  summarize_all(function(x) length(which(x != 0))/length(x)) %>%
  pivot_longer(cols=-c("manualAnnot", "Disease_Status"), names_to = "feature", values_to = "pctexp") %>%
  pivot_wider(names_from = "Disease_Status", values_from = "pctexp") %>%
  select(manualAnnot, feature, ND, 'AAb+', T1D)
hmplt_meanExp = hmplt_data %>%
  select(-c(cell)) %>%
  group_by(Disease_Status, manualAnnot) %>%
  summarize_all(mean) %>%
  as.data.frame
hmplt_meanExp

ds = c("AAb+", "ND", "T1D")
bc = nk_clust$cluster
bc
hmplt_scaledExp = lapply(bc, function(xx) {
                           out = hmplt_meanExp[which(hmplt_meanExp$manualAnnot == xx),]
                           #rnames = hmplt_meanExp[keep,2]
                           #print(rnames)
                           rownames(out) = sprintf("(%s) %s", out$Disease_Status, out$manualAnnot)
                           out = out[,-c(1,2)]
                           #out = out[,keep]
                           #colnames(out) = rnames
                           out = t(scale(as.matrix(out)))
                           attr(out, "cluster") = xx
                           out
})
rna_pvalues = lapply(hmplt_scaledExp, fill_pvalues, modality="RNA")
rna_pvalues

#### ADT scaling
Seurat::DefaultAssay(so) = "adt"
adt_genes = get_genes(compres, so, "adt", nk_clust$cluster)
adt_genes

adt_genes
hmplt = Seurat::GetAssayData(so, slot="data")[which(rownames(so) %in% adt_genes),]
hmplt = hmplt[,which(colnames(so) %in% rownames(so[[]])[which(so[["manualAnnot"]][,1] %in% nk_clust$cluster)])]
hmplt_data = as_tibble(t(as.matrix(hmplt)), rownames="cell") %>% left_join(as_tibble(so[[c("manualAnnot", "Disease_Status")]], rownames="cell")) 
hmplt_data
hmplt_meanExp = hmplt_data %>%
  select(-c(cell)) %>%
  group_by(Disease_Status, manualAnnot) %>%
  summarize_all(mean) %>%
  as.data.frame
adt_pctexp = hmplt_data %>%
  select(-c(cell)) %>%
  group_by(Disease_Status, manualAnnot) %>%
  summarize_all(function(x) length(which(x != 0))/length(x)) %>%
  pivot_longer(cols=-c("manualAnnot", "Disease_Status"), names_to = "feature", values_to = "pctexp") %>%
  pivot_wider(names_from = "Disease_Status", values_from = "pctexp") %>%
  select(manualAnnot, feature, ND, 'AAb+', T1D)
hmplt_meanExp

ds = c("AAb+", "ND", "T1D")
bc = nk_clust$cluster
bc
adt_scaledExp = lapply(bc, function(xx) {
                           out = hmplt_meanExp[which(hmplt_meanExp$manualAnnot == xx),]
                           #rnames = hmplt_meanExp[keep,2]
                           #print(rnames)
                           rownames(out) = sprintf("(%s) %s", out$Disease_Status, out$manualAnnot)
                           out = out[,-c(1,2)]
                           #out = out[,keep]
                           #colnames(out) = rnames
                           out = t(scale(as.matrix(out)))
                           attr(out, "cluster") = xx
                           out
})
adt_scaledExp = lapply(adt_scaledExp, function(x) {
         rownames(x) = sapply(rownames(x), function(xx) goi[goi$adt_number == xx,"gene"])
         x
})
adt_pvalues = lapply(adt_scaledExp, fill_pvalues, modality="adt", goi=goi[goi$adt_number != "",])
adt_pvalues
#######
todrop = lapply(hmplt_scaledExp, function(x) rownames(x)[which(is.na(x[,1]))])
todrop
length(hmplt_scaledExp)
if(length(hmplt_scaledExp) > 1) {
  hmplt_scaledExp = mapply(function(x, y) x[-which(rownames(x) %in% y),], hmplt_scaledExp, todrop)
  rna_pvalues=mapply(function(x, y) x[-which(rownames(x) %in% y),], rna_pvalues, todrop)
} else {
  hmplt_scaledExp = lapply(hmplt_scaledExp, function(x, y) x[-which(rownames(x) %in% y),], y=todrop[[1]])
  rna_pvalues = lapply(rna_pvalues, function(x, y) x[-which(rownames(x) %in% y),], y=todrop[[1]])
}
colmax = max(do.call(rbind, c(hmplt_scaledExp, adt_scaledExp)))
colmin = min(do.call(rbind, c(hmplt_scaledExp, adt_scaledExp)))
pctexp = list()
colscale = list()
rna_pvalues
for (i in 1:length(hmplt_scaledExp)) {
  pctexp[[i]] = rna_pctexp %>% filter(!feature %in% todrop[[i]])
  colscale[[i]] = circlize::colorRamp2(breaks = c(colmin,0,colmax), rev(hcl.colors(n=24, palette = "RdBu"))[c(1,12,24)])
}
rna_pvalues
pctexp
rna_plts = mapply(make_hm, hmplt_scaledExp, lapply(rna_pvalues, function(x) apply(x, 2, function(xx) -log10(xx+1e-300))), colmin= colmin, colmax=colmax, rect_label="RNA", colscale=colscale, pctexp = pctexp, rannot_labels = FALSE)  

for (i in 1:length(hmplt_scaledExp)) {
  pctexp[[i]] = adt_pctexp
  colscale[[i]] = circlize::colorRamp2(breaks = c(colmin,0,colmax), rev(hcl.colors(n=24, palette = "YlOrBr"))[c(1,12,24)])
}
adt_plts = mapply(make_hm, adt_scaledExp, lapply(adt_pvalues, function(x) apply(x, 2, function(xx) -log10(xx+1e-300))), colmin, colmax, "ADT", colscale, pctexp = pctexp, rannot_labels = FALSE)  
#grid.grabExpr(draw(adt_plts[[1]]))
rna_plts
hm_list = lapply(seq_along(hmplt_scaledExp), function(x) {
                   draw(rna_plts[[x]] %v% adt_plts[[x]], padding=unit(c(1,0,1,0), "in"))
                   ComplexHeatmap::decorate_annotation("pvalues", {grid::grid.xaxis(at=1:3, label=colnames(rna_pvalues[[1]]), main=FALSE, edits = gEdit(gPath="labels", rot=90, y=unit(1.01,"npc"), just = "left"))})
                   ComplexHeatmap::decorate_annotation("percent_expressing", {grid::grid.xaxis(at=1:3, label=colnames(pctexp[[1]])[3:5], main=FALSE, edits = gEdit(gPath="labels", rot=90, y=unit(1.01,"npc"), just = "left"))})
                   print("meh")
                   grid.grab()
  })
hm_list
#cairo_pdf("figures/nk_diseaseStatus_hm_cluster_groups.pdf", width = 5*length(hmplt_scaledExp), height = 30)
cairo_pdf("figures/nk_diseaseStatus_hm_individual_clusters.pdf", width = 5*length(hmplt_scaledExp), height = 30)
patchwork::wrap_plots(hm_list, nrow=1)
dev.off()

#### rank plots
nd_t1d_se = compres$RNA__ND_T1D__2
g1 = c("GRB2", "CD8A", "LAT", "NKG7", "CD7", "MICA", "HCST", "CD247", "ULBP3", "SH3BP2", "HLA-A", "ITGB2", "SRGN")
g2 = c()
attr(nd_t1d_se, "results")[which(rownames(attr(nd_t1d_se, "results")) %in% ),]














########################################
########################################
Seurat::DefaultAssay(so) = "RNA"
rna_genes = get_genes(compres, so, "RNA", nk_clust$cluster)
rna_genes

hmplt = Seurat::GetAssayData(so, slot="data")[which(rownames(so) %in% rna_genes),]
hmplt = hmplt[,which(colnames(so) %in% rownames(so[[]])[which(so[["manualAnnot"]][,1] %in% nk_clust$cluster)])]
hmplt_data = as_tibble(t(as.matrix(hmplt)), rownames="cell") %>% left_join(as_tibble(so[[c("manualAnnot", "Disease_Status")]], rownames="cell")) 
hmplt_data
hmplt_meanExp = hmplt_data %>%
  select(-c(cell)) %>%
  group_by(Disease_Status, manualAnnot) %>%
  summarize_all(mean) %>%
  as.data.frame
hmplt_meanExp

ds = c("AAb+", "ND", "T1D")
bc = nk_clust$cluster
bc

hmplt_meanExp

hmplt_scaledExp = lapply(bc, function(xx) {
                           out = hmplt_meanExp[which(hmplt_meanExp$manualAnnot == xx),]
                           #rnames = hmplt_meanExp[keep,2]
                           #print(rnames)
                           rownames(out) = sprintf("(%s) %s", out$Disease_Status, out$manualAnnot)
                           out = out[,-c(1,2)]
                           #out = out[,keep]
                           #colnames(out) = rnames
                           out = t(scale(as.matrix(out)))
                           out
})

# Get rid of NAs from counts that are all 0
hmplt_scaledExp = lapply(hmplt_scaledExp, function(x) {
                           x[which(is.na(x[,1])),] = 0.0
                           x
})
colmax = max(do.call(rbind, c(hmplt_scaledExp, adt_scaledExp)))
colmin = min(do.call(rbind, c(hmplt_scaledExp, adt_scaledExp)))
#### ADT scaling
Seurat::DefaultAssay(so) = "adt"
adt_genes = get_genes(compres, so, "adt", nk_clust$cluster)
adt_genes

adt_genes
hmplt = Seurat::GetAssayData(so, slot="data")[which(rownames(so) %in% adt_genes),]
hmplt = hmplt[,which(colnames(so) %in% rownames(so[[]])[which(so[["manualAnnot"]][,1] %in% nk_clust$cluster)])]
hmplt_data = as_tibble(t(as.matrix(hmplt)), rownames="cell") %>% left_join(as_tibble(so[[c("manualAnnot", "Disease_Status")]], rownames="cell")) 
hmplt_data
hmplt_meanExp = hmplt_data %>%
  select(-c(cell)) %>%
  group_by(Disease_Status, manualAnnot) %>%
  summarize_all(mean) %>%
  as.data.frame
hmplt_meanExp

ds = c("AAb+", "ND", "T1D")
bc = nk_clust$cluster
bc
adt_scaledExp = lapply(bc, function(xx) {
                           out = hmplt_meanExp[which(hmplt_meanExp$manualAnnot == xx),]
                           #rnames = hmplt_meanExp[keep,2]
                           #print(rnames)
                           rownames(out) = sprintf("(%s) %s", out$Disease_Status, out$manualAnnot)
                           out = out[,-c(1,2)]
                           #out = out[,keep]
                           #colnames(out) = rnames
                           out = t(scale(as.matrix(out)))
                           out
})
adt_scaledExp 
goi
adt_scaledExp = lapply(adt_scaledExp, function(x) {
         rownames(x) = sapply(rownames(x), function(xx) goi[goi$adt_number == xx,"gene"])
         x
})
#######

#max_fold_change = lapply(hmplt_meanExp, function(xx) as.matrix(apply(xx, 1, function(x) log2(max(x)/min(x)))))
make_rna_hm = function(scaled_data, bottom_scaled_data, print_annot_labels, colmin, colmax) {
  #print(dim(as.data.frame(bottom_scaled_data)))
  #print(as.data.frame(t(bottom_scaled_data)))
  #print(dim(as.data.frame(t(bottom_scaled_data))))
  #print(ncol(as.matrix(scaled_data)))
  column_order = unname(sapply(colnames(scaled_data), function(cn) which(grepl(strsplit(cn, " ")[[1]][1], c("ND", "AAb+", "T1D")))))
  colscale = circlize::colorRamp2(c(colmin,0,colmax), rev(hcl.colors(n=24, palette = "RdBu"))[c(1,12,24)])
  colscale2 = circlize::colorRamp2(c(colmin,0,colmax), rev(hcl.colors(n=24, palette = "YlOrBr")[c(1,12,24)]))
  bottom_annot_col = lapply(rownames(bottom_scaled_data), function(x) colscale2)
  names(bottom_annot_col) = as.character(rownames(bottom_scaled_data))
  #print(bottom_annot_col)

  ComplexHeatmap::Heatmap(scaled_data,
          na_col = "grey90", 
          cluster_rows=T, 
          cluster_columns=F, 
          rect_gp = gpar(col = "white", lwd = 2),
          col = colscale,
          row_title = "RNA", 
          row_title_gp = gpar(fill = "black", col = "white", border = "black"),
          #cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
          #  grid.text(round(scaled_data, 3)[i, j], x, y)
          #u},
          show_heatmap_legend = FALSE,
          border = FALSE,"
          column_order = column_order,
          bottom_annotation = ComplexHeatmap::columnAnnotation(df=as.data.frame(t(bottom_scaled_data)), 
                                                               col = bottom_annot_col,
                                                               show_annotation_name=print_annot_labels,
                                                               show_legend=FALSE,
                                                               simple_anno_size=unit(.2,"in"),
                                                               gp= gpar(col = "white", lwd = 2),
                                                              border = FALSE)
)
}
plts = mapply(make_rna_hm, hmplt_scaledExp, adt_scaledExp, c(rep(FALSE, length(hmplt_scaledExp)-1), TRUE), colmin, colmax)
length(plts)
nk_hm_list = plts[[1]] + plts[[2]] 
nk_hm_list = plts[[1]] + plts[[2]] + plts[[3]] + plts[[4]] + plts[[5]] + plts[[6]] +plts[[7]] + plts[[8]] + plts[[9]] + plts[[10]] + plts[[11]] + plts[[12]] 
pdf("figures/nk_diseaseStatus_hm.pdf", width = 30, height = 20)
draw(nk_hm_list, ht_gap = unit("0.5", "cm"))
dev.off()







#max_fold_change = lapply(hmplt_meanExp, function(xx) as.matrix(apply(xx, 1, function(x) log2(max(x)/min(x)))))
make_adt_hm = function(scaled_data, max_fold_change) {
  column_order = unname(sapply(colnames(scaled_data), function(cn) which(grepl(strsplit(cn, " ")[[1]][1], c("ND", "AAb+", "T1D")))))
  Heatmap(scaled_data,
          na_col = "grey90", 
          cluster_rows=T, 
          cluster_columns=F, 
          rect_gp = gpar(col = "white", lwd = 2), 
          col = rev(hcl.colors(n=25, palette = "Tropic")),
          row_title = "ADT",
          row_title_gp = gpar(fill = "black", col = "white", border = "black"),
          cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
            grid.text(round(scaled_data, 3)[i, j], x, y)
          },
          column_order = column_order,
          #left_annotation = rowAnnotation(MaxFC=max_fold_change)
  )
}

adt_plts = lapply(hmplt_scaledExp, make_adt_hm)


#####################
bcell_rna = deg %>%
  filter(modality == "RNA") %>%
  filter(cluster == "all B cell clusters")
bcell_rna
Seurat::DefaultAssay(so) = "RNA"
hmplt = Seurat::GetAssayData(so, slot="scale.data")[which(rownames(so) %in% bcell_rna$f(comparisons[["RNA__ND_T1D__13"]]eature),]
hmplt_meanExp = as_tibble(t(as.matrix(hmplt)), rownames="cell") %>% left_join(as_tibble(so[[c("manualAnnot", "Disease_Status")]], rownames="cell")) %>%
  filter(manualAnnot %in% clusters[["all B cell clusters"]]) %>%
  select(-c(cell, manualAnnot)) %>%
  group_by(Disease_Status) %>%
  summarize_all(mean) %>%
  as.data.frame
rownames(hmplt_meanExp)=hmplt_meanExp[,1]
hmplt_meanExp=hmplt_meanExp[,-1]

hmplt_meanExp_scaled_rna = t(scale(as.matrix(hmplt_meanExp)))
max_fold_change = apply(hmplt_meanExp, 2, function(x) log2(max(x)/min(x)))
bcell_rna_hm = Heatmap(hmplt_meanExp_scaled_rna, 
        na_col = "grey90", 
        cluster_rows=T, 
        cluster_columns=F, 
        rect_gp = gpar(col = "white", lwd = 2), 
        col = rev(hcl.colors(n=25, palette = "RdBu")),
        row_title = "RNA", 
        row_title_gp = gpar(fill = "black", col = "white", border = "black"),
        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
          grid.text(round(hmplt_meanExp_scaled_rna, 3)[i, j], x, y)
        },
        column_order = c("ND", "AAb+", "T1D"),
        left_annotation = rowAnnotation(MaxFC=max_fold_change)
        
)


bcell_adt = deg %>%
  filter(modality == "ADT") %>%
  filter(cluster == "all B cell clusters")
bcell_adt
Seurat::DefaultAssay(so) = "adt"
hmplt = Seurat::GetAssayData(so, slot="data")[which(rownames(so) %in% bcell_adt$feature),]
hmplt_meanExp = as_tibble(t(as.matrix(hmplt)), rownames="cell") %>% left_join(as_tibble(so[[c("manualAnnot", "Disease_Status")]], rownames="cell")) %>%
  filter(manualAnnot %in% clusters[["all B cell clusters"]]) %>%
  select(-c(cell, manualAnnot)) %>%
  group_by(Disease_Status) %>%
  summarize_all(mean) %>%
  as.data.frame
rownames(hmplt_meanExp)=hmplt_meanExp[,1]
hmplt_meanExp=hmplt_meanExp[,-1]

hmplt_meanExp_scaled_adt = t(scale(as.matrix(hmplt_meanExp)))
sig = bcell_adt %>%
  separate(comparison, into = c("id1", "id2"), sep = "_")

max_fold_change = apply(hmplt_meanExp, 2, function(x) log2(max(x)/min(x)))
bcell_adt_hm = Heatmap(hmplt_meanExp_scaled_adt, 
        na_col = "grey90", 
        cluster_rows=T, 
        cluster_columns=F, 
        rect_gp = gpar(col = "white", lwd = 2), 
        col = hcl.colors(n=25, palette = "Tropic"),
        row_title = "ADT", 
        row_title_gp = gpar(fill = "black", col = "white", border = "black"),
        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
          grid.text(round(hmplt_meanExp_scaled_adt, 3)[i, j], x, y, just="left")
          grid.text(round(hmplt_meanExp_scaled_adt, 5)[i, j], x, y, just="right")
        },
        column_order = c("ND", "AAb+", "T1D"),
        left_annotation = rowAnnotation(MaxFC=max_fold_change)
        
)
pdf("figures/disease_state_hm.pdf", width = 8.5, height = 11)
bcell_rna_hm %v% bcell_adt_hm
dev.off()


###########################################
###########################################


nk_rna = deg %>%
  filter(modality == "RNA") %>%
  filter(cluster == "NK, NK/ILC")
nk_rna
Seurat::DefaultAssay(so) = "RNA"
hmplt = Seurat::GetAssayData(so, slot="data")[which(rownames(so) %in% nk_rna$feature),]
hmplt_meanExp = as_tibble(t(as.matrix(hmplt)), rownames="cell") %>% left_join(as_tibble(so[[c("manualAnnot", "Disease_Status")]], rownames="cell")) %>%
  filter(manualAnnot %in% clusters[["all B cell clusters"]]) %>%
  select(-c(cell, manualAnnot)) %>%
  group_by(Disease_Status) %>%
  summarize_all(mean) %>%
  as.data.frame
rownames(hmplt_meanExp)=hmplt_meanExp[,1]
hmplt_meanExp=hmplt_meanExp[,-1]

hmplt_meanExp_scaled_rna = t(scale(as.matrix(hmplt_meanExp)))
max_fold_change = apply(hmplt_meanExp, 2, function(x) log2(max(x)/min(x)))
nk_rna_hm = Heatmap(hmplt_meanExp_scaled_rna, 
        na_col = "grey90", 
        cluster_rows=T, 
        cluster_columns=F, 
        rect_gp = gpar(col = "white", lwd = 2), 
        col = hcl.colors(n=25, palette = "Purple-Yellow"),
        row_title = "RNA", 
        row_title_gp = gpar(fill = "black", col = "white", border = "black"),
        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
          grid.text(round(hmplt_meanExp_scaled_rna, 3)[i, j], x, y, gp =gpar(fontsize = 7))
        },
        column_order = c("ND", "AAb+", "T1D"),
        left_annotation = rowAnnotation(MaxFC=max_fold_change)
        
)


nk_adt = deg %>%
  filter(modality == "ADT") %>%
  filter(cluster == "NK, NK/ILC")
nk_adt
Seurat::DefaultAssay(so) = "adt"
hmplt = Seurat::GetAssayData(so, slot="data")[which(rownames(so) %in% nk_adt$feature),]
hmplt_meanExp = as_tibble(t(as.matrix(hmplt)), rownames="cell") %>% left_join(as_tibble(so[[c("manualAnnot", "Disease_Status")]], rownames="cell")) %>%
  filter(manualAnnot %in% clusters[["all B cell clusters"]]) %>%
  select(-c(cell, manualAnnot)) %>%
  group_by(Disease_Status) %>%
  summarize_all(mean) %>%
  as.data.frame
rownames(hmplt_meanExp)=hmplt_meanExp[,1]
hmplt_meanExp=hmplt_meanExp[,-1]

hmplt_meanExp_scaled_adt = t(scale(as.matrix(hmplt_meanExp)))
sig = nk_adt %>%
  separate(comparison, into = c("id1", "id2"), sep = "_")

max_fold_change = apply(hmplt_meanExp, 2, function(x) log2(max(x)/min(x)))
nk_adt_hm = Heatmap(hmplt_meanExp_scaled_adt, 
        na_col = "grey90", 
        cluster_rows=T, 
        cluster_columns=F, 
        rect_gp = gpar(col = "white", lwd = 2), 
        col = hcl.colors(n=25, palette = "Heat"),
        row_title = "ADT", 
        row_title_gp = gpar(fill = "black", col = "white", border = "black"),
        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
          grid.text(round(hmplt_meanExp_scaled_adt, 3)[i, j], x, y, gp =gpar(fontsize = 7))
        },
        column_order = c("ND", "AAb+", "T1D"),
        left_annotation = rowAnnotation(MaxFC=max_fold_change)
        
)
pdf("figures/nk_disease_state_hm.pdf", width = 8.5, height = 11)
nk_rna_hm %v% nk_adt_hm
dev.off()
