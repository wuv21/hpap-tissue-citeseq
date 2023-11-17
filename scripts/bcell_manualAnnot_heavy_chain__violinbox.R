library(Seurat)
library(dplyr)
library(magrittr)
library(tidyr)
library(ggplot2)
library(ggpubr)
source("scripts/so_helpers.R")
source("scripts/transform_expression.R")
source("scripts/generic.R")
source("figures/genericFigureSettings.R")

DISEASESTATUSCOLORS

ASSAY="RNA"
COMPAREVAR="Disease_Status"
ANNOTVAR="manualAnnot"
ANNOTVARINCLUDE = c("B class-switched memory #1", "B IgM+ memory/marginal zone like #2 ")

so_pln_only = readRDS("rds/so_pln_only.rds")
compres = readRDS("rds/wuv_compres_rna_bcellsep_genelist_V3.rds")

clusters = unique(so_pln_only[["manualAnnot"]])
bcell_clusters = clusters[which(startsWith(clusters[,1], "B")),1]
nk_clusters = c("NK", "NK/ILC")
convert = c(rep("All B cells combined", length(bcell_clusters)), rep("All NK Cells combined", length(nk_clusters)))
names(convert)=c(bcell_clusters, nk_clusters)
so_pln_only = soAddGroupedAnnotVar(so_pln_only, "manualAnnot", "groupedAnnot", convert)

#goiV3 = tibble::as_tibble(read.table("HPAP_CITEseq_gene_list_V3.csv", sep = ',', header=TRUE, row.names=NULL, fill=NA))
#goiV3$gene
#genes_of_interest = goiV3[goiV3$modality == "RNA",]$gene[-c(1, 11, 12, 13, 14)]
genes_of_interest = c("IGHG1", "IGHM", "IGHA1", "IGHA2", "MYC", "KLF2")
stopifnot(all(genes_of_interest %in% rownames(so_pln_only)))
unique(so_pln_only[[ANNOTVAR]][,1])
stopifnot(all(ANNOTVARINCLUDE %in% unique(so_pln_only[[ANNOTVAR]][,1])))

compres = lapply(compres, function(x, goi) x[x$gene %in% goi, ], genes_of_interest)
Seurat::DefaultAssay(so_pln_only) = ASSAY

pltData = so_pln_only[genes_of_interest,]

pltData = subset(pltData, !!sym(ANNOTVAR) %in% ANNOTVARINCLUDE)
dim(pltData)

pltData = seuratObjMetaTibble(pltData, assay = ASSAY)
attr(pltData, "datacol") = seq_along(genes_of_interest)+1
rna_pctexp = percent_expressing(pltData, COMPAREVAR, ANNOTVAR, zero = 0.0)
rna_meanExp = as.data.frame(mean_expression(pltData, COMPAREVAR, ANNOTVAR))
rna_stderr_bars = as.data.frame(expression_stderr_bars(pltData, COMPAREVAR, ANNOTVAR, method = mean, bar_method = sd))
rna_scaledExp = scale_expression(rna_meanExp, compareVar = COMPAREVAR, by = ANNOTVAR)
rna_scaledExp

rna_pvalues = lapply(rna_scaledExp, fillPvalsFromFindFeaturesCombinatorial, compres, modality="RNA")
names(rna_pvalues) = sapply(rna_scaledExp, function(x) attr(x, "cluster"))
rna_pvalues

vioplt = pltData %>% 
  select("cell", COMPAREVAR, ANNOTVAR, colnames(pltData)[attr(pltData, "datacol")]) %>%
  pivot_longer(cols = !c("cell", COMPAREVAR, ANNOTVAR), names_to = "Gene", values_to = "normalized_expression") %>%
  mutate(Disease_Status = factor(Disease_Status, levels = c("ND", "AAb+", "T1D")))
vioplt
vioplt[[ANNOTVAR]] = factor(vioplt[[ANNOTVAR]], levels=ANNOTVARINCLUDE) 
vioplt[[ANNOTVAR]] 

bracketx = function(offset, ds) {
  coords = c("ND" = offset+(1-0.3), "AAb+" = offset+(1.0), "T1D" = offset+(1+0.3))
  coords[ds]
}

pdf("/srv/http/betts/hpap/figures/bcell_sep_heavy_chain.violins.pdf", width = 11, height=8)
for (g in genes_of_interest) {
  pltdata =vioplt %>% filter(Gene == g)
  pltp = pvalues2pubr(rna_pvalues, g, levels(vioplt[[ANNOTVAR]]), bracketx)
  pltp
  #pltp = pvalues2pubr(rna_pvalues, g, bracketx)
  maxy = max(pltdata$normalized_expression)
  print(ggplot(data = pltdata, aes(x = !!sym(ANNOTVAR), y = normalized_expression, fill = Disease_Status)) +
    geom_violin(color = "black", position = "dodge", scale = "width") +
    geom_boxplot(color = "black", position = position_dodge(0.9), alpha = 0.5, width=0.2, notch=TRUE) + 
    geom_bracket(data = pltp, aes(xmin = group1x, xmax = group2x, label = p.sym), y.position = maxy*1.2, step.increase=0.1, inherit.aes = FALSE) +
    scale_fill_manual(values=unname(DISEASESTATUSCOLORS)) +
    ylab(sprintf("Normalized Expression (%s)", g)) +
    theme_bw() +
    theme(
          axis.title.x = element_blank(),
          #axis.text.x = element_blank(),
          #axis.ticks.x = element_blank(),
          #axis.ticks.length.x = unit(0, "in"),
          axis.text.x = element_text(angle= 90),
          panel.grid = element_blank()
    ))
}
dev.off()
