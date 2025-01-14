library(magrittr)
library(ggplot2)
library(languageserver)
library(dplyr)
library(tidyr)
source("scripts/dimPlots.R")
source("scripts/so_helpers.R")
source("scripts/parse_comparisons.R")
source("figures/genericFigureSettings.R")
source("scripts/transform_expression.R")
library(Seurat)
library(grid)

#################################################################################
# %%
#################################################################################

FLOW_RDS_PATH="./figures/greg_flow_data/rds/dfLineageFilter.rds"
SC_RDS_PATH="/home/ubuntu/projmnt/betts/hpap/rds/seuMergedPostHSP_forFigures_2025-01-12_04-07-24.rds"

#################################################################################
# %% A - B cell differences in pLN
################################################################################
fd = readRDS("./figures/greg_flow_data/rds/dfLineageFilter.rds") %>%
  filter(LN_type == "pLN" & metric == "B cells CD27+") %>%
  mutate(`Disease Status` = factor(`Disease Status`, levels = c("ND", "AAb+", "T1D")))
figA = ggplot(fd, aes(`Disease Status`, value, color = `Disease Status`)) +
  ylim(0,max(fd$value)*1.25)+
  labs(
    y = "% of B cells",
    title = "pLN B cells") +
  geom_boxplot(
    width = 0.8,
    outlier.shape = NA,
    show.legend = FALSE) +
  geom_point(
    size = 1,
    stroke = 0.2,
    alpha = 0.4,
    show.legend = FALSE,
    position = position_jitterdodge(jitter.width = 1, dodge.width = 0.8)) +
  scale_color_manual(values = COLORS$disease) +
  theme_classic() +
  subplotTheme +
  theme(
    plot.margin = unit(c(5,0,0,8), "pt"),
    axis.title.x = element_blank(),
    plot.title = element_text(size = 6, hjust = 0.5),
    plot.title.position = "panel",
    axis.text.x = element_text(size = 5, angle = 45, vjust = 1, hjust = 1, color = "#000000"),
    axis.text.y = element_text(size = 5, color = "#000000"),
    axis.title.y = element_text(size = 5))

#################################################################################
# %% B - Empty frame for flow traces
################################################################################

figB = ggplot() +
  ggtitle("CD27xCD69 B cell flow traces")+
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size=6)
  )


#################################################################################
# %% Data prep for C-F
################################################################################

seu = readRDS(SC_RDS_PATH)
so_pln_only = subset(seu, Tissue %in% c("pLN-H", "pLN-T"))

clusters = unique(so_pln_only[["manualAnnot"]])
bcell_clusters = clusters[which(startsWith(clusters[,1], "B")),1]

convert = named_vector(
                       bcell_clusters,
                       c(rep("All B cells combined", length(bcell_clusters)))
                       )
convert

so_pln_only = soAddGroupedAnnotVar(so_pln_only, "manualAnnot", "groupedAnnot", convert)

col_fun_DiseaseHm <- circlize::colorRamp2(c(-0.75, 0, 0.75), c("blue", "white", "red"))

COMPAREVAR="Disease_Status"

goi = tibble::as_tibble(read.table("miscellaneous_gene_lists/HPAP_CITEseq_gene_list_V3.csv", sep = ',', header=TRUE, row.names=NULL, fill = NA))
coi = unique(unname(convert))
Seurat::DefaultAssay(so_pln_only) = "RNA"
igh_rna = goi[goi$modality == "RNA" & grepl("^IGH", goi$gene),]$gene

#################################################################################
# %% DEG for B cells
#################################################################################

# All B cells together, v3
compsout_rna = findMarkersCombinatorial(subset(so_pln_only, groupedAnnot == "All B cells combined"), combVar = "Disease_Status", assay = "RNA")
compsout_rna = compsout_rna[compsout_rna$gene %in% goi[goi$Heatmap == 2,]$gene,]
saveRDS(compsout_rna, "rds/wuv_compres_allbcells_rna_genelist_V3.rds")

# Bcells v3 separate clusters
annotVar = "manualAnnot"
compsout_rna = lapply(unique(bcell_clusters), function(xx) {
                        df = findMarkersCombinatorial(subset(so_pln_only, !!sym(annotVar) == xx), combVar = "Disease_Status", assay = "RNA")
                        attr(df, "cluster") = xx
                        print(df)
})
names(compsout_rna) = bcell_clusters
names(compsout_rna)
unique(so_pln_only[["manualAnnot"]][,1])

saveRDS(compsout_rna, "rds/wuv_compres_bcellsep_rna_genelist_V3.rds")

#################################################################################
# %% Some figure variables
#################################################################################

HMLABSIZE=3
HMEXPWID=0.5
HMPVALWID=0.12

#################################################################################
# %% Ci - Total B cells heavy chain heatmap
################################################################################

ANNOTVAR = "groupedAnnot"

compres = readRDS("./rds/wuv_compres_allbcells_rna_genelist_V3.rds")
totalb_compres = list()
totalb_compres[["All B cells combined"]] = compres
write.table(totalb_compres[["All B cells combined"]] %>% filter(gene %in% c("IGHA1", "IGHM", "IGHG1")), "./outs/all_b_cells_IGs_figSY_stats.csv", sep = ",", row.names = FALSE, quote = FALSE)

hmplt = so_pln_only[igh_rna,]
hmplt = subset(hmplt, !!sym(ANNOTVAR) %in% coi)

hmplt_data = seuratObjMetaTibble(hmplt, assay = "RNA")
attr(hmplt_data, "datacol") = seq_along(igh_rna)+1
rna_meanExp = as.data.frame(mean_expression(hmplt_data, COMPAREVAR, ANNOTVAR))
rna_scaledExp = scale_expression(rna_meanExp, compareVar = COMPAREVAR, by = ANNOTVAR)
totalb_compres[["All B cells combined"]]$p_val_adj_all_sym = pValSymnum(totalb_compres[["All B cells combined"]]$p_val_adj_all)
rna_pvalues = lapply(rna_scaledExp, fillPvalsFromFindFeaturesCombinatorial, totalb_compres, modality="RNA", pvalue_column = "p_val_adj_all_sym", na_fill_value = "ns")

colnames(rna_scaledExp[[1]]) = gsub('[ ].*', '', colnames(rna_scaledExp[[1]]))
colnames(rna_scaledExp[[1]]) = gsub('[()]', '', colnames(rna_scaledExp[[1]]))

hm_rna_scaled= ComplexHeatmap::Heatmap(
                                       rna_scaledExp[[1]],
                                       na_col = "grey90",
                                       cluster_rows=T,
                                       cluster_columns=F,
                                       show_row_names=T,
                                       rect_gp = gpar(col = "white", lwd = 1),
                                       col = col_fun_DiseaseHm,
                                       row_title_gp = gpar(fill = "black", col = "white", border = "black"),
                                       row_names_side='left',
                                       row_names_gp = gpar(fontsize=HMLABSIZE),
                                       column_names_gp = gpar(fontsize=HMLABSIZE),
                                       show_heatmap_legend = FALSE,
                                       border = FALSE,
                                       width=unit(HMEXPWID, "in"),
                                       column_order = c(2,1,3)
)
hm_rna_pval = ComplexHeatmap::Heatmap(
                                      rna_pvalues[[1]],
                                      na_col = "grey90",
                                      cluster_rows=T,
                                      cluster_columns=F,
                                      show_row_names=F,
                                      rect_gp = gpar(col = "white", lwd = 1),
                                      col = COLORS[["pval-heatmap"]],
                                      row_title_gp = gpar(fill = "black", col = "white", border = "black"),
                                      row_names_gp = gpar(fontsize=HMLABSIZE),
                                      column_names_gp = gpar(fontsize=HMLABSIZE),
                                      show_heatmap_legend = FALSE,
                                      border = FALSE,
                                      width=unit(HMPVALWID, "in"),
                                      column_order=c(1,3,2)
)
figCi = grid.grabExpr(draw(hm_rna_scaled + hm_rna_pval, row_dend_side='left', row_dend_width=unit(1, "mm"), gap = unit(0.5, "mm"), column_title = "Total B Cells", column_title_gp = gpar(fontsize = 4), padding = unit(c(0, 10, 0, 0), "pt")))

#################################################################################
# %% Cii/Ciii- Class Switched Memory #1 / IgM+ Marginal Zone-like #2 heavy chain heatmaps
################################################################################

ANNOTVAR = "manualAnnot"

compres = readRDS("./rds/wuv_compres_bcellsep_rna_genelist_V3.rds")
names(compres)
write.table(compres[["B class-switched memory #1"]] %>% filter(gene %in% c("IGHA1", "IGHM", "IGHG1")), "./outs/B_class_switched_memory_1_IGs_figSY_stats.csv", sep = ",", row.names = FALSE, quote = FALSE)
write.table(compres[["B IgM+ memory/marginal zone like"]] %>% filter(gene %in% c("IGHA1", "IGHM", "IGHG1")), "./outs/B_igm_memory_marginal_zone_like_IGs_figSY_stats.csv", sep = ",", row.names = FALSE, quote = FALSE)

hmplt = subset(so_pln_only[igh_rna,], !!sym(ANNOTVAR) %in% bcell_clusters)

hmplt_data = seuratObjMetaTibble(hmplt, assay = "RNA")
attr(hmplt_data, "datacol") = seq_along(igh_rna)+1
rna_meanExp = as.data.frame(mean_expression(hmplt_data, COMPAREVAR, ANNOTVAR))
rna_scaledExp = scale_expression(rna_meanExp, compareVar = COMPAREVAR, by = ANNOTVAR)
names(rna_scaledExp) = sapply(rna_scaledExp, function(x) attr(x, "cluster"))
compres = lapply(compres, function(x) {
                           x$p_val_adj_all_sym = pValSymnum(x$p_val_adj_all)
                           x
})
rna_pvalues = lapply(rna_scaledExp, fillPvalsFromFindFeaturesCombinatorial, compres, modality="RNA", pvalue_column = "p_val_adj_all_sym", na_fill_value = "ns")
names(rna_pvalues) = names(rna_scaledExp)

hms=list()
i=1
for (clust in c("B class-switched memory #1", "B IgM+ memory/marginal zone like")) {
  colnames(rna_scaledExp[[clust]]) = gsub('[ ].*', '', colnames(rna_scaledExp[[clust]]))
  colnames(rna_scaledExp[[clust]]) = gsub('[()]', '', colnames(rna_scaledExp[[clust]]))

  hm_rna_scaled= ComplexHeatmap::Heatmap(
                                         rna_scaledExp[[clust]],
                                         na_col = "grey90",
                                         cluster_rows=T,
                                         cluster_columns=F,
                                         show_row_names=T,
                                         rect_gp = gpar(col = "white", lwd = 1),
                                         col = col_fun_DiseaseHm,
                                         row_title_gp = gpar(fill = "black", col = "white", border = "black"),
                                         row_names_side='left',
                                         row_names_gp = gpar(fontsize=HMLABSIZE),
                                         column_names_gp = gpar(fontsize=HMLABSIZE),
                                         show_heatmap_legend = FALSE,
                                         border = FALSE,
                                         width=unit(HMEXPWID, "in"),
                                         column_order = c(2,1,3)
  )
  hm_rna_pval = ComplexHeatmap::Heatmap(
                                        rna_pvalues[[clust]],
                                        na_col = "grey90",
                                        cluster_rows=T,
                                        cluster_columns=F,
                                        show_row_names=F,
                                        rect_gp = gpar(col = "white", lwd = 1),
                                        col = COLORS[["pval-heatmap"]],
                                        row_title_gp = gpar(fill = "black", col = "white", border = "black"),
                                        row_names_gp = gpar(fontsize=HMLABSIZE),
                                        column_names_gp = gpar(fontsize=HMLABSIZE),
                                        show_heatmap_legend = FALSE,
                                        border = FALSE,
                                        width=unit(HMPVALWID, "in"),
                                        column_order=c(1,3,2)
  )
  hms[[i]] = grid.grabExpr(draw(hm_rna_scaled + hm_rna_pval, row_dend_side='left', row_dend_width=unit(1, "mm"), gap = unit(0.5, "mm"), column_title = clust, column_title_gp = gpar(fontsize = 4), padding = unit(c(0, 0, 0, 0), "pt")))
  i=i+1
}

figCii = hms[[1]]
figCiii = hms[[2]]

lgds = list("expr" = ComplexHeatmap::Legend(
                                            col_fun=circlize::colorRamp2(c(-0.75, 0, 0.75), c("blue", "white", "red")),
                                            direction = "horizontal",
                                            title_position = "topcenter",
                                            title = "Avg. Expression",
                                            title_gp = gpar(fontsize = 4, fontface = "plain", hjust = 0.0),
                                            title_gap = unit(0.75, "mm"),
                                            labels_gp = gpar(fontsize = 4),
                                            grid_height = unit(1, "mm"),
                                            legend_width = unit(0.625, "in"),
                                            legend_height = unit(1, "mm")
                                            ),
            "pvalues" = ComplexHeatmap::Legend(labels = names(COLORS[["pval-heatmap"]]),
                                               title = "Significance",
                                               title_position = "topcenter",
                                              title_gap = unit(0.5, "mm"),
                                               direction="horizontal",
                                               legend_gp = gpar(fill = unname(COLORS[["pval-heatmap"]])),
                                              grid_height = unit(1, "mm"),
                                              grid_width = unit(2, "mm"),
                                                legend_width = unit(0.625, "in"),
                                                legend_height = unit(0.75, "mm"),
                                               title_gp = gpar(fontsize=4),
                                               labels_gp = gpar(fontsize=4, fontface = "plain"),
                                               nrow=1
                                               )
            )


hm_lgd = grid.grabExpr(draw(lgds[["expr"]], x = unit(0.30, "npc"), y = unit(0, "npc"), just='bottom') +
                       draw(lgds[["pvalues"]], x = unit(0.70, "npc"), y = unit(0.15, "npc"), just='bottom')
                     )

#################################################################################
# %% D/E/Fi - Total B cells heavy chain violin/box plots
################################################################################
ANNOTVAR="groupedAnnot"

genes_of_interest = c("IGHA1", "IGHM", "IGHG1")
plttotalb = subset(so_pln_only[genes_of_interest,], groupedAnnot == "All B cells combined")
plttotalb = seuratObjMetaTibble(plttotalb, assay = "RNA")
attr(plttotalb, "datacol") = seq_along(genes_of_interest)+1
rna_meanExp = as.data.frame(mean_expression(plttotalb, COMPAREVAR, ANNOTVAR))
rna_scaledExp = scale_expression(rna_meanExp, compareVar = COMPAREVAR, by = ANNOTVAR)

totalb_vioplt = plttotalb %>%
  select("cell", COMPAREVAR, ANNOTVAR, colnames(plttotalb)[attr(plttotalb, "datacol")]) %>%
  pivot_longer(cols = !c("cell", COMPAREVAR, ANNOTVAR), names_to = "Gene", values_to = "normalized_expression") %>%
  mutate(Disease_Status = factor(Disease_Status, levels = c("ND", "AAb+", "T1D")))

totalb_plots = list()
i = 1
for (g in genes_of_interest) {
  pltdata = totalb_vioplt %>% filter(Gene == g)
  maxy = max(pltdata$normalized_expression)
   p = ggplot(data = pltdata, aes(x = Disease_Status, y = normalized_expression, fill = Disease_Status)) +
     geom_violin(color = "black", scale = "width", lwd=0.1) +
    geom_boxplot(color = "black", alpha = 0.5, width=0.2, notch=TRUE,  outlier.size=0.05, outlier.alpha=0.25, lwd=0.1) +
    scale_y_continuous(labels=scales::label_number(accuracy=0.1)) +
    scale_fill_manual(values=COLORS[["disease"]]) +
    ylab(sprintf("%s\nNorm. Expr.",g)) +
    guides(fill = "none") +
    theme_classic() +
    subplotTheme +
    theme(
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(size=4, hjust=0.5, vjust=5),
          plot.title.position = "panel",
          panel.grid = element_blank(),
          axis.text.y = element_text(size=5),
          axis.title.y = element_text(size=5),
    )
    topm=0
    botm=0
    if (i == 1) {
      p = p +
        ggtitle("Total B Cells")
    }
    if (i == 3) {
      p = p + theme(
        axis.text.x = element_text(size = 5, angle = 45, vjust = 1, hjust = 1, color = "#000000")
        )
    } else {
      botm=10
    }
    totalb_plots[[i]] = p+
     theme(
           plot.margin = unit(c(0,0,0,3), "pt")
      )

    i = i+1
}

#################################################################################
# %% D/E/Fii - B class-switched memory #1 and B IgM+ memory/marginal zone #2 cells heavy chain violin/box plots
################################################################################
ANNOTVAR="manualAnnot"

genes_of_interest = c("IGHA1", "IGHM", "IGHG1")
pltData = so_pln_only[genes_of_interest,]

pltsepb = subset(pltData,  !!sym(ANNOTVAR) %in% bcell_clusters)

pltsepb = seuratObjMetaTibble(pltsepb, assay = "RNA")
attr(pltsepb, "datacol") = seq_along(genes_of_interest)+1
rna_pctexp = percent_expressing(pltsepb, COMPAREVAR, ANNOTVAR, zero = 0.0)
rna_meanExp = as.data.frame(mean_expression(pltsepb, COMPAREVAR, ANNOTVAR))
rna_scaledExp = scale_expression(rna_meanExp, compareVar = COMPAREVAR, by = ANNOTVAR)
names(rna_scaledExp) = sapply(rna_scaledExp, function(x) attr(x, "cluster"))

sepb_vioplt = pltsepb %>%
  select("cell", COMPAREVAR, ANNOTVAR, colnames(pltsepb)[attr(pltsepb, "datacol")]) %>%
  pivot_longer(cols = !c("cell", COMPAREVAR, ANNOTVAR), names_to = "Gene", values_to = "normalized_expression") %>%
  mutate(Disease_Status = factor(Disease_Status, levels = c("ND", "AAb+", "T1D")))

sepb_plots = list()
for (clust in c("B class-switched memory #1", "B IgM+ memory/marginal zone like")) {
  sepb_plots[[clust]] = list()
  i = 1
  for (g in genes_of_interest) {
    pltdata = sepb_vioplt %>% filter(manualAnnot == clust & Gene == g)
    maxy = max(pltdata$normalized_expression)
    p = ggplot(data = pltdata, aes(x = Disease_Status, y = normalized_expression, fill = Disease_Status)) +
      geom_violin(color = "black", scale = "width", lwd=0.1) +
      geom_boxplot(color = "black", alpha = 0.5, width=0.2, notch=TRUE, outlier.size=0.01, outlier.alpha=0.25, outlier.shape=21, lwd=0.1) +
      scale_fill_manual(values=COLORS[["disease"]]) +
      scale_y_continuous(labels=scales::label_number(accuracy=0.1)) +
      ylab(sprintf("%s\nNorm. Expr.",g)) +
      guides(fill = "none") +
      theme_classic() +
      subplotTheme +
      theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            plot.title = element_text(size=4, hjust=0.5, color = "black", vjust=5),
            plot.title.position = "panel",
            panel.grid = element_blank(),
            axis.text.y = element_text(size=5, color="black"),
            axis.title.y = element_blank(),
      )
    topm=0
    botm=0
    if (i == 1) {
      p = p +
        ggtitle(clust)
    } else {}

    if (i == 3) {
      p = p + theme(
        axis.text.x = element_text(size = 5, angle = 45, vjust = 1, hjust = 1, color = "#000000")
        )
    } else {
      botm = 10
    }
   sepb_plots[[clust]][[i]] = p +
     theme(
      plot.margin = unit(c(0,0,0,3), "pt")
      )
    i = i+1
  }
}

# %%
figSYlayout <- c(
  patchwork::area(1,1,7,2), # a
  patchwork::area(1,3,7,6), # b
  patchwork::area(8,1,17,6), # c
  patchwork::area(18,1,21,6), # d
  patchwork::area(22,1,26,6), # e
  patchwork::area(27,1,29,6) # f
)

plot = wrap_elements(full=figA) + wrap_elements(plot=figB) +
wrap_elements(full=wrap_elements(full=wrap_plots(list(figCi,figCii,figCiii))) +
                              wrap_elements(full=hm_lgd)+
                                plot_layout(design=c(
                                                  area(1,1,38,1),
                                                  area(39,1,40,1)
                                                  )
                              )
                              ) +
  wrap_plots(list(totalb_plots[[1]], sepb_plots[["B class-switched memory #1"]][[1]], sepb_plots[["B IgM+ memory/marginal zone like"]][[1]])) +
  wrap_plots(list(totalb_plots[[2]], sepb_plots[["B class-switched memory #1"]][[2]], sepb_plots[["B IgM+ memory/marginal zone like"]][[2]])) +
  wrap_plots(list(totalb_plots[[3]], sepb_plots[["B class-switched memory #1"]][[3]], sepb_plots[["B IgM+ memory/marginal zone like"]][[3]])) +
  patchwork::plot_layout(design=figSYlayout) +
  patchwork::plot_annotation(tag_levels = list(c(LETTERS[1:3], LETTERS[4], rep("",2), LETTERS[5], rep("",2), LETTERS[6], rep("",2))))

saveFinalFigure(plot=plot,
                prefixDir = "/srv/http/betts/hpap/final_figures/",
                fn = "figSY_v1",
                devices = c("pdf_base"),
                addTimestamp = FALSE,
                gheight=5.50,
                gwidth=3.75)
# %%

