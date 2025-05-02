# %%
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
source("scripts/deg.R")
library(Seurat)
library(grid)

#################################################################################
# %%
#################################################################################

FLOW_RDS_PATH="./figures/greg_flow_data/rds/dfLineageFilter.rds"
SC_RDS_PATH="/home/ubuntu/projmnt/betts/coculture/hpap/seuMergedPostHSP_forFigures_2025-01-12_04-07-24.rds"
set.seed(42)

#################################################################################
# %% A - B cell differences in NK cells
################################################################################
fd = readRDS(FLOW_RDS_PATH) %>%
  filter(LN_type == "pLN" & metric == "NK CD56dimCD16+") %>%
  mutate(`Disease Status` = factor(`Disease Status`, levels = c("ND", "AAb+", "T1D")))
figA =  ggplot(fd, aes(`Disease Status`, value, color = `Disease Status`)) +
  ylim(0,max(fd$value)*1.25)+
  labs(
    y = "% of NK cells",
    title = "pLN NK cells") +
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
        axis.title.x = element_blank(),
        plot.title = element_text(size = 6, hjust = 0.5),
        plot.title.position = "panel",
        plot.margin = unit(c(5,0,0,10), "pt"),
        axis.text.x = element_text(size = 5, angle = 45, vjust = 1, hjust = 1, color = "#000000"),
        axis.text.y = element_text(size = 5, color = "#000000"),
        axis.title.y = element_text(size = 5)
  )

#################################################################################
# %% B - Empty frame for flow traces
################################################################################
  
figB = ggplot() +
  ggtitle("CD56xCD16 NK cells flow")+
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size=6)
  )


################################################################################
# %% Data prep for C-G
################################################################################
COMPAREVAR="Disease_Status"
ANNOTVAR="groupedAnnot"

seu = readRDS(SC_RDS_PATH)
so_pln_only = subset(seu, Tissue %in% c("pLN-H", "pLN-T"))

clusters = unique(so_pln_only[["manualAnnot"]])
bcell_clusters = clusters[which(startsWith(clusters[,1], "B")),1]
nk_clusters = c("NK", "NK/ILC")

convert = c(rep("All B cells combined", length(bcell_clusters)), rep("All NK Cells combined", length(nk_clusters)))
names(convert)=c(bcell_clusters, nk_clusters)
convert

so_pln_only = soAddGroupedAnnotVar(so_pln_only, "manualAnnot", "groupedAnnot", convert)

totalnk = subset(so_pln_only, !!sym(ANNOTVAR) == "All NK Cells combined")

#################################################################################
# %% C - NK cells ranked bar plot
################################################################################
compsout_nk_rna = findMarkersCombinatorial(subset(so_pln_only, groupedAnnot == "All NK Cells combined"), combVar = "Disease_Status", assay = "RNA")

goirb = read.table("miscellaneous_gene_lists/NK_list_231025.csv", sep=",", header=TRUE)
goirb = goirb[goirb$group == "rank",]$gene
sig_genes = compsout_nk_rna[compsout_nk_rna$gene %in% goirb & compsout_nk_rna$matchup == "T1D_vs_ND" & compsout_nk_rna$p_val_adj_all <= 0.05,]

sig_genes$gene = factor(sig_genes$gene, levels=sig_genes$gene[order(sig_genes$avg_log2FC, decreasing=FALSE)])
sig_genes$pvalsymm = pValSymnum(sig_genes$p_val_adj_all)

sig_genes
compsout_nk_rna[compsout_nk_rna$gene %in% c("GZMB", "KLRB1"),]
fig6_de_stats = compsout_nk_rna[compsout_nk_rna$gene %in% c("GZMB", "KLRB1"),c("gene", "avg_log2FC", "p_val_adj_all", "pct.1", "pct.2", "upregulated", "matchup")]
fig6_de_stats
write.table(fig6_de_stats, file="/srv/http/betts/hpap/final_figures/amsesk/stats/fig6_de_stats.csv", sep=",", quote=F, row.names=F)

# %%
figC = ggplot(data=sig_genes, aes(x=gene, y=avg_log2FC, fill=pvalsymm)) +
  geom_bar(stat="identity", color = "black", lwd=0.15) +
  scale_fill_manual(values = COLORS[["pval-heatmap"]]) +
  geom_text(data = subset(sig_genes, avg_log2FC >= 0), aes(label=gene, x=gene, color=pvalsymm), y=0.005, size = 1.0, hjust=0, family="sans", fontface="bold")+
  geom_text(data = subset(sig_genes, avg_log2FC < 0), aes(label=gene, x=gene, color=pvalsymm), y=-0.005, size = 1.0, hjust=1, family="sans", fontface="bold")+
  scale_color_manual(values=c("ns" = "black", "*" = "black", "**" = "black", "***" = "white")) +
  geom_hline(yintercept = 0, color = "black") +
  guides(color="none") +
  coord_flip()+
  ylab("log2(Fold Change)")+
  theme_classic() +
  subplotTheme +
  theme(
        plot.background = element_blank(),
        panel.background = element_blank(),
        rect = element_blank(),
        legend.position = c(0.85,0.20),
        legend.text = element_text(size=4),
        legend.key.height=unit(1.5,"mm"),
        legend.key.width =unit(1.5,"mm"),
        plot.title = element_text(size=5, hjust=0.5),
        plot.title.position = "panel",
        plot.margin = margin(0.254,0,0.15,0, "cm"),
        panel.spacing = unit(0, "pt"),
        panel.grid = element_blank(),
        axis.text.x = element_text(size=4, color="black"),
        axis.title.x = element_text(family="sans", size=5, color="black",vjust=2),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length.y = unit(0,"cm"),
        legend.title = element_blank(),
        axis.line.y = element_blank()
  )
  ylims = c(min(sig_genes$avg_log2FC-0.01),max(sig_genes$avg_log2FC)+0.01)
  ylims
  figC = figC +
    scale_y_continuous(limits=ylims, expand=c(0,0))

#################################################################################
# %% D/E - NK cells violin/boxes
################################################################################

ANNOTVAR="groupedAnnot"

genes_of_interest = c("GZMB", "KLRB1")

pltData = totalnk[genes_of_interest,]

pltData = seuratObjMetaTibble(pltData, assay = "RNA")
attr(pltData, "datacol") = seq_along(genes_of_interest)+1
rna_meanExp = as.data.frame(mean_expression(pltData, COMPAREVAR, ANNOTVAR))
rna_stderr_bars = as.data.frame(expression_stderr_bars(pltData, COMPAREVAR, ANNOTVAR, method = mean, bar_method = sd))
rna_scaledExp = scale_expression(rna_meanExp, compareVar = COMPAREVAR, by = ANNOTVAR)
rna_scaledExp

nkvioplt = pltData %>% 
  select("cell", COMPAREVAR, ANNOTVAR, colnames(pltData)[attr(pltData, "datacol")]) %>%
  pivot_longer(cols = !c("cell", COMPAREVAR, ANNOTVAR), names_to = "Gene", values_to = "normalized_expression") %>%
  mutate(Disease_Status = factor(Disease_Status, levels = c("ND", "AAb+", "T1D")))

nkplts=list()
i=1
for (g in genes_of_interest) {
  pltdata =nkvioplt %>% filter(Gene == g)
  maxy = max(pltdata$normalized_expression)
  p = ggplot(data = pltdata, aes(x = Disease_Status, y = normalized_expression, fill = Disease_Status)) +
    geom_violin(color = "black", scale = "width", lwd=0.1) +
    geom_boxplot(color = "black", alpha = 0.5, width=0.2, notch=TRUE,  outlier.size=0.05, outlier.alpha=0.25, lwd=0.1) + 
    ylim(0,maxy*1.25) +
    scale_fill_manual(values=COLORS[["disease"]]) +
    ylab("Norm. Expr.") +
    guides(fill = "none") +
    ggtitle(sprintf("pLN %s", g)) +
    theme_classic() +
    subplotTheme +
    theme(
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(size=6, hjust=0.5),
          plot.title.position = "panel",
          panel.grid = element_blank(),
          axis.text.y = element_text(size=5),
          axis.title.y = element_text(size=5),
          plot.margin = unit(c(0,1,1,1), "mm")
    )
  if (i == 2) {
    p = p + theme(
      axis.text.x = element_text(size = 5, angle = 45, vjust = 1, hjust = 1, color = "#000000")
      )
  }
  nkplts[[i]] = p
  i = i+1
}

#################################################################################
# %% E/F - Heagmaps of differentially expressed genes between GZMB+/- NK cells
################################################################################
lod = 0.25
totalnk[["expresses_gzmb"]] = ifelse(FetchData(totalnk, "GZMB") >= lod, TRUE, FALSE)
DefaultAssay(totalnk) = "RNA"
Idents(totalnk) = "expresses_gzmb"
de_rna = FindMarkers(totalnk, ident.1 = TRUE, ident.2 = FALSE, assay = "RNA", logfc.threshold = 0.1)
de_rna[de_rna$p_val_adj <=0.05,]

de_rna = de_rna[de_rna$p_val_adj <=0.05,]
max_i = max(seq_along(rownames(de_rna)))
if (max_i > 20) { 
  max_i = 20
}
top20ish = de_rna[which(de_rna$avg_log2FC >= sort(de_rna$avg_log2FC, decreasing = T)[max_i] & de_rna$avg_log2FC > 0),]
top20ish
top20ish$dir = "up"
bottom20ish = de_rna[which(de_rna$avg_log2FC <= sort(de_rna$avg_log2FC, decreasing = F)[max_i] & de_rna$avg_log2FC < 0),]
bottom20ish$dir = "down"

# Stop if any of the gene names are shared between top and bottom 20ish
stopifnot(!any(rownames(top20ish) %in% rownames(bottom20ish)))

de_rna = rbind(top20ish, bottom20ish)
de_rna$pvalsymm = pValSymnum(de_rna$p_val_adj)

rna = seuratObjMetaTibble(totalnk, assay = "RNA")
rna[c("Disease_Status", "GZMB")] %>% group_by(Disease_Status) %>% summarize(GZMB = mean(GZMB))
attr(rna, "datacol") = seq_along(rownames(totalnk))+1
dim(rna)
rna_meanExp = as.data.frame(mean_expression(rna, compareVar = "expresses_gzmb", annotVar = "groupedAnnot"))

meanne = t(rna_meanExp[,c("expresses_gzmb", rownames(de_rna))])[-1,]
colnames(meanne) = c("GZMB-", "GZMB+")

col_fun_DiseaseHm_up <- circlize::colorRamp2(c(min(meanne[rownames(de_rna[de_rna$dir == "up",]),]), mean(meanne[rownames(de_rna[de_rna$dir == "up",]),]), max(meanne[rownames(de_rna[de_rna$dir == "up",]),])), c("blue", "white", "red"))
col_fun_DiseaseHm_down <- circlize::colorRamp2(c(min(meanne[rownames(de_rna[de_rna$dir == "down",]),]), mean(meanne[rownames(de_rna[de_rna$dir == "down",]),]), max(meanne[rownames(de_rna[de_rna$dir == "down",]),])), c("blue", "white", "red"))

de_rna_up = de_rna[de_rna$dir == "up",]
de_rna_down = de_rna[de_rna$dir == "down",]
figE = ComplexHeatmap::Heatmap(
                        meanne[rownames(de_rna[de_rna$dir == "up",]),],
                        na_col = "grey90", 
                        cluster_rows=F, 
                        cluster_columns=F, 
                        show_row_names=T,
                        column_names_rot=0,
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          if (j == 2) {
                            grid.text(de_rna_up[i,"pvalsymm"], x, y, gp = gpar(fontsize=4), vjust=0.75)
                          } #else {
                            #grid.text(round(meanne[i,j], 2), x, y, gp = gpar(fontsize=4), vjust=0.75)
                          #}
                        },

                        rect_gp = gpar(col = "white", lwd = 1), 
                        col = col_fun_DiseaseHm_up,
                        row_title_gp = gpar(fill = "black", col = "white", border = "black"),
                        row_names_side='left',
                        row_names_gp = gpar(fontsize=4),
                        column_names_gp = gpar(fontsize=4),
                        column_names_centered=TRUE,
                        show_heatmap_legend = TRUE,
                        heatmap_legend_param=list(title = "Norm. Expr.", 
                                                  title_gp = gpar(fontsize=3), 
                                                  title_gap = unit(1, "mm"),
                                                  labels_gp = gpar(fontsize=3), 
                                                  title_position="topcenter", 
                                                  legend_height=unit(0.3, "in"),
                                                  legend_width=unit(0.01, "in"),
                                                  grid_width=unit(1,"mm"),
                                                  grid_height=unit(1,"mm")
                        ),
                        border = FALSE,
                        width=unit(0.8, "in"),
                        )
figE = grid.grabExpr(draw(figE, column_title = "", column_title_gp = gpar(fontsize = 5), padding = unit(c(0, 0, 0, 0), "pt"), gap=unit(0,"mm")))
figF = ComplexHeatmap::Heatmap(
                        meanne[rownames(de_rna[de_rna$dir == "down",]),],
                        na_col = "grey90", 
                        cluster_rows=F, 
                        cluster_columns=F, 
                        show_row_names=T,
                        column_names_rot=0,
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          if (j == 1) {
                            grid.text(de_rna_down[i,"pvalsymm"], x, y, gp = gpar(fontsize=4), vjust=0.75)
                          }
                        },
                        rect_gp = gpar(col = "white", lwd = 1), 
                        col = col_fun_DiseaseHm_down,
                        row_title_gp = gpar(fill = "black", col = "white", border = "black"),
                        row_names_side='left',
                        row_names_gp = gpar(fontsize=3.5),
                        column_names_gp = gpar(fontsize=3.5),
                        column_names_centered=TRUE,
                        show_heatmap_legend = TRUE,
                        column_order = c(2,1),
                        heatmap_legend_param=list(title = "Norm. Expr.", 
                                                  title_gp = gpar(fontsize=3), 
                                                  labels_gp = gpar(fontsize=3), 
                                                  title_position="topcenter", 
                                                  legend_height=unit(0.3, "in"),
                                                  legend_width=unit(0.01, "in"),
                                                  grid_width=unit(1,"mm"),
                                                  grid_height=unit(1,"mm")
                        ),
                        border = FALSE,
                        width=unit(0.8, "in"),
)
figF = grid.grabExpr(draw(figF, column_title = "", column_title_gp = gpar(fontsize = 5), padding = unit(c(0, 0, 0, 0), "pt"), gap=unit(0,"mm")))

# %%
fig6layout <- c(
  patchwork::area(1,1,30,12), #a
  patchwork::area(1,13,30,30), #b
  patchwork::area(31,1,90,15), #c
  patchwork::area(31,16,60,30), #d
  patchwork::area(61,16,90,30), #e
  patchwork::area(91,1,150,14), #f
  patchwork::area(91,15,150,30) #g
  )
plot = wrap_elements(full=figA) + figB + 
  wrap_elements(full=figC) + nkplts[[1]] + nkplts[[2]] +
  wrap_elements(full=figE) + wrap_elements(full=figF) +
  patchwork::plot_layout(design=fig6layout) +
  patchwork::plot_annotation(tag_levels = list(letters[1:7])) &
  plotTagTheme

pdf("/srv/http/betts/hpap/final_figures/amsesk/pdf/fig6_v3_final.pdf", width=3.75, height=5.50, family="sans")
print(plot)
dev.off()

ylims
# %%
# saveFinalFigure(plot=plot,
#                 # prefixDir = "figures/outs",
#                 prefixDir = "/srv/http/betts/hpap/final_figures",
#                 fn = "fig6_v2",
#                 devices = c("Cairo"),
#                 addTimestamp = FALSE,
#                 gheight=5.50,
#                 gwidth=3.75)

# %%
pdf("/srv/http/betts/hpap/final_figures/dotty.pdf", width=8.5, height=8.5)
p=Seurat::DotPlot(seu, features = c("NR4A1", "NR4A2", "BCL2", "FOS", "JUN"), assay="RNA", group.by = "manualAnnot")
print(p)
dev.off()

# %%

