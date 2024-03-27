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
library(Seurat)
library(grid)

# %% Main object
so = readRDS("../seuMergedPostHSP_forFigures_2023-09-17_09-03-10.rds")

################################################################################
################################################################################
# mLN
################################################################################
################################################################################

################################################################################
# %% Data prep
################################################################################

COMPAREVAR="Disease_Status"
ANNOTVAR="groupedAnnot"

so_mln_only = subset(so, Tissue == "MES")

clusters = unique(so_mln_only[["manualAnnot"]])
bcell_clusters = clusters[which(startsWith(clusters[,1], "B")),1]
nk_clusters = c("NK", "NK/ILC")

convert = rep("All NK Cells combined", length(nk_clusters))
names(convert)=nk_clusters
convert

so_mln_only = soAddGroupedAnnotVar(so_mln_only, "manualAnnot", "groupedAnnot", convert)

totalnk = subset(so_mln_only, !!sym(ANNOTVAR) == "All NK Cells combined")

################################################################################
# %% Violins
################################################################################

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

mln_plts=list()
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
    ggtitle(sprintf("mLN %s", g)) +
    theme_classic() +
    subplotTheme +
    theme(
          axis.title.x = element_blank(),
          plot.title = element_text(size=6, hjust=0.5),
          plot.title.position = "panel",
          panel.grid = element_blank(),
          axis.text.y = element_text(size=5),
          axis.title.y = element_text(size=5),
          plot.margin = unit(c(0,1,1,1), "mm"),
          axis.text.x = element_text(size = 5, angle = 45, vjust = 1, hjust = 1, color = "#000000")
    )
  mln_plts[[i]] = p
  i=i+1
}

################################################################################
################################################################################
# %% Spleen
################################################################################
################################################################################

################################################################################
# %% Data prep
################################################################################

COMPAREVAR="Disease_Status"
ANNOTVAR="groupedAnnot"

so_spleen_only = subset(so, Tissue == "Spleen")

clusters = unique(so_spleen_only[["manualAnnot"]])
bcell_clusters = clusters[which(startsWith(clusters[,1], "B")),1]
nk_clusters = c("NK", "NK/ILC")

convert = rep("All NK Cells combined", length(nk_clusters))
names(convert)=nk_clusters
convert

so_spleen_only = soAddGroupedAnnotVar(so_spleen_only, "manualAnnot", "groupedAnnot", convert)

totalnk = subset(so_spleen_only, !!sym(ANNOTVAR) == "All NK Cells combined")


################################################################################
# %% Violins
################################################################################

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

spleen_plts=list()
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
  ggtitle(sprintf("Spleen %s", g)) +
  theme_classic() +
  subplotTheme +
  theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 5, angle = 45, vjust = 1, hjust = 1, color = "#000000"),
        plot.title = element_text(size=6, hjust=0.5),
        plot.title.position = "panel",
        panel.grid = element_blank(),
        axis.text.y = element_text(size=5),
        axis.title.y = element_text(size=5),
        plot.margin = unit(c(0,1,1,1), "mm")
  )
  spleen_plts[[i]] = p
  i=i+1
}

# %%
pdf("/srv/http/betts/hpap/figures/s6_pieces.pdf", width = 4.25, height = 4.25)
mln_plts[[1]] + mln_plts[[2]] +
spleen_plts[[1]] + spleen_plts[[2]] + plot_layout(nrow = 2)
dev.off()
# %%

