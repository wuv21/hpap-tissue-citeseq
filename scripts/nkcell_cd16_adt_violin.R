library(magrittr)
library(ggplot2)
library(patchwork)
library(dplyr)
source("scripts/so_helpers.R")
source("scripts/generic.R")

goi = tibble::as_tibble(read.table("miscellaneous_gene_lists/B_NK_CITEseq_gene_list_allComparisons.csv", sep = ',', header=TRUE, row.names=NULL))
goi$comparison = gsub("AAb", "AAb+", goi$comparison)
convert=goi[goi$gene %in% c("CD16", "CD56"),][,c("gene", "adt_number")]
convert=named_vector(convert$gene, convert$adt_number)
so_pln_only = readRDS("rds/so_pln_only.rds")

nk_cells = Cells(subset(so_pln_only, manualAnnot %in% c("NK", "NK/ILC")))
so_nk = so_pln_only[,nk_cells]
Idents(so_nk) = "Disease_Status"
DefaultAssay(so_nk) = "adt"

convert
toplt = GetAssayData(so_nk, "data")
toplt = as_tibble(cbind(so_nk[["Disease_Status"]][,1], toplt[which(rownames(toplt) == convert[["CD16"]]), ]))
colnames(toplt) = c("Disease_Status", "CD16")
toplt = toplt %>% 
  mutate(CD16 = as.numeric(CD16))
toplt

for (ds in unique(toplt$Disease_Status)) {
  print(ds)
  print(dim(toplt %>% filter(Disease_Status == ds & CD16 >= 1.2)))
  print(dim(toplt %>% filter(Disease_Status == ds & CD16 < 1.2)))

}

pdf("/srv/http/betts/hpap/figures/nkcells_cd16_violin.pdf", width = 11, height=8)
pltdata = toplt
print(ggplot(data = pltdata, aes(x = 1, y = CD16, fill = Disease_Status)) +
      geom_violin(color = "black", position = "dodge", scale = "width") +
      geom_boxplot(color = "black", position = position_dodge(0.9), alpha = 0.5, width=0.2, notch=TRUE) + 
      scale_fill_manual(values=unname(DISEASESTATUSCOLORS)) +
      geom_hline(yintercept = 1.2, linetype = "dashed", color = "black") +
      ylab("Normalized Expression (CD16)") +
      theme_bw() +
      theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.length.x = unit(0, "in"),
            #axis.text.x = element_text(angle= 90),
            panel.grid = element_blank()
            ))
dev.off()
as_tibble(toplt)

plts = list()
for (ds in names(subs)) {
  plt_d = subs[[ds]] %>% mutate_all(~ifelse(!!sym(convert[["CD16"]]) == 0 & !!sym(convert[["CD56"]]) == 0, NA, .))
  x_over = subs[[ds]] %>% filter(!!sym(convert[["CD16"]]) >= 1.2) %>% pull(!!sym(convert[["CD16"]])) %>% length
  x_under = subs[[ds]] %>% filter(!!sym(convert[["CD16"]]) < 1.2) %>% pull(!!sym(convert[["CD16"]])) %>% length
  y_over = subs[[ds]] %>% filter(!!sym(convert[["CD56"]]) >= 2.3) %>% pull(!!sym(convert[["CD56"]])) %>% length
  y_under = subs[[ds]] %>% filter(!!sym(convert[["CD56"]]) < 2.3) %>% pull(!!sym(convert[["CD56"]])) %>% length
  perc = sapply(c(x_over, x_under, y_over, y_under), function(x) round(x/dim(plt_d)[1], 4)*100)
  print(perc)
  plts[[ds]] = ggplot(plt_d, aes(x=!!sym(adtnums[2]), y=!!sym(adtnums[1]))) +
    stat_bin2d(bins=100) +
    scale_fill_distiller(palette = "Spectral") +
    geom_vline(xintercept = 1.2, linetype="dashed")+
    geom_text(y=3.5, x=1.45, label=sprintf("%s%%", perc[3]))+
    geom_text(y=3.5, x=0.95, label=sprintf("%s%%", perc[4]))+
    geom_hline(yintercept = 2.3, linetype="dashed")+
    geom_text(y=2.55, x=3.0, label=sprintf("%s%%", perc[1]))+
    geom_text(y=2.05, x=3.0, label=sprintf("%s%%", perc[2]))+
    xlab("CD16")+
    ylab("CD56")+
    ggtitle(ds)
}
ggsave(file="figures/CD16xCD56_scatter.tiff", plot = wrap_plots(plts, nrow=1), width=28, height=8.5, device="tiff", dpi=150)

Seurat::FeatureScatter(so, feature1="A0047", feature2="A0083", cells=nk_cells, group.by="Disease_Status", pt.size=0.5)
