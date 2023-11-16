library(ggplot2)
source("figures/genericFigureSettings.R")
source("scripts/parse_comparisons.R")

pValSymnum <- function(x, showNs = TRUE) {
  tmp <- sapply(x, function(y) {
                  if (is.na(y)) {
                    return(NA)
                  }

                  if (y < 0.001) {
                    return("***")
                  } else if (y >= 0.001 & y < 0.01) {
                    return("**")
                  } else if (y >= 0.01 & y < 0.05) {
                    return("*")
                  } else {
                    return(ifelse(showNs, "ns", ""))
                  }
})

  return(tmp)
}

so_pln_only = readRDS("rds/so_pln_only.rds")

clusters = unique(so_pln_only[["manualAnnot"]])
bcell_clusters = clusters[which(startsWith(clusters[,1], "B")),1]
nk_clusters = c("NK", "NK/ILC")

convert = c(rep("All B cells combined", length(bcell_clusters)), rep("All NK Cells combined", length(nk_clusters)))
names(convert)=c(bcell_clusters, nk_clusters)
convert

so_pln_only = soAddGroupedAnnotVar(so_pln_only, "manualAnnot", "groupedAnnot", convert)
colnames(so_pln_only[[]])

sig_genes = readRDS("rds/wuv_compres_rna_genelist_V1.rds")
goirb = read.table("NK_list_231025.csv", sep=",", header=TRUE)
goirb = goirb[goirb$group == "rank",]$gene
colnames(sig_genes)
sig_genes=sig_genes[sig_genes$gene %in% goirb & sig_genes$matchup == "T1D_vs_ND" & sig_genes$p_val_adj_all <= 0.05,]
dim(sig_genes)

colnames(sig_genes)

sig_genes$gene = factor(sig_genes$gene, levels=sig_genes$gene[order(sig_genes$avg_log2FC, decreasing=FALSE)])
sig_genes$pvalsymm = pValSymnum(sig_genes$p_val_adj_all)

pdf("/srv/http/betts/hpap/figures/nk_ranked_bar.pdf", width=8.5, height=8.5)
ggplot(data=sig_genes, aes(x=gene, y=avg_log2FC, fill=pvalsymm)) +
  geom_bar(stat="identity", color = "black") +
  scale_fill_manual(values = COLORS[["pval-heatmap"]]) +
  geom_text(data = subset(sig_genes, avg_log2FC >= 0), aes(label=gene, x=gene), y=0.009, size = 2.65, hjust=0)+
  geom_text(data = subset(sig_genes, avg_log2FC < 0), aes(label=gene, x=gene), y=-0.015, size = 2.65, hjust=1)+
  geom_hline(yintercept = 0, color = "black") +
  coord_flip()+
  ylab("log2(Fold Change)")+
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.length.y = unit(0,"cm"),
    legend.title = element_blank(),
    axis.line.x = element_line(color="black")
  )
dev.off()


range(-log10(sig_genes$p_val_adj_all))
