
library(ggplot2)

so_pln_only = readRDS("rds/so_pln_only.rds")


clusters = unique(so[["manualAnnot"]])
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
sig_genes=sig_genes[sig_genes$gene %in% goirb & sig_genes$matchup == "T1D_vs_ND", & sig_genes$p_val_adj_all >= 5]
dim(sig_genes)

colnames(sig_genes)
sig_genes
colscale = hcl.colors(n=21, palette = "Earth")[c(1,6,11,16,21)]
sig_genes$gene = factor(sig_genes$gene, levels=sig_genes$gene[order(sig_genes$avg_log2FC, decreasing=FALSE)])
pdf("/srv/http/betts/hpap/figures/nk_ranked_bar_e-10.pdf", width=8.5, height=16.5)
ggplot(data=sig_genes, aes(x=gene, y=avg_log2FC, fill=-log10(p_val_adj_all))) +
  geom_bar(stat="identity", color = "black") +
  scale_fill_gradient2(low = colscale[1], mid = colscale[3], high = colscale[5], midpoint = 11) +
  geom_text(data = subset(sig_genes, avg_log2FC >= 0), aes(label=gene, x=gene), y=0.009, size = 2.65, hjust=0)+
  geom_text(data = subset(sig_genes, avg_log2FC < 0), aes(label=gene, x=gene), y=-0.015, size = 2.65, hjust=1)+
  geom_hline(yintercept = 0, color = "black") +
  coord_flip()+
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.length.y = unit(0,"cm"),
    axis.line.x = element_line(color="black")
  )
dev.off()


range(-log10(sig_genes$p_val_adj_all))
