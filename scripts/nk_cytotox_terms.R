library(fgsea)
library(msigdbr)
library(enricheR)

sig_genes= readRDS("rds/wuv_compres_rna_genelist_V1.rds")
ctg = read.table("NK_list_231025.csv", sep=",", header=TRUE)
ctg = ctg[ctg$group == "GO",]$gene
sig_genes=sig_genes[sig_genes$gene %in% ctg & sig_genes$matchup == "T1D_vs_ND" & sig_genes$p_val_adj_all >= 0.05,]
dim(sig_genes)

msigdbr_df = msigdbr::msigdbr()
sub = msigdbr_df[msigdbr_df$gene_symbol %in% sig_genes$gene,]
msigdbr_list = split(x = sub$gene_symbol, f = sub$gs_name)
fgsea(pathways = msigdbr_list, stats=)


install.packages("enrichR")
