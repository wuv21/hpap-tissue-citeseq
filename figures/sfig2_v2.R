# note that this code is written to be run from the project base directory
renv::load("/data/hpap-citeseq/hpap-citeseq-analysis")

source("figures/genericFigureSettings.R")
source("scripts/dimPlots.R")
source("scripts/deg.R")
library(Seurat)
library(cowplot)
library(forcats)
library(presto)
library(rstatix) #dplyr-compatible stats package
library(ggpubr) #adds p-values annotations to ggplots
library(WRS2) #contains robust statistical tests
library(multcomp) #allows for multiple comparisons within ANCOVA
library(ggrepel)

set.seed(42)


################################################################################
# S2C - hsp filter cutoff
################################################################################
load("rds/preManualAnnotSubset.rData")

seuMerged$Disease_Status <- factor(seuMerged$Disease_Status, levels = c("ND", "AAb+", "T1D"))

seuFiltDf <- FetchData(
  object = seuMerged,
  vars = c("TissueCondensed", "heatShockProgram1", "Disease_Status", "DonorID")) %>%
  mutate(Disease_Status = factor(Disease_Status, levels = c("ND", "AAb+", "T1D")))

hspCutoff <- quantile(seuMerged$heatShockProgram1, probs = 0.95)
p2 <- seuFiltDf %>%
  ggplot(aes(x = heatShockProgram1, color = TissueCondensed)) +
  geom_density() +
  theme_classic() +
  geom_vline(xintercept = hspCutoff, color = "black", linetype = "dotted") +
  labs(x = "HSP score", y = "Density", color = "Tissue") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  subplotTheme +
  theme(
    legend.key.size = unit(3, "pt"),
    legend.text = element_text(size = BASEPTFONTSIZE),
    legend.position = "bottom",
    legend.title = element_text(size = BASEPTFONTSIZE),
    legend.margin = margin(t = -5, b = 0))

rm(seuMerged)
gc()

source("figures/genericFigureSettings.R")
source("scripts/dimPlots.R")
source("scripts/deg.R")

################################################################################
# S2A - researcher differences
# modified code from @Greg
################################################################################
dfFresh <- readRDS("figures/greg_flow_data/rds/dfFresh.rds")

dfScientistStat <- dfFresh %>%
  filter(Staining_Flow_Researcher == "D" | Staining_Flow_Researcher == "A") %>% #these researchers did the vast majority of the data acquisition, and are the only ones considered for this comparison
  filter(LN_type == "pLN") %>%
  filter(`Disease Status` == "ND") %>%
  filter(!str_detect(metric, "Tn CD127.CD27.")) %>% #all filtered values are 0, stats test would fail
  filter(!str_detect(metric, "Tn CD127.")) %>% #all filtered values are 100, stats test would fail
  filter(!str_detect(metric, "Tn CD27.")) %>% #all filtered values are 100, stats test would fail
  group_by(metric) %>%
  wilcox_test(value ~ Staining_Flow_Researcher) %>%
  adjust_pvalue() %>%
  filter(p.adj < 0.05) %>%
  mutate(pvalSymbol = case_when(p.adj < 0.001 ~ "***",
    p.adj >= 0.001 & p.adj < 0.01 ~ "**",
    p.adj >= 0.01 & p.adj < 0.05 ~ "*")) %>%
  arrange(p.adj)

figA <- dfFresh %>%
  filter(Staining_Flow_Researcher != "C") %>%
  filter(metric %in% dfScientistStat$metric) %>%
  filter(`Disease Status` == "ND") %>%
  mutate(metric = fct_relevel(`metric`, dfScientistStat$metric)) %>%
  ggplot(aes(Staining_Flow_Researcher, value, color = Staining_Flow_Researcher)) +
  geom_boxplot(outlier.size = 0.3) +
  stat_pvalue_manual(data = dfScientistStat, label = "pvalSymbol", 
    y.position = 105, xmin = "group1", xmax = "group2",
    tip.length = 0, label.size = 2) +
  ylab("% of Parent") +
  guides(color = guide_legend(title = "Researcher")) +
  theme_classic() +
  scale_color_brewer(palette = "Dark2") +
  subplotTheme +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
  ) +
  facet_wrap(~ factor(metric), nrow = 1, strip.position = "bottom") +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = BASEPTFONTSIZE),
    legend.margin = margin(t = -15, b = 0),
    strip.text.x.bottom = element_text(angle = 90, hjust = 1, size = 6),
    strip.background = element_rect(fill = "#FFFFFF00", color = "#FFFFFF00"),
    strip.placement = "outside")


################################################################################
# S2B - hsp score
# 
# note that this is a modified approach from my main_merge.R where i first find
# the hsp signatures from certain samples. this was used to inform the
# approach used in main_merge.R to filter out based on HSP signatures.
################################################################################
# unhashed sample run info
unhashedSampleMeta <- read.csv("metadata/unhashed_sample_metadata.csv")

seuListFns <- rdsCheckAndRun(
  fn = "rds/allSeuList_v1.rds",
  f = initialSeuratListMaker,
  version = "Merge",
  sampleIds = unhashedSampleMeta$id,
  sampleRuns = unhashedSampleMeta$run,
  rnaOutDir = rnaOut,
  adtOutDir = adtOut,
  htoOutDir = htoOut,
  plotDir = "outs/qc")

seuList <- lapply(seuListFns, function(x) {
  return(readRDS(x))
})

names(seuList) <- unhashedSampleMeta$id

# remove the following wells because poor QC
seuList[["20220531_hpap_2_Pool_Run1_Well2"]] <- NULL
seuList[["20220531_hpap_2_Pool_Run2_Well14"]] <- NULL

sampleMeta <- read.csv("metadata/sample_metadata.csv")
sampleMetaDict <- sampleMeta %>%
  group_split(Run)

names(sampleMetaDict) <- sapply(sampleMetaDict, function(x) {return(as.character(x[1, "Run"]))})

for (i in seq(1, length(seuList))) {
  xName <- names(seuList)[i]
  
  message(paste0("working on ", xName))
  
  seuList[[i]] <- subset(seuList[[i]], subset = hto_classification.global == "Singlet")
  
  scrubletOutDir <- "scrubletScores"
  # call doublets using scrublet
  # REMEMBER that this script needs to be run while in the hpap-cite conda env.
  scrubletScoreFn <- paste0(scrubletOutDir, "/", xName, "_doubletScores.txt")
  if (!file.exists(scrubletScoreFn)) {
    message("running scrublet")
    
    # export matrix
    matrixFn <- paste0("matrixForDoubletScoring/", xName, ".mtx")
    Matrix::writeMM(seuList[[i]]@assays$RNA@counts, file = matrixFn)
    
    system(
      glue("python scripts/doubletScorer.py --mm {matrixFn} --sampleName {xName} --outDir {scrubletOutDir}"),
      wait = TRUE
    )
    
    stopifnot(file.exists(scrubletScoreFn))
  }
  
  # read scrublet output and add to seurat object
  message("filtering doublets.")
  
  scores <- read.table(scrubletScoreFn, header = FALSE)[, 1]
  stopifnot(length(scores) == length(seuList[[i]]$orig.ident))
  
  seuList[[i]] <- AddMetaData(seuList[[i]], scores, col.name = "scrubletScore")
  
  # filter out scrublet called doublets
  seuList[[i]] <- subset(seuList[[i]], subset = scrubletScore < 0.25)

  # add metadata
  runN <- as.character(seuList[[i]]$runN[1])
  
  message(paste0("loading metadata: ", runN))
  
  seuDf <- data.frame(
    HTO_DNA_ID = as.character(seuList[[i]]$hash.ID)
  ) %>%
    left_join(sampleMetaDict[[runN]]) %>%
    dplyr::select(-HTO_DNA_ID, -SampleID, -Run_VW, -Run, -rna_fastq_dir, -adthto_fastq_dir, -notes)
  
  rownames(seuDf) <- colnames(seuList[[i]])
  seuList[[i]] <- AddMetaData(seuList[[i]], seuDf)
}

heatshockGenes <- read.csv("miscellaneous_gene_lists/hsp_genes.tsv", sep = "\t")
seuListWithScore <- lapply(seuList, function(x) {
  y <- AddModuleScore(
    object = x,
    features = list(heatshockGenes$Approved.symbol),
    name = "heatShockProgram",
    seed = 42)
  
  return(y)
})


seuDfCompiled <- lapply(seuListWithScore, function(x) {
  df <- FetchData(
    object = x,
    vars = c("Tissue", "heatShockProgram1", "Disease_Status", "DonorID")) %>%
    mutate(Disease_Status = factor(Disease_Status, levels = c("ND", "AAb+", "T1D")))
  
  return(df)
})

seuDf <- bind_rows(seuDfCompiled)

seuDf <- seuDf %>%
  mutate(tissueCondensed = case_when(
    grepl("pLN", Tissue) ~ "pLN",
    grepl("(MES|SMA)", Tissue) ~ "mesLN",
    TRUE ~ Tissue))

p1 <- seuDf %>%
  ggplot(aes(x = DonorID, y = heatShockProgram1, fill = tissueCondensed)) +
  geom_boxplot(position = position_dodge2(width = 0.8, preserve = "single"),
    outlier.color = "#00000050",
    outlier.size = 0.7) +
  theme_bw() +
  facet_wrap(~ Disease_Status, scales = "free_x") +
  subplotTheme +
  theme(
    legend.position = "bottom",
    legend.margin = margin(t = -5, b = 0),
    legend.text = element_text(size = BASEPTFONTSIZE),
    legend.title = element_text(size = BASEPTFONTSIZE),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.direction = "horizontal") +
  geom_hline(yintercept = 0.2, linetype = "dotted") +
  labs(
    x = "Donor",
    y = "HSP Score",
    fill = "Tissue")



################################################################################
# S2D - cold ischemia
# 
# code taken from @Greg
################################################################################
dfLineageFilter <- readRDS("figures/greg_flow_data/rds/dfLineageFilter.rds")

MetaFactors <- c("HPAP Donor", "Disease Status", "Tissue", "cold_ischemia", "Sex", "Age", "Staining_Flow_Researcher")
dfMetaCheck <- dfLineageFilter %>%
  dplyr::select(all_of(MetaFactors)) %>%
  filter(Tissue != "Spleen") %>%
  distinct()


coldIsDisease <- dfMetaCheck %>%
  mutate(`Disease Status` = fct_relevel(`Disease Status`, "ND", "AAb+", "T1D")) %>%
  ggplot(aes(`Disease Status`, cold_ischemia)) +
  geom_boxplot() +
  ylab("Cold Ischemia Time (hr)") +
  theme_classic() +
  subplotTheme +
  theme(
    axis.title.x = element_blank())


################################################################################
# S2E - cold ischemia ranked
# 
# code taken from @Greg
################################################################################
#generate a stats-friendly dataframe for all immune populations
dfAllPopsFreqStats <- dfLineageFilter %>%
  pivot_wider(names_from = metric, values_from = value)
colnames(dfAllPopsFreqStats) <- gsub("\\+ \\(\\%T ?cells\\)", "", colnames(dfAllPopsFreqStats))
colnames(dfAllPopsFreqStats) <- gsub(" ", "_", colnames(dfAllPopsFreqStats))

#inputs and function that generates correlation p values for cold ischemia v all immune populations
coldIsFreqStats <- filter(dfAllPopsFreqStats, !is.na(cold_ischemia)) %>%
  filter(Disease_Status == "ND") %>%
  filter(Tissue != "Spleen")
coldIsPops <- colnames(coldIsFreqStats[15:ncol(coldIsFreqStats)])

dfcoldIsCorrStats <- lapply(coldIsPops, function(x) {
  corrx <- cor_test(coldIsFreqStats, x, cold_ischemia, method = "spearman")
}) %>%
  bind_rows()

#for visualizing plots of interest (correlation value >= 0.3 or <= -0.3)
dfcoldIsPops2 <- dfcoldIsCorrStats %>%
  filter(cor >= 0.3 | cor <= -0.3) %>%
  arrange(desc(cor))
coldIsPops2 <- dfcoldIsPops2$var1

dfcoldIsPopsLebel <- dfcoldIsPops2 %>%
  filter(!str_detect(var1, "CD._T.*_CD127.CD27.")) %>%
  filter(!str_detect(var1, "CD._Tn_CD27.")) %>%
  filter(!str_detect(var1, "CD._Tn_CD127.")) %>%
  filter(!str_detect(var1, "CD._T.*_CD127\\+CD27\\+")) %>%
  filter(!str_detect(var1, "^CD._Mem_CD.*$")) %>%
  filter(!str_detect(var1, "^CD._Mem_.*\\+$")) %>%
  filter(cor >= 0.3 | cor <= -0.3)


#ranked correlation plot of cold ischemia effects on immune population frequency
coldIsRankedPlot <- dfcoldIsCorrStats %>%
  filter(!is.na(cor)) %>%
  filter(!str_detect(var1, "CD._T.*_CD127.CD27.")) %>%
  filter(!str_detect(var1, "CD._Tn_CD27.")) %>%
  filter(!str_detect(var1, "CD._Tn_CD127.")) %>%
  filter(!str_detect(var1, "CD._T.*_CD127\\+CD27\\+")) %>%
  filter(!str_detect(var1, "^CD._Mem_CD.*$")) %>%
  filter(!str_detect(var1, "^CD._Mem_.*\\+$")) %>%
  ggplot(aes(reorder(var1, -cor), cor)) +
  geom_hline(yintercept = c(0.3, -0.3),
    linetype = "dashed",
    color = "grey") +
  geom_point(size = 0.3, alpha = 0.75) + 
  geom_point(data = dfcoldIsPopsLebel, aes(var1, cor), color = "blue", size = 0.3) +
  geom_text_repel(data = filter(dfcoldIsPopsLebel), 
    aes(var1, cor, label = var1), size = 2, seed = 42, color = "#000000", 
    segment.size = 0.1, segment.alpha = 0.8, max.overlaps = 50) +
  ylab("Spearmann Corr. Value") +
  xlab("Immune Population") +
  theme_classic() +
  subplotTheme +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.length.x = unit(0, "in"),
  )

################################################################################
# Final layout and plot all
################################################################################
layout <- c(
  patchwork::area(1, 1, 3, 10), # a
  patchwork::area(4, 1, 6, 10), # b
  patchwork::area(7, 1, 8, 5), # c
  patchwork::area(7, 6, 8, 10), # d
  patchwork::area(9, 1, 11, 7) # b
)

p <- wrap_elements(plot = figA) +
  wrap_elements(plot = p1) +
  wrap_elements(plot = p2) +
  wrap_elements(plot = coldIsDisease) +
  wrap_elements(plot = coldIsRankedPlot) +
  plot_annotation(tag_levels = list(LETTERS[1:5])) +
  plot_layout(design = layout) &
  plotTagTheme

saveFinalFigure(
  plot = p,
  prefixDir = "figures/outs",
  fn = "sfig2_v2_final",
  devices = c("pdf", "png"),
  addTimestamp = TRUE,
  gwidth = 8,
  gheight = 11)

