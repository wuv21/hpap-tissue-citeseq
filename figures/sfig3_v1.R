# note that this code is written to be run from the project base directory
# renv::load("/data/hpap-citeseq/hpap-citeseq-analysis")

source("figures/genericFigureSettings.R")
source("scripts/dimPlots.R")
source("scripts/deg.R")
library(Seurat)
library(cowplot)
library(presto)

set.seed(42)

################################################################################
# S3A - hsp score
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
  
  seuDf <- data.frame(
    HTO_DNA_ID = as.character(seuList[[i]]$hash.ID)
  ) %>%
    left_join(sampleMetaDict[[runN]]) %>%
    select(-HTO_DNA_ID, -SampleID, -Run_VW, -Run, -rna_fastq_dir, -adthto_fastq_dir, -notes)
  
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
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.direction = "horizontal") +
  geom_hline(yintercept = 0.2, linetype = "dotted") +
  labs(
    x = "Donor",
    y = "HSP Score",
    fill = "Tissue")


################################################################################
# S3B - hsp filter cutoff
################################################################################
load("rds/preManualAnnotSubset.rData")

seuMerged$Disease_Status <- factor(seuMerged$Disease_Status, levels = c("ND", "AAb+", "T1D"))

seuFiltDf <- FetchData(
  object = seuMerged,
  vars = c("TissueCondensed", "heatShockProgram1", "Disease_Status", "DonorID")) %>%
  mutate(Disease_Status = factor(Disease_Status, levels = c("ND", "AAb+", "T1D")))

p2 <- seuFiltDf %>%
  ggplot(aes(x = heatShockProgram1, color = TissueCondensed)) +
  geom_density() +
  theme_classic() +
  geom_vline(xintercept = quantile(seuMerged$heatShockProgram1, probs = 0.95), color = "black", linetype = "dotted") +
  labs(x = "HSP score", y = "Density", color = "Tissue") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))

# ################################################################################
# # Final layout and plot all
# ################################################################################
layout <- c(
  patchwork::area(1, 1, 3, 6), # a
  patchwork::area(4, 1, 6, 6) # b
)

p <- wrap_elements(plot = p1) +
  wrap_elements(plot = p2) +
  plot_annotation(tag_levels = list(LETTERS[1:2])) +
  plot_layout(design = layout) &
  plotTagTheme

saveFinalFigure(
  plot = p,
  prefixDir = "figures/outs",
  fn = "sfig3_final",
  devices = c("pdf", "png"),
  addTimestamp = TRUE,
  gwidth = 7.5,
  gheight = 7)

