# load libraries
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(Seurat)
library(harmony)
library(glue)
library(patchwork)
library(DropletUtils)
library(clustree)

source("scripts/adtProcessing.R")
source("scripts/generic.R")
source("scripts/initalSeuratListMaker.R")
source("scripts/dimPlots.R")

# create list of seurat objects
rnaOut <- "../hpap-citeseq-data/rna"
adtOut <- "../hpap-citeseq-data/adt"
htoOut <- "../hpap-citeseq-data/hto"

# load tsa catalog
tsaCatalog <- readRDS("rds/tsa_catalog.rds")

# unhashed sample run info
unhashedSampleMeta <- read.csv("metadata/unhashed_sample_metadata.csv")

# load mered dataset if it exists
seuMergedFn <- "rds/seuMerged.rds"
if (file.exists(seuMergedFn)) {
  seuMerged <- readRDS(seuMergedFn)
} else {
  # prepare initial seurat list with all three modalities
  seuList <- rdsCheckAndRun(
    fn = "rds/allSeuList_v1.rds",
    f = initialSeuratListMaker,
    version = "Merge",
    sampleIds = unhashedSampleMeta$id,
    sampleRuns = unhashedSampleMeta$run,
    rnaOutDir = rnaOut,
    adtOutDir = adtOut,
    htoOutDir = htoOut,
    plotDir = "outs/qc")
  
  names(seuList) <- unhashedSampleMeta$id
  seuList[["20220531_hpap_2_Pool_Run1_Well2"]] <- NULL
  
  sampleMeta <- read.csv("metadata/sample_metadata.csv")
  sampleMetaDict <- sampleMeta %>% 
    group_split(Run)
  
  names(sampleMetaDict) <- sapply(sampleMetaDict, function(x) {return(as.character(x[1, "Run"]))})
  
  for (i in seq(1, length(seuList))) {
    xName <- names(seuList)[i]
    
    seuList[[i]] <- subset(seuList[[i]], subset = hto_classification.global == "Singlet")
    
    seuList[[i]] <- FindVariableFeatures(seuList[[i]])
    
    VariableFeatures(seuList[[i]], assay = "adt") <- rownames(seuList[[i]][["adt"]])
    seuList[[i]] <- NormalizeData(seuList[[i]], assay = "adt", normalization.method = 'CLR', margin = 2)
    
    # add metadata
    runN <- as.character(seuList[[i]]$runN[1])
    
    seuDf <- data.frame(
      HTO_DNA_ID = as.character(seuList[[i]]$hash.ID)
    ) %>%
      left_join(sampleMetaDict[[runN]]) %>%
      select(-HTO_DNA_ID, -SampleID)
    
    rownames(seuDf) <- colnames(seuList[[i]])
    seuList[[i]] <- AddMetaData(seuList[[i]], seuDf)
  }
  
  print('merging datasets')
  seuMerged <- merge(seuList[[1]], seuList[2:length(seuList)])
  VariableFeatures(seuMerged, assay = "RNA") <- SelectIntegrationFeatures(
    object.list = seuList,
    assay = rep("RNA", length(seuList)))
  
  seuMerged <- ScaleData(seuMerged, assay = "RNA") %>%
    RunPCA(assay = "RNA")
  
  saveRDS(seuMerged, seuMergedFn)
}


# harmonize merged seurat
harmonizedFn <- "rds/seuMergedAndHarmonized.rds"
if (file.exists(harmonizedFn)) {
  seuMerged <- readRDS(harmonizedFn)
} else {
  # ensure that all metadata is character
  # for some reason, harmony doesn't like integer metadata
  harmonyGroupingVars <- c("runN", "DonorID")
  
  for (v in harmonyGroupingVars) {
    tmpHolder <- as.character(seuMerged[[v]][, 1])
    seuMerged <- AddMetaData(seuMerged, tmpHolder, col.name = v)
    
    stopifnot(sum(tmpHolder == seuMerged[[v]][, 1]) == length(tmpHolder))
  }
  
  print('harmonizing RNA')
  seuMerged <- RunHarmony(
    object = seuMerged,
    reduction.save = "harmonyRNA",
    reduction = "pca",
    assay.use = "RNA",
    project.dim = FALSE,
    dims.use = 1:50,
    verbose = TRUE,
    group.by.vars = harmonyGroupingVars)
  
  print("working on adt")
  VariableFeatures(seuMerged, assay = "adt") <- rownames(seuMerged[["adt"]])
  
  print("scaling adt")
  seuMerged <- ScaleData(seuMerged, assay = "adt")
  
  print("calculating PCA for adt")
  seuMerged <- RunPCA(seuMerged, assay = "adt", reduction.name = "apca")
  
  print('harmonizing ADT')
  seuMerged <- RunHarmony(
    object = seuMerged,
    assay.use = "adt",
    reduction = "apca",
    group.by.vars = harmonyGroupingVars,
    reduction.save = "harmonyADT")
  
  
  print('performing WNN')
  
  DefaultAssay(seuMerged) <- "RNA"
  seuMerged <- FindMultiModalNeighbors(
    seuMerged,
    reduction.list = list("harmonyRNA", "harmonyADT"),
    dims.list = list(1:30, 1:30), modality.weight.name = c("RNA.weight", "ADT.weight"))
  
  saveRDS(seuMerged, harmonizedFn)
}

seuMerged <- RunUMAP(seuMerged, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

DimPlotCustom(seuMerged,
  groupBy = c("DonorID", "Tissue", "runN"),
  groupByTitles = c("Donor", "Tissue", "Run"),
  nCols = 3)

seuMerged <- FindClusters(seuMerged, graph.name = "wsnn", algorithm = 3, resolution = c(0.5, 0.75, 1), verbose = TRUE)
saveRDS(seuMerged, "rds/seuMergedAndHarmonized_withClusters.rds")


clustreeRes <- clustree(seuMerged, prefix = "wsnn_res.", node_colour = "sc3_stability")

p1 <- FeaturePlotCustom(
  seuMerged,
  markers = c("A0034","A0046","A0083", "A0100", "A0053", "A0081"),
  modality = "adt",
  tsa_catalog = tsaCatalog,
  max.cutoff = 5
)

p2 <- FeaturePlotCustom(
  seuMerged,
  markers = c("CD3E", "CD8A", "NKG7", "MS4A1", "CST3", "MS4A7"), 
  modality = "RNA",
  tsa_catalog = tsaCatalog
)

p1 / p2

# table of cell counts...
seuDf <- data.frame(
  runN = seuMerged$runN,
  tissue = seuMerged$Tissue,
  disease = seuMerged$Disease_Status,
  donorID = seuMerged$DonorID,
  wsnn_cluster0.5 = seuMerged$wsnn_res.0.5
) %>%
  mutate(disease = factor(disease, levels = c("ND", "AAb+", "T1D")))


seuDf %>%
  group_by(disease, tissue, donorID, wsnn_cluster0.5) %>%
  summarize(Count = n()) %>%
  group_by(donorID, tissue) %>%
  mutate(prop = Count / sum(Count)) %>%
  mutate(wsnn_cluster0.5 = as.numeric(wsnn_cluster0.5)) %>%
  ggplot(aes(x = wsnn_cluster0.5, y = prop, color = disease)) +
  geom_point() +
  facet_wrap(~ tissue, ncol = 1)


seuDf %>%
  group_by(disease, tissue, donorID) %>%
  summarize(Count = n()) %>%
  ggplot(aes(y = Count, x = tissue, color = disease)) +
  geom_point(size = 3, alpha = 0.7, position = position_dodge(width = 0.5)) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_line(color = "#eeeeee50", size = 30),
    panel.grid.minor.y = element_blank()
  ) +
  scale_color_manual(values = c("#0000ff", "#880000", "#ff0000")) +
  labs(y = "# of cells",
    x = "Tissue",
    color = "Disease status")

data.frame(
  runN = seuMerged$runN,
  tissue = seuMerged$Tissue,
  disease = seuMerged$Disease_Status,
  donorID = seuMerged$DonorID,
  nCountHTO = seuMerged$nCount_hto,
  nCountRNA = seuMerged$nCount_RNA,
  nCountADT = seuMerged$nCount_adt,
  nFeatureRNA = seuMerged$nFeature_RNA
) %>%
  group_by(runN, tissue, donorID) %>%
  dplyr::summarize(
    avgHTOCount = mean(nCountHTO),
    avgRNACount = mean(nCountRNA),
    avgADTCount = mean(nCountADT),
    avgRNAFeatures = mean(nFeatureRNA)
  ) %>%
  pivot_longer(cols = all_of(starts_with("avg")), names_to = "modality", values_to = "avg") %>%
  ggplot(aes(x = modality, y = avg, color = donorID, shape = tissue)) +
  geom_jitter() +
  labs(y = "Avg value (over all qc passing cells)") +
  facet_wrap(~ runN, nrow = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# use propeller/speckle for downstream analysis.
  