# load libraries
silentLoadLibrary <- function(x) {
  suppressWarnings(suppressMessages(library(x, character.only = TRUE)))
}

libraries <- c(
  "dplyr",
  "tidyr",
  "stringr",
  "ggplot2",
  "Seurat",
  "harmony",
  "glue",
  "patchwork",
  "DropletUtils",
  "clustree",
  "future"
)

invisible(lapply(libraries, silentLoadLibrary))

# reticulate::use_condaenv("hpap-cite", required = TRUE)

source("scripts/generic.R")
source("scripts/adtProcessing.R")
source("scripts/initalSeuratListMaker.R")
source("scripts/clustering.R")
source("scripts/dimPlots.R")

# create list of seurat objects
rnaOut <- "../hpap-citeseq-data/rna"
adtOut <- "../hpap-citeseq-data/adt"
htoOut <- "../hpap-citeseq-data/hto"

# load tsa catalog
tsaCatalog <- readRDS("rds/tsa_catalog.rds")

# unhashed sample run info
unhashedSampleMeta <- read.csv("metadata/unhashed_sample_metadata.csv")

# load merged dataset if it exists
seuMergedFn <- "rds/seuMerged.rds"
harmonizedFn <- "rds/seuMergedAndHarmonized.rds"
umapRdsFn <- "rds/seuMergedAndHarmonized_withRNAClusters3_UMAP.rds"

if (file.exists(seuMergedFn) & file.exists(harmonizedFn)) {
  message("harmonized and merged file already exists, so skipping initial merging")
  
} else if (file.exists(seuMergedFn)) {
  message("loading previously merged rds file")
  seuMerged <- readRDS(seuMergedFn)
  
} else {
  # prepare initial seurat list with all three modalities
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
    
    # export matrix
    matrixFn <- paste0("matrixForDoubletScoring/", xName, ".mtx")
    Matrix::writeMM(seuList[[i]]@assays$RNA@counts, file = matrixFn)
    
    # call doublets using scrublet
    message("running scrublet")
    scrubletOutDir <- "scrubletScores"
    
    # REMEMBER that this script needs to be run while in the hpap-cite conda env.
    system(
      glue("python scripts/doubletScorer.py --mm {matrixFn} --sampleName {xName} --outDir {scrubletOutDir}"),
      wait = TRUE
    )
    
    # read scrublet output and add to seurat object
    message("filtering doublets.")

    scrubletScoreFn <- paste0(scrubletOutDir, "/", xName, "_doubletScores.txt")
    stopifnot(file.exists(scrubletScoreFn))
    
    scores <- read.table(scrubletScoreFn, header = FALSE)[, 1]
    stopifnot(length(scores) == length(seuList[[i]]$orig.ident))
    
    seuList[[i]] <- AddMetaData(seuList[[i]], scores, col.name = "scrubletScore")
  
    # filter out scrublet called doublets
    seuList[[i]] <- subset(seuList[[i]], subset = scrubletScore < 0.25)
    
    message("finding variable features.")
    seuList[[i]] <- FindVariableFeatures(seuList[[i]])
    
    # TODO consider running this without the control abs
    VariableFeatures(seuList[[i]], assay = "adt") <- rownames(seuList[[i]][["adt"]])
    seuList[[i]] <- NormalizeData(seuList[[i]], assay = "adt", normalization.method = 'CLR', margin = 2)
    
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
  
  message('merging datasets')
  seuMerged <- merge(seuList[[1]], seuList[2:length(seuList)])
  VariableFeatures(seuMerged, assay = "RNA") <- SelectIntegrationFeatures(
    object.list = seuList,
    assay = rep("RNA", length(seuList)))
  
  seuMerged <- ScaleData(seuMerged, assay = "RNA") %>%
    RunPCA(assay = "RNA")
  
  saveRDS(seuMerged, seuMergedFn)
}


# harmonize merged seurat
if (file.exists(harmonizedFn) & file.exists(umapRdsFn)) {
  message("next rds (clustered + umap) already found so skipping import of harmonized rds")
  
} else if (file.exists(harmonizedFn)) {
  message("loading previously merged and harmonized rds")
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
  
  message('harmonizing RNA')
  seuMerged <- RunHarmony(
    object = seuMerged,
    reduction.save = "harmonyRNA",
    reduction = "pca",
    assay.use = "RNA",
    project.dim = FALSE,
    dims.use = 1:50,
    verbose = TRUE,
    group.by.vars = harmonyGroupingVars)
  
  message("working on adt")
  VariableFeatures(seuMerged, assay = "adt") <- rownames(seuMerged[["adt"]])
  
  message("scaling adt")
  seuMerged <- ScaleData(seuMerged, assay = "adt")
  
  message("calculating PCA for adt")
  seuMerged <- RunPCA(seuMerged, assay = "adt", reduction.name = "apca")
  
  message("harmonizing ADT")
  seuMerged <- RunHarmony(
    object = seuMerged,
    assay.use = "adt",
    reduction = "apca",
    group.by.vars = harmonyGroupingVars,
    reduction.save = "harmonyADT")
  
  
  message("performing WNN")
  
  DefaultAssay(seuMerged) <- "RNA"
  seuMerged <- FindMultiModalNeighbors(
    seuMerged,
    # prune.SNN = 0,
    reduction.list = list("harmonyRNA", "harmonyADT"),
    dims.list = list(1:30, 1:30), modality.weight.name = c("RNA.weight", "ADT.weight"))
  
  saveRDS(seuMerged, harmonizedFn)
}

# find clusters based on rna only
if (!file.exists(umapRdsFn)) {
  seuMerged <- FindNeighbors(seuMerged, assay = "RNA", reduction = "harmonyRNA", dims = 1:30)
  seuMerged <- FindClusters(
    seuMerged,
    algorithm = 4,
    method = "igraph",
    resolution = c(0.5, 0.75, 1),
    verbose = TRUE)
  
  
  seuMerged <- RunUMAP(seuMerged,
    dims = 1:30,
    reduction.name = "rna.umap",
    reduction.key = "rnaUMAP_",
    reduction = "harmonyRNA")
  
  saveRDS(seuMerged, umapRdsFn)
  
} else {
  seuMerged <- readRDS(umapRdsFn)
}


smallMultipleUmaps(seu = seuMerged, parameter = "RNA_snn_res.1")

smallMultipleUmaps(seu = seuMerged, parameter = "DonorID", ncol = 6, height = 5)
smallMultipleUmaps(seu = seuMerged, parameter = "Tissue", ncol = 5, height = 2.5)
smallMultipleUmaps(seu = seuMerged, parameter = "runN", ncol = 7, height = 1.5)
smallMultipleUmaps(seu = seuMerged, parameter = "Disease_Status", ncol = 3, width = 4, height = 1.5)
smallMultipleUmaps(seu = subset(seuMerged, subset = Tissue != "Spleen"),
  parameter = "Disease_Status", filename = "Disease_Status_noSpleen", ncol = 3, width = 4, height = 1.5)

p <- DimPlotCustom(seuMerged,
  groupBy = c("DonorID", "Tissue", "runN"),
  groupByTitles = c("Donor", "Tissue", "Run"),
  nLegendCols = 3,
  nCols = 3) &
  textSizeOnlyTheme
savePlot(plot = p, fn = "umap_all", devices = "png", gwidth = 9, gheight = 4)


clusterOrder <- unique(seuMerged$RNA_snn_res.1)
clusterOrder <- factor(clusterOrder, levels = str_sort(clusterOrder, numeric = TRUE))

labellerAdt <- tsaCatalog$cleanName
names(labellerAdt) <- paste0("adt_", tsaCatalog$DNA_ID)


# load in manual annotations
manualAnnot <- read.csv("manualClusterAnnotations/RNAsnn1.csv", header = FALSE)
manualAnnotDict <- manualAnnot[, 2]
names(manualAnnotDict) <- as.character(manualAnnot[, 1])

manualAnnot <- manualAnnotDict[as.character(seuMerged$RNA_snn_res.1)]
names(manualAnnot) <- names(seuMerged$RNA_snn_res.1)

seuMerged <- AddMetaData(
  object = seuMerged,
  metadata = manualAnnot,
  col.name = 'manualAnnot'
)
seuMerged <- subset(seuMerged, subset = manualAnnot != "")


tissueCondensed <- case_when(
  grepl("pLN", seuMerged$Tissue) ~ "pLN",
  grepl("(MES|SMA)", seuMerged$Tissue) ~ "mesLN",
  TRUE ~ seuMerged$Tissue
)
names(tissueCondensed) <- names(seuMerged$Tissue)
seuMerged <- AddMetaData(
  object = seuMerged,
  metadata = tissueCondensed,
  col.name = 'TissueCondensed'
)

rm(manualAnnot)
rm(tissueCondensed)

# overview graphs
smallMultipleUmaps(seu = seuMerged, parameter = "manualAnnot")

seuDf <- FetchData(object = seuMerged, vars = c("TissueCondensed", "Disease_Status", "manualAnnot", "DonorID")) %>%
  mutate(Disease_Status = factor(Disease_Status, levels = c("ND", "AAb+", "T1D")))

# cluster proportion graph
seuDf %>%
  group_by(Disease_Status, TissueCondensed, DonorID, manualAnnot) %>%
  summarize(Count = n()) %>%
  group_by(DonorID, TissueCondensed) %>%
  mutate(prop = Count / sum(Count)) %>%
  ggplot(aes(x = manualAnnot, y = prop, fill = Disease_Status)) +
  geom_point(alpha = 0.9, color = "#555555", pch = 21, aes(group = Disease_Status), position = position_dodge(width = 0.8)) +
  facet_wrap(~ TissueCondensed, ncol = 1) +
  scale_fill_manual(values = c("#4daf4a", "#377eb8", "#e41a1c")) +
  labs(y = "Proportion within tissue of speciifc donor",
    fill = "Status") +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.margin = margin(t = -20, b = 0),
    plot.margin = margin(10, 50, 10, 50, unit = "pt"),
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_blank(),
    panel.grid.major.y = element_line(color = "#eeeeee"),
    panel.grid.major.x = element_line(color = "#eeeeee", size = 0.8))


# number of cells for each cluster
seuDf %>%
  group_by(Disease_Status, TissueCondensed, manualAnnot) %>%
  summarize(Count = n()) %>%
  ggplot(aes(x = Count, y = manualAnnot, fill = Disease_Status)) +
  geom_point(alpha = 0.9, color = "#555555", pch = 21, aes(group = Disease_Status), position = position_dodge(width = 0.8)) +
  facet_wrap(~ TissueCondensed, nrow = 1) +
  scale_x_continuous(labels = scales::label_number(suffix = "K", scale = 1e-3, big.mark = ",")) +
  theme_classic() +
  labs(x = "Number of cells", fill = "Disease Status") +
  theme(
    legend.position = "bottom",
    plot.margin = margin(10, 30, 10, 30, unit = "pt"),
    strip.background = element_blank(),
    panel.spacing = unit(15, units = "pt"),
    axis.title.y = element_blank(),
    text = element_text(size = BASEPTFONTSIZE + 2),
    panel.grid.major.y = element_line(color = "#dddddd80", size = 6),
    panel.grid.major.x = element_line(color = "#55555580", size = 0.8, linetype = "dotted"))

# TODO make a manual annot sorter
