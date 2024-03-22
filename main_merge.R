################################################################################
# load libraries and custom scripts
################################################################################

# uncomment if running in RProject
renv::load("/data/hpap-citeseq/hpap-citeseq-analysis")

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
  "future",
  "presto",
  "msigdbr",
  "fgsea",
  "ComplexHeatmap",
  "WGCNA",
  "hdWGCNA"
)

# reticulate::use_condaenv("hpap-cite")

invisible(lapply(libraries, silentLoadLibrary))

# general scripts
source("scripts/generic.R")
source("scripts/adtProcessing.R")
source("scripts/initalSeuratListMaker.R")
source("scripts/clustering.R")

# plotting scripts
source("scripts/dimPlots.R")
source("scripts/deg.R")

set.seed(42)

################################################################################
# preliminary setup and loading of meta/data files
################################################################################
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
    scrubletScoreFn <- paste0(scrubletOutDir, "/", xName, "_doubletScores.txt")
    if (!file.exists(scrubletScoreFn)) {
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
  
  message('merging datasets')
  # filter out specific donors due to poor viability during thaw and initial qc
  seuList <- lapply(seuList, function(x) {
    x <- x[, x@meta.data$DonorID != "HPAP-043" & x@meta.data$DonorID != "HPAP-024" & !(x@meta.data$DonorID == "HPAP-098" & x@meta.data$Tissue == "Spleen")]
    
    message("finding variable features.")
    seuList[[i]] <- FindVariableFeatures(seuList[[i]])
    
    VariableFeatures(seuList[[i]], assay = "adt") <- rownames(seuList[[i]][["adt"]])
    seuList[[i]] <- NormalizeData(seuList[[i]], assay = "adt", normalization.method = 'CLR', margin = 2)
    return(x)
  })
  
  seuMerged <- merge(seuList[[1]], seuList[2:length(seuList)])
  VariableFeatures(seuMerged, assay = "RNA") <- SelectIntegrationFeatures(
    object.list = seuList,
    assay = rep("RNA", length(seuList)))
  
  seuMerged <- ScaleData(seuMerged, assay = "RNA") %>%
    RunPCA(assay = "RNA")
  
  saveRDS(seuMerged, seuMergedFn)
  
  rm(seuList)
  gc()
}

################################################################################
# harmony batch effect correction against run number and donor ID.
################################################################################
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
  # normalize since this info isn't saved after merging for some reason...
  seuMerged <- NormalizeData(seuMerged, assay = "adt", normalization.method = 'CLR', margin = 2)
  
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

################################################################################
# cluster finding by RNA modality only
################################################################################
# plan("multisession", workers = 2)
# options(future.globals.maxSize = 175000 * 1024^2)
if (!file.exists(umapRdsFn)) {
  seuMerged <- FindNeighbors(seuMerged, assay = "RNA", reduction = "harmonyRNA", dims = 1:30)
  seuMerged <- FindClusters(
    seuMerged,
    algorithm = 4,
    method = "igraph",
    resolution = c(0.75, 1),
    verbose = TRUE)
  
  saveRDS(seuMerged, "rds/seuMerged_withClusters")

seuMerged <- RunUMAP(seuMerged,
  dims = 1:30,
  reduction.name = "rna.umap",
  reduction.key = "rnaUMAP_",
  n.neighbors = 50,
  min.dist = 0.2,
  reduction = "harmonyRNA")
  
  saveRDS(seuMerged, umapRdsFn)
  
} else {
  seuMerged <- readRDS(umapRdsFn)
}

clusterOrder <- unique(seuMerged$RNA_snn_res.1)
clusterOrder <- factor(clusterOrder, levels = str_sort(clusterOrder, numeric = TRUE))

labellerAdt <- tsaCatalog$cleanName
names(labellerAdt) <- paste0("adt_", tsaCatalog$DNA_ID)

smallMultipleUmaps(seu = seuMerged, parameter = "RNA_snn_res.1")
smallMultipleUmaps(seu = seuMerged, parameter = "DonorID", ncol = 6, height = 5)
smallMultipleUmaps(seu = seuMerged, parameter = "Tissue", ncol = 4, height = 5)
smallMultipleUmaps(seu = seuMerged, parameter = "runN", ncol = 4, height = 5)
smallMultipleUmaps(seu = seuMerged, parameter = "Disease_Status", ncol = 3, width = 4, height = 1.5)

p <- DimPlotCustom(seuMerged,
  groupBy = c("DonorID", "Tissue", "runN"),
  groupByTitles = c("Donor", "Tissue", "Run"),
  nLegendCols = 3,
  nCols = 3) &
  textSizeOnlyTheme
savePlot(plot = p, fn = "umap_all", devices = "png", gwidth = 9, gheight = 4)

# the giant umap
DimPlot(seuMerged, reduction = "rna.umap", group.by = "RNA_snn_res.1") +
  theme(
    plot.title = element_blank(),
    legend.position = "bottom",
    legend.justification = "center"
  ) +
  scale_color_discrete(limits = levels(clusterOrder))

# add in condensed tissue state
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
rm(tissueCondensed)

# add heatshock protein score
heatshockGenes <- read.csv("hsp_genes.tsv", sep = "\t")
seuMerged <- AddModuleScore(
  object = seuMerged,
  features = list(heatshockGenes$Approved.symbol),
  name = "heatShockProgram",
  seed = 42
)

seuDf <- FetchData(object = seuMerged,
  vars = c("TissueCondensed", "heatShockProgram1", "Disease_Status", "RNA_snn_res.1", "DonorID")) %>%
  mutate(Disease_Status = factor(Disease_Status, levels = c("ND", "AAb+", "T1D"))) %>%
  mutate(cell_id = Cells(seuMerged))

seuDf %>%
  ggplot(aes(x = DonorID, y = heatShockProgram1, fill = TissueCondensed)) +
  geom_boxplot(position = position_dodge2(width = 0.8, preserve = "single")) +
  theme_bw() +
  geom_hline(yintercept = 0.2, linetype = "dotted") +
  facet_wrap(~ Disease_Status, scales = "free_x")

rm(seuDf)
save.image("rds/preManualAnnotSubset.rData")

baseMarkers <- c(
  "A0034", #CD3
  "A0072", #CD4
  "A0046", #CD8
  "A0081", #CD14
  "A0083", #CD16
  "A0047", #CD56
  "A0050", #CD19
  "A0087", #CD45RO
  "A0063", #CD45RA
  "A0154", #CD27
  "A0386", #CD28
  "A0156", #CD95
  "A0390", #CD127
  "A0146", #CD69
  "A0085", #CD25
  "A0089", #TIGIT
  "A0088", #PD1
  "A0141", #CCR5
  "A0144", #CXCR5
  "A0149"  #CD161
)

# base adt dot style
p <- DotPlot(seuMerged, features = paste0("adt_", baseMarkers), group.by = "RNA_snn_res.1") +
  scale_y_discrete(limits = levels(clusterOrder)) +
  textSizeOnlyTheme +
  theme(
    panel.grid.major.y = element_line(color = "#eeeeee"),
    axis.text.x = element_text(angle = 45, size = BASEPTFONTSIZE, hjust = 1)) +
  scale_x_discrete(labels = labellerAdt)


basePanelRNA <- c(
  "IL7R", # cd127
  "CCR7",
  "CD14",
  "CSF3R", #cd114
  "LYZ",
  "S100A4",
  "MS4A1",
  "CD8A",
  "FCGR3A",
  "MS4A7",
  "GNLY",
  "NKG7",
  "FCER1A",
  "CST3",
  "PPBP", #cxcl7 (platelet basic protein)
  "CXCR5",
  "TOX",
  "TCF7",
  "TBX21",
  "CD69",
  "CD3D",
  "CD3G",
  "CD40LG",
  "ICOS",
  "FOXP3",
  "GATA3",
  "IKZF2", # helios
  "RORC",
  "BCL6",
  "PRDM1", #blimp
  "CD38",
  "TNFRSF17" #BCMA
)

# base rna dot style
p <- DotPlot(seuMerged, features = basePanelRNA, group.by = "RNA_snn_res.1") +
  scale_y_discrete(limits = levels(clusterOrder)) +
  textSizeOnlyTheme +
  theme(
    panel.grid.major.y = element_line(color = "#eeeeee"),
    axis.text.x = element_text(angle = 45, size = BASEPTFONTSIZE, hjust = 1))

# # cluster DE for adt and rna
# # get differential expression of various ADT markers
# getClusterDifferences(
#   seu = seuMerged,
#   clusterName = "RNA_snn_res.1",
#   tsvFn = "outs/tsv/DE_ADT_cluster.tsv",
#   assay = "adt",
#   tsa_catalog = tsaCatalog,
#   findMarkerMethod = "wilcox",
#   only.pos = TRUE,
#   densify = TRUE,
#   max.cells.per.ident = 100000,
#   random.seed = 42)
# 
# options(future.globals.maxSize = 10000 * 1024^2)
# plan("multisession", workers = 8)
# # get differential expression of various RNA markers
# getClusterDifferences(
#   seu = seuMerged,
#   clusterName = "RNA_snn_res.1",
#   assay = "RNA",
#   tsvFn = "outs/tsv/DE_RNA_cluster.tsv",
#   parallelize = TRUE,
#   findMarkerMethod = "wilcox",
#   only.pos = TRUE,
#   densify = TRUE,
#   max.cells.per.ident = 100000,
#   random.seed = 42)
# 


################################################################################
# export for panc db processing
################################################################################
# FetchData(object = seuMerged,
#   vars = c("DonorID", "Tissue", "runN", "well", "hash.ID")) %>%
#   mutate(barcode = stringr::str_match(rownames(.), "[ATGC]+$")[, 1]) %>%
#   dplyr::rename(
#     donorID = DonorID,
#     tissue = Tissue,
#     hashID = hash.ID,
#     projectWellID = well
#   ) %>%
#   write.table(x = ., file = "outs/csv/panc-db_barcodes.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
# 
# sampleMeta %>%
#   select(DonorID, Tissue, Run, HTO_DNA_ID) %>%
#   dplyr::rename(
#     donorID = DonorID,
#     tissue = Tissue,
#     runN = Run,
#     hashID = HTO_DNA_ID,
#   ) %>%
#   write.table(x = ., file = "outs/csv/panc-db_sampleMeta.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
# 
# unhashedSampleMeta %>%
#   mutate(runID = stringr::str_match(id, "(\\d{8}_hpap_\\d)_(.*)")[, 2:3]) %>%
#   mutate(project = runID[, 1]) %>%
#   mutate(well = runID[, 2]) %>%
#   select(-runID) %>%
#   dplyr::rename(
#     projectWellID = id,
#     runN = run
#   ) %>%
#   left_join(sampleMeta %>% select(Run, rna_fastq_dir, adthto_fastq_dir) %>% distinct(), by = c("runN" = "Run")) %>%
#   mutate(rna_fastq_dir = gsub(
#     "/project/bettslab/citeseq_raw/",
#     "/project/betts_shescott_lab/hpap_atlas_citeseq/rna_fastq",
#     rna_fastq_dir)) %>%
#   mutate(adthto_fastq_dir = gsub(
#     "/project/bettslab/citeseq_raw/",
#     "/project/betts_shescott_lab/hpap_atlas_citeseq/adthto_fastq/",
#     adthto_fastq_dir)) %>%
#   write.table(x = ., file = "outs/csv/panc-db_runs.csv", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)


################################################################################
# load in manual annotations and create condensed labeling for cluster and
# phenotypes and look at HSP signatures
################################################################################
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
rm(manualAnnot)

phenotypeCondensed <- case_when(
  grepl("^CD4", seuMerged$manualAnnot) ~ "CD4 T cell",
  grepl("^CD8", seuMerged$manualAnnot) ~ "CD8 T cell",
  grepl("^B", seuMerged$manualAnnot) ~ "B cell",
  grepl("NK", seuMerged$manualAnnot) ~ "NK/ILC cell",
  TRUE ~ "Mono/APC cell"
)
names(phenotypeCondensed) <- names(seuMerged$manualAnnot)
seuMerged <- AddMetaData(
  object = seuMerged,
  metadata = phenotypeCondensed,
  col.name = 'phenotypeCondensed'
)

rm(phenotypeCondensed)
save.image("rds/postManualAnnotSubset.rData")

### START HERE WHEN LOADING PREVIOUS RDATA
load("rds/postManualAnnotSubset.rData")

################################################################################
# post annotation graphs
################################################################################
# HSP signature graphs
# seuDf <- FetchData(object = seuMerged,
#   vars = c("TissueCondensed", "heatShockProgram1", "Disease_Status", "manualAnnot", "RNA_snn_res.1", "DonorID")) %>%
#   mutate(Disease_Status = factor(Disease_Status, levels = c("ND", "AAb+", "T1D")))

# seuDf %>%
#   ggplot(aes(x = DonorID, y = heatShockProgram1, fill = TissueCondensed)) +
#   geom_boxplot(position = position_dodge2(width = 0.8, preserve = "single")) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   geom_hline(yintercept = quantile(seuMerged$heatShockProgram1, probs = 0.90), color = "red") +
#   geom_hline(yintercept = quantile(seuMerged$heatShockProgram1, probs = 0.95), color = "blue")
# 
# 
# seuDf %>%
#   ggplot(aes(x = manualAnnot, y = heatShockProgram1)) +
#   geom_boxplot(position = position_dodge2(width = 0.8, preserve = "single")) +
#   theme_bw() +
#   textSizeOnlyTheme +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   scale_color_discrete(limits = levels(manualClusterOrder)) +
#   theme(plot.margin = margin(3, 3, 3, 3, "lines")) +
#   geom_hline(yintercept = quantile(seuMerged$heatShockProgram1, probs = 0.90), color = "red") +
#   geom_hline(yintercept = quantile(seuMerged$heatShockProgram1, probs = 0.95), color = "blue")
# 
# seuDf %>%
#   ggplot(aes(x = heatShockProgram1, color = TissueCondensed)) +
#   geom_density() +
#   theme_classic() +
#   geom_vline(xintercept = quantile(seuMerged$heatShockProgram1, probs = 0.90), color = "red") +
#   geom_vline(xintercept = quantile(seuMerged$heatShockProgram1, probs = 0.95), color = "blue")

# rm(seuDf)

##### START run this again after loading manual annot rdata
hspCutoff <- quantile(seuMerged$heatShockProgram1, probs = 0.95)

customSortAnnotation <- function(x) {
  priority <- c("B", "CD4", "CD8", "NK", "DC", "Monocyte")

  x <- sort(unique(x))
  newX <- c()
  
  for (p in priority) {
    pIndices <- grepl(p, x)
    newX <- append(newX, x[pIndices])
    x <- x[!pIndices]
  }
  
  if (length(newX) > 0) {
    newX <- append(newX, x)
  }
  
  
  return(newX)
}

manualClusterOrder <- unique(seuMerged$manualAnnot)
manualClusterOrder <- factor(manualClusterOrder, levels = customSortAnnotation(manualClusterOrder))
#### END run this again after loading manual annot rdata


################################################################################
# wgcna analysis for pLN only
################################################################################
wcgnaCheckpointFile <- "rds/postNetworkPostModule_pLN_ND_T1D_v2.rds"
wcgnaCheckpointImage <- "rds/postNetworkPostModule_pLN_ND_T1D_v2.RData"
if (!file.exists(wcgnaCheckpointFile)) {
  # filter dataset to remove cells with >95 percentile of hsp gene signatures
  seuMergedSmol <- subset(seuMerged, subset = heatShockProgram1 < hspCutoff & TissueCondensed == "pLN" & Disease_Status %in% c("ND", "T1D"))
  gc() # need to clear from RAM after the previous step
  
  seuMergedSmol <- SetupForWGCNA(
    seuMergedSmol,
    gene_select = "fraction", # the gene selection approach
    fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
    wgcna_name = "hpap_citeseq_diseaseSubset_pLN" # the name of the hdWGCNA experiment
  )
  
  seuMergedSmol <- MetacellsByGroups(
    seurat_obj = seuMergedSmol,
    group.by = c("DonorID", "Disease_Status", "phenotypeCondensed"), # specify the columns in seurat_obj@meta.data to group by
    reduction = 'harmonyRNA', # select the dimensionality reduction to perform KNN on
    assay = "RNA",
    k = 25, # nearest-neighbors parameter
    max_shared = 10, # maximum number of shared cells between two metacells
    ident.group = 'Disease_Status', # set the Idents of the metacell seurat object,
    max_iter = 2500,
    verbose = TRUE
  )
  
  seuMergedSmol <- NormalizeMetacells(seuMergedSmol)
  
  seuMergedSmol <- SetDatExpr(
    seuMergedSmol,
    group_name = c("ND", "T1D"),
    group.by = 'Disease_Status',
    assay = "RNA"
  )
  
  # Test different soft powers:
  seuMergedSmol <- TestSoftPowers(
    seuMergedSmol,
    networkType = 'signed'
  )
  
  # plot the results:
  plot_list <- PlotSoftPowers(seuMergedSmol)
  wrap_plots(plot_list, ncol = 2)
  
  # construct co-expression network:
  seuMergedSmol <- ConstructNetwork(
    seuMergedSmol,
    soft_power = 9, #adjust this based on the threshold identified from previous graph...
    setDatExpr = FALSE,
    tom_name = 'nd_t1d_pln_topological' # name of the topological overlap matrix written to disk
  )
  
  PlotDendrogram(seuMergedSmol, main='ND/T1D pLN hdWGCNA Dendrogram')
  gc()
  
  seuMergedSmol <- ScaleData(seuMergedSmol, features = VariableFeatures(seuMergedSmol), assay = "RNA")
  save.image("rds/postNetworkPreModule_pLN_ND_T1D_v2.RData")
  
  # this step is massive in RAM usage...
  seuMergedSmol <- ModuleEigengenes(
    seuMergedSmol,
    group.by.vars = c("DonorID")
  )
  
  save.image(wcgnaCheckpointImage)
  
  seuMergedSmol <- ModuleConnectivity(
    seuMergedSmol,
    group.by = 'Disease_Status', group_name = c("ND", "T1D")
  )
  
  seuMergedSmol <- ResetModuleNames(
    seuMergedSmol,
    new_name = "T1D-M"
  )
  
  PlotKMEs(seuMergedSmol, ncol = 7, text_size = 3.5)
  
  modulesInterest <- GetModules(seuMergedSmol)
  hubInterest <- GetHubGenes(seuMergedSmol, n_hubs = 25)
  head(hubInterest)
  
  seuMergedSmol <- ModuleExprScore(
    seuMergedSmol,
    n_genes = 25,
    method = 'Seurat'
  )
  
  # get hMEs from seurat object
  MEs <- GetMEs(seuMergedSmol, harmonized = TRUE)
  mods <- colnames(MEs)
  mods <- mods[mods != 'grey']
  
  # add hMEs to Seurat meta-data:
  seuMergedSmol@meta.data <- cbind(seuMergedSmol@meta.data, MEs)
  seuMergedSmol$Disease_Status <- factor(seuMergedSmol$Disease_Status, levels = c("ND", "T1D"))
  
  saveRDS(seuMergedSmol, file = wcgnaCheckpointFile)
} else {
  load(wcgnaCheckpointImage)
  seuMergedSmol <- readRDS(wcgnaCheckpointFile)
}

# compute correlations to traits
seuMergedSmol <- ModuleTraitCorrelation(
  seuMergedSmol,
  traits = c("Disease_Status"),
  group.by = "manualAnnot"
)

moduleTraitCorRes <- GetModuleTraitCorrelation(seuMergedSmol)

corMatrix <- do.call("rbind", moduleTraitCorRes$cor)
fdrMatrix <- do.call("rbind", moduleTraitCorRes$fdr)

tmp <- Heatmap(
  matrix = corMatrix,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_row_names = TRUE,
  row_title_rot = 0,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  row_gap = unit(5, "pt"),
  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
    fdrVal <- fdrMatrix[i,j]
    fdrText <- symnum(
      fdrVal,
      corr = FALSE, na = FALSE,
      cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
      symbols = c("***", "**", "*", " ", ""))

    grid.text(fdrText, x, y, gp = gpar(fontsize = 8))
  },
  name = "Module Trait Correlation",
  heatmap_legend_param = list(
    direction = "vertical"),
  row_title_gp = gpar(fontsize = 10)
)
draw(tmp, padding = unit(c(1, 3, 1, 1), "lines"))

tmp <- which(fdrMatrix < 0.05 & abs(corMatrix) > 0.3)
arrayInd(tmp, .dim = dim(corMatrix))


# modules of interest from Greg
gregModulesInterest <- read.csv("modules_genes_followup.csv")

seuMergedHspCutoff <- subset(seuMerged, subset = heatShockProgram1 < hspCutoff & TissueCondensed == "pLN")
seuMergedHspCutoff$Disease_Status <- factor(seuMergedHspCutoff$Disease_Status, levels = c("ND", "AAb+", "T1D"))
seuMergedHspCutoff$manualAnnot <- factor(seuMergedHspCutoff$manualAnnot, levels = levels(manualClusterOrder))

# look at all pop comparison
genesAllPopComparison <- gregModulesInterest[gregModulesInterest$compare_across_all_cell_pops, "gene_name"]
genesAllPopComparison <- unique(genesAllPopComparison)

genesAllPopComparisonAvgMat <- AverageExpression(
  object = seuMergedHspCutoff,
  assays = "RNA",
  features = genesAllPopComparison,
  group.by = c("manualAnnot", "Disease_Status"),
  slot = "data")

genesAllPopComparisonAvgMat <- genesAllPopComparisonAvgMat$RNA

colGroup <- rep(levels(seuMergedHspCutoff$Disease_Status), times = length(levels(seuMergedHspCutoff$manualAnnot)))
colGroup <- factor(colGroup, levels = levels(seuMergedHspCutoff$Disease_Status))

colSplit <- rep(levels(seuMergedHspCutoff$manualAnnot), each = length(levels(seuMergedHspCutoff$Disease_Status)))

p <- Heatmap(
  matrix = t(scale(t(genesAllPopComparisonAvgMat))),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 12),
  column_names_gp = gpar(fontsize = 7),
  column_names_rot = 45,
  column_split = colSplit,
  bottom_annotation = columnAnnotation(
    `Disease\nStatus` = colGroup,
    col = list(`Disease\nStatus` = DISEASESTATUSCOLORS),
    annotation_legend_param = list(`Disease\nStatus` = list(direction = "horizontal"))),
  name = "Scaled Average\nExpression",
  heatmap_legend_param = list(
    direction = "horizontal",
    title_position = "topcenter"),
  row_title_gp = gpar(fontsize = 8),
  column_title_side = "bottom",
  column_title_rot = 90,
  column_title_gp = gpar(fontsize = 10),
)

draw(p,
  merge_legend = TRUE,
  padding = unit(c(5, 3, 3, 1), "lines"),
  annotation_legend_side = "bottom",
  heatmap_legend_side = "bottom")


# look at specific modules of interest
tmp <- gregModulesInterest %>%
  filter(cell_population != "all") %>%
  group_by(module) %>%
  summarize(popInterest = unique(cell_population))

uniqModules <- unique(tmp$module)
modulePopsDict <- lapply(uniqModules, function(x) {
  y <- unlist(tmp[tmp$module == x, "popInterest"])
  names(y) <- NULL
  
  return(y)
})
names(modulePopsDict) <- uniqModules

moduleGeneDict <- lapply(uniqModules, function(x) {
  y <- unlist(gregModulesInterest[gregModulesInterest$module == x, "gene_name"])
  y <- unique(y)
  names(y) <- NULL
  
  return(y)
})
names(moduleGeneDict) <- uniqModules

moduleVlnPlots <- lapply(uniqModules, function(x) {
  popsOfInterest <- modulePopsDict[[x]]
  genesOfInterest <- moduleGeneDict[[x]]
  
  if (length(popsOfInterest) > 5) {
    ncols <- 2
  } else {
    ncols <- 4
  }
  
  seu <- subset(seuMergedHspCutoff, subset = manualAnnot %in% popsOfInterest)
  p <- VlnPlot(seu,
    features = genesOfInterest,
    group.by = "manualAnnot",
    split.by = "Disease_Status",
    pt.size = 0,
    ncol = ncols
  ) &
    scale_fill_manual(values = DISEASESTATUSCOLORS)
})

################################################################################
# generating subplots for figures
################################################################################

# post-hsp subset
seuMergedPostHSP <- subset(seuMerged,
  subset = heatShockProgram1 < hspCutoff)
gc()

saveRDS(seuMergedPostHSP,
  paste0("outs/rds/seuMergedPostHSP_forFigures_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".rds"))

print('done')
