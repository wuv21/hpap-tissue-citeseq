source("scripts/generic.R")
source("scripts/adtProcessing.R")


.processFromCellranger <- function(
  dataDir,
  sampleID,
  runN,
  plotDir,
  minCells = 3,
  minFeatures = 200
) {
  
  seuTmp <- Read10X(data.dir = dataDir, strip.suffix = TRUE)
  seuTmp <- CreateSeuratObject(
    counts = seuTmp,
    min.cells = minCells,
    min.features = minFeatures)
  
  seuTmp <- RenameCells(seuTmp, add.cell.id = sampleID)
  seuTmp <- AddMetaData(seuTmp, runN, col.name = "runN")
  seuTmp <- AddMetaData(seuTmp, sampleID, col.name = "well")
  
  seuTmp[["percent.mt"]] <- PercentageFeatureSet(seuTmp, pattern = "^MT-")
  ggsave(
    VlnPlot(seuTmp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3),
    filename = paste0(plotDir, "/", sampleID, "_seuOrigVlnPlts.png"),
    height = 5,
    width = 8)
  
  return(seuTmp)
}

.processHTOFromKallisto <- function(
  smpl,
  seuTmp,
  adt,
  htoOutDir
) {
  
  hto <- createMatFromKallisto(
    sample = smpl,
    dir = htoOutDir,
    dataMode = "hto",
    sampleSuffix = ""
  )
  
  colnames(hto) <- paste0(smpl, "_", colnames(hto))
  
  seuTmp <- subset(seuTmp, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 12.5)
  allowableCbcs <- Reduce(intersect, list(colnames(seuTmp), colnames(adt), colnames(hto)))
  seuTmp <- subset(seuTmp, cells = allowableCbcs)
  
  seuTmp <- NormalizeData(seuTmp)
  
  seuTmp[["hto"]] <- CreateAssayObject(counts = hto[, allowableCbcs])
  seuTmp[["adt"]] <- CreateAssayObject(counts = adt[, allowableCbcs])
  
  seuTmp <- NormalizeData(seuTmp, assay = "hto", normalization.method = "CLR")
  seuTmp <- HTODemux(seuTmp, assay = "hto", positive.quantile = 0.99)
  
  return(seuTmp)
}


initialSeuratListMaker <- function(
  version,
  sampleIds,
  sampleRuns,
  rnaOutDir,
  adtOutDir,
  htoOutDir,
  plotDir = "outs/qc"
) {
  
  seuList <- lapply(seq_along(sampleIds), function(i) {
    x <- sampleIds[i]
    runN <- sampleRuns[i]
    
    print(paste0("working on: ", x))

    seuFileRdsFn <- paste0("rds/seuOrig", version, "/", x, "_seuOrig.rds")
    seuDemuxRdsFn <- paste0("rds/seuDemux", version, "/", x, "_seuDemux.rds")
    
    if (!file.exists(seuDemuxRdsFn)) {
      seu <- rdsCheckAndRun(
        fn = seuFileRdsFn,
        f = .processFromCellranger,
        dataDir = paste0(rnaOutDir, "/", x),
        sampleID = x,
        runN = runN,
        plotDir = plotDir)
      
      # adt creation
      adt <- createMatFromKallisto(
        sample = x,
        dir = adtOutDir,
        dataMode = "adt",
        outRDSPath = paste0("rds/ADT/", x, "_ADT.rds"),
        sampleSuffix = "")
      
      colnames(adt) <- paste0(x, "_", colnames(adt))
      
      seu <- rdsCheckAndRun(
        fn = seuDemuxRdsFn,
        f = .processHTOFromKallisto,
        smpl = x,
        seuTmp = seu,
        adt = adt,
        htoOutDir = htoOutDir)
      
      return(seuDemuxRdsFn)
    } else {
      # seu <- readRDS(seuDemuxRdsFn)
      return(seuDemuxRdsFn)
    }
  })
  
  return(seuList)  
}

  