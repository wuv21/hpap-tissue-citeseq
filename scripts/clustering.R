source("scripts/generic.R")

clusterMerger <- function(
  seuratObj, 
  clusterName,
  newClusterName = paste0("merged", "_", clusterName),
  nnSlot = "weighted.nn",
  minClusterSize = 100
) {
  stopifnot(sum(rownames(seuratObj[[clusterName]]) == seuratObj[[nnSlot]]@cell.names) == nrow(seuratObj[[clusterName]]))
  
  origClusters <- as.character(seuratObj[[clusterName]][, 1])
  
  resTabulated <- table(origClusters)
  idxOfInterest <- origClusters %in% names(resTabulated[resTabulated < minClusterSize])
  
  nnIdx <- seuratObj[[nnSlot]]@nn.idx[idxOfInterest, ]
  nnDist <- seuratObj[[nnSlot]]@nn.dist[idxOfInterest, ]
  
  nnIdxCluster <- apply(nnIdx, MARGIN = 2, function(x) {
    return(origClusters[x])
  })
  
  mergedCluster <- sapply(c(1: nrow(nnIdxCluster)), function(i) {
    x <- nnIdxCluster[i, ]
    dists <- nnDist[i, ]
    
    counts <- table(x)
    
    avgDistByCluster <- aggregate(dists, by = list(x), FUN = mean)
    avgDistByCluster2 <- avgDistByCluster$x
    names(avgDistByCluster2) <- avgDistByCluster$Group.1
    
    countMax <- counts[which(counts == max(counts))]
    if (length(countMax) != 1) {
      avgDistByCluster2 <- avgDistByCluster2[names(avgDistByCluster2) %in% names(countMax)]
      avgDistByCluster2 <- base::sort(avgDistByCluster2, decreasing = FALSE)
    
      message(paste0(seuratObj[[nnSlot]]@cell.names[i], " (index ", i, ") has >1 top NN clusters"))
      message(paste(names(avgDistByCluster2), avgDistByCluster2, sep = ":", collapse = " vs "))
      
      countMax <- avgDistByCluster2[1]
    }
    
    return(names(countMax))
  })
  
  seuratObj[[newClusterName]] <- origClusters
  seuratObj[[newClusterName]][idxOfInterest, 1] <- mergedCluster
  
  return(seuratObj)
}


getClusterDifferences <- function(
  seu,
  clusterName,
  tsvFn,
  tsaCatalog = NULL,
  assay = "adt",
  featuresToUse = NULL,
  parallelize = FALSE,
  nCores = 1,
  findMarkerMethod = "wilcox",
  ...
) {
  
  if (!is.null(tsaCatalog) & !is.null(featuresToUse)) {
    featuresToUse <- rownames(tsaCatalog[!tsaCatalog$isCtrl, "DNA_ID"])
  }
  
  Idents(seu) <- clusterName
  
  if (parallelize) {
    message("Parallelization enabled. Running FindMarkers in parallel.")
    
    # derived from https://gist.github.com/diazdc/1735102c243cd16acb1b1f3fd09a26e1
    
    startTime <- Sys.time()
    mcFindMarkers <- function(i) {
      ident1 <- i
      
      message(paste0("Comparing ", ident1, " and other cells\n"))
      res <- FindMarkers(
        seu,
        assay = assay,
        ident.1 = as.character(ident1),
        test.use = findMarkerMethod,
        features = featuresToUse,
        ...)
      
      endTime <- Sys.time()
      timeTaken <- round(endTime - startTime, digits = 1)
      
      message(paste0("Finished in ", timeTaken, " seconds."))
      
      print(head(res))
      if (nrow(res) != 0) {
        res$gene <- rownames(res)
        res$cluster <- rep(as.character(i), nrow(res))
        
        return(res)
      } else {
        return(NULL)
      }
    }
    
    uniqIdents <- unique(Idents(seu))
    # markers <- parallel::mclapply(uniqIdents, mcFindMarkers, mc.cores = nCores)
    markers <- lapply(unique(Idents(seu)), mcFindMarkers)
    names(markers) <- uniqIdents
    
    clusterMarkers <- dplyr::bind_rows(markers)
    
  } else {
    message("Parallelization disabled Running FindAllMarkers in serial mode.\n")
    
    clusterMarkers <- FindAllMarkers(
      seu,
      assay = assay,
      slot = "data",
      test.use = findMarkerMethod,
      features = featuresToUse,
      ...)
  }

  
  if (!("gene" %in% names(clusterMarkers))) {
    clusterMarkers$gene <- rownames(clusterMarkers)
  }
  
  if (!is.null(tsaCatalog)) {
    clusterMarkers <- clusterMarkers %>%
      left_join(tsaCatalog %>% dplyr::select(DNA_ID, cleanName), by = c("gene" = "DNA_ID"))
  }
  
  clusterMarkersFilt <- clusterMarkers %>%
    filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    filter(avg_log2FC > 0) %>%
    mutate(piScore = avg_log2FC * -log10(p_val_adj)) %>%
    dplyr::arrange(desc(piScore), .by_group = TRUE)
  
  write.table(clusterMarkersFilt, file = tsvFn, row.names = FALSE, quote = FALSE, sep = "\t")
}