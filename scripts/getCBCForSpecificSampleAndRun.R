sample <- "2_3"
run <- "20220531_hpap_2_Pool_Run1_Well5"

# sample <- "4_4"
# run <- "20220801_hpap_4_Pool_Well2"
# 
# sample <- "5_4"
# run <- "20220923_hpap_5_Pool_Well2"

tmpSeuMerged <- readRDS("rds/seuMergedAndHarmonized.rds")
tmpSeuMerged$combinedMeta <- paste(tmpSeuMerged$runN, tmpSeuMerged$HTO_antibody, sep = "_")

cbcOfInterest <- names(tmpSeuMerged$combinedMeta[tmpSeuMerged$combinedMeta == sample])
cbcOfInterest <- cbcOfInterest[grepl(run, cbcOfInterest)]

cbcOfInterest <- stringr::str_match(cbcOfInterest, "_([ATGC]+$)")[, 2]
write.table(
  x = paste0("CB:Z:", cbcOfInterest, "-1"),
  file = paste0(run, "___", sample, ".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
