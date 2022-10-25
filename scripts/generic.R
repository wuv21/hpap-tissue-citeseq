rdsCheckAndRun <- function(fn, f, ...) {
  if (!file.exists(fn)) {
    print(paste0("File (", fn, ") does not exist. Generating now."))
    
    res <- do.call(f, list(...))
    
    print(paste0("Saving file (", fn, ")..."))
    saveRDS(res, fn)
    
  } else {
    print(paste0("File (", fn, ") exists. Loading now."))
    res <- readRDS(fn)
  }
  
  return(res)
}