library(affy)
library(frma)

# NEED TO ADD IN A QC CHECK AND ERROR IF GIVEN BAD ARRAY

pre_proc <- function(direcs, file_names, crossProbes) {
  batch <- ReadAffy(filenames = direcs)
  eset <- exprs(frma(batch))
  colnames(eset) <- file_names
  if (!is.null(crossProbes)) {
    hyb <- read.csv(crossProbes, stringsAsFactors = FALSE, header = FALSE)
    hyb <- as.list(hyb[, 1])
    eset <- eset[!(rownames(eset) %in% hyb), , drop = FALSE]
  }
  return(eset)
}
