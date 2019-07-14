library(affy)
library(frma)

#' Preprocess
#'
#' Applies fRMA preprocessing to input cel files
#'
#' @param direcs the directories or datapaths of cel files
#' @param file_names the names of the cel files
#' @param cross_probes the datapath for csv file of probes to remove
#' @return dataframe of fRMA preprocessed data (all samples)
#'
Preprocess <- function(direcs, file_names, cross_probes) {
  batch <- ReadAffy(filenames = direcs)
  eset <- exprs(frma(batch))
  colnames(eset) <- file_names
  if (!is.null(cross_probes)) {
    hyb <- read.csv(cross_probes, stringsAsFactors = FALSE, header = FALSE)
    hyb <- as.list(hyb[, 1])
    eset <- eset[!(rownames(eset) %in% hyb), , drop = FALSE]
  }
  return(eset)
}
