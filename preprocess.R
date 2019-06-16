library(affy)
library(frma)

#NEED TO ADD IN A QC CHECK AND ERROR IF GIVEN BAD ARRAY

# get_file_name <- function(direc) {
#   out <- strsplit(direc,'\\\\')
#   folder_depth <- length(out[[1]])
#   file_name <- out[[1]][folder_depth]
# } # SEE IF THIS BIT IS ACTUALLY NEEDED NOW

pre_proc <- function(direcs,file_names) {
  batch <- ReadAffy(filenames = direcs)
  eset <- exprs(frma(batch))
  # file_names <- lapply(direcs,get_file_name)
  colnames(eset) <- file_names
  return(eset)
}



