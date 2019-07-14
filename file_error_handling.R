#' celFileError
#'
#' Determines if there is a problem with uploaded cel files
#'
#' @param flist a list of input files datapath
#' @return The resulting error or NULL if no error
#'
celFileError <- function(flist) {
  if (is.null(flist)) {
    return("Error: Please upload files to process")
  } else {
    for (i in 1:length(flist)) {
      len <- nchar(flist[[i]])
      last4 <- substr(flist[[i]], len - 3, len)
      if (!identical(last4, ".CEL")) {
        return("Error: Upload includes file(s) not ending in .CEL ")
      }
      cdf_name <- whatcdf(flist[[i]])
      if (!identical(cdf_name, "HG-U133_Plus_2")) {
        return("Error: Input file(s) are for wrong platform. Use CEL files
               for HG-U133_Plus_2 only.")
      }
      }
    }
  return(NULL)
  }

#' csvFileError
#'
#' Determines if there is a problem with uploaded csv files
#'
#' @param fl a file datapath for list of crosshybridized probes
#' @return The resulting error or NULL if no error
#'
csvFileError <- function(fl) {
  if (is.null(fl)) {
    return(NULL)
  } else {
    len <- nchar(fl)
    last4 <- substr(fl, len - 3, len)
    if (!identical(last4, ".csv")) {
      return("Error: Uploaded probe list not a .csv")
    }
    dimension <- ncol(read.csv(fl))
    if (dimension > 1) {
      return("Error: Uploaded probe list has more than 1 column. Ensure list
             is a single column of Affymetrix probe symbols.")
    }
    }
  return(NULL)
  }