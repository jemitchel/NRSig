
source("C:/Users/jonat/Documents/R/NRSig-app/preprocess.R")
source("C:/Users/jonat/Documents/R/NRSig-app/compute_enriched_NRs.R")

pipeline <- function(files,rem) {
  # preprocesses cel files
  fBaseNames <- unname(sapply(files,basename))
  samples_matrix <- pre_proc(files,fBaseNames,rem)

  # computes enrichment results
  results <- CalcEnrich(samples_matrix,rem)
  return(results)
}


# pipeline("C:/Users/jonat/Documents/Research/app test folder/c2_(HG-U133_Plus_2).CEL",NULL)

