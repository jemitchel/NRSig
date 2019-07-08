
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


# res <- pipeline("C:/Users/jonat/Documents/Research/app test folder/ER.CEL",NULL)
# saveRDS(res,"C:/Users/jonat/Documents/R/NRSig-app/data/exres2.rds")

