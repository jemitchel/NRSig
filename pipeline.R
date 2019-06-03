
pipeline <- function(files) {
  source("C:/Users/jonat/Documents/R/NRSig-app/preprocess.R")
  samples_matrix <- pre_proc(files)
  
  # Computes enrichment
  source("C:/Users/jonat/Documents/R/NRSig-app/compute_enriched_NRs.R")
  results <- CalcEnrich(samples_matrix)
  return(results[1])
}

pipeline("C:\\Users\\jonat\\Documents\\Research\\app test folder\\c2_(HG-U133_Plus_2).CEL")
